#' Aggregate observed data based on a tree
#'
#' Aggregate observed values based on a column (sample) tree, e.g. for
#' differential state analysis. The returned object will contain one
#' abundance matrix for each node in the tree, with columns corresponding
#' to sample IDs and rows corresponding to the same features as the rows of
#' the input object. The nodes correspond to
#' either the original sample clusters, or larger metaclusters encompassing
#' several of the original clusters (defined by the provided column tree).
#'
#' @author Ruizhu Huang, Charlotte Soneson
#' @export
#'
#' @param TSE A \code{TreeSummarizedExperiment} object. Rows represent
#'     variables (e.g., genes) and columns represent observations (e.g., cells).
#'     The object must contain a column tree, whose tips represent initial
#'     cell clusters (the \code{cluster_id} annotation indicates which of
#'     these clusters a cell belongs to). The internal nodes represent
#'     increasingly coarse partitions of the cells obtained by successively
#'     merging the original clusters according to the provided column tree.
#' @param assay The name or index of the assay from \code{TSE} to aggregate
#'     values from.
#' @param sample_id A character scalar indicating the column of
#'     \code{colData(TSE)} that corresponds to the sample ID. These will be the
#'     columns of the output object.
#' @param group_id A character scalar indicating the column of
#'     \code{colData(TSE)} that corresponds to the group/condition. This
#'     information will be propagated to the aggregated object.
#' @param cluster_id A character scalar indicating the column of
#'     \code{colData(TSE)} that corresponds to the initial cluster ID for
#'     each cell.
#' @param FUN The aggregation function.
#' @param message A logical scalar, indicating whether progress messages
#'     should be printed to the console.
#'
#' @returns A \code{SummarizedExperiment} object. Each assay represents the
#'     aggregated values for one node in the tree.
#'
#' @importFrom SummarizedExperiment SummarizedExperiment colData assays
#'     assayNames
#' @importFrom TreeSummarizedExperiment findDescendant convertNode showNode
#'     colTree
#' @importFrom dplyr select distinct mutate
#' @importFrom utils flush.console
#'
#' @examples
#' suppressPackageStartupMessages({
#'     library(TreeSummarizedExperiment)
#'     library(ape)
#'     library(ggtree)
#' })
#'
#' set.seed(1L)
#' tr <- rtree(3, tip.label = LETTERS[seq_len(3)])
#' ggtree(tr) +
#'     geom_text(aes(label = node), hjust = -1, vjust = 1) +
#'     geom_text(aes(label = label), hjust = -1, vjust = -1)
#'
#' cc <- matrix(rpois(60, 10), nrow = 6)
#' rownames(cc) <- paste0("gene", seq_len(6))
#' colnames(cc) <- paste0("cell", seq_len(10))
#' cd <- data.frame(sid = rep(seq_len(2), each = 5),
#'                  gid = rep(letters[seq_len(2)], each = 5),
#'                  cid = sample(LETTERS[seq_len(3)], size = 10,
#'                               replace = TRUE),
#'                  stringsAsFactors = FALSE)
#' tse <- TreeSummarizedExperiment(assays = list(counts = cc),
#'                                 colTree = tr,
#'                                 colNodeLab = cd$cid,
#'                                 colData = cd)
#'
#' out <- aggDS(TSE = tse, assay = "counts", sample_id = "sid",
#'              group_id = "gid", cluster_id = "cid")
#'
#' ## Aggregated counts for the node 5
#' SummarizedExperiment::assay(out, "alias_5")
#' ## This is equal to the sum of the counts from nodes 1 and 2
#' SummarizedExperiment::assay(out, "alias_1")
#' SummarizedExperiment::assay(out, "alias_2")
#'
aggDS <- function(TSE, assay = "counts", sample_id = "sample_id",
                  group_id = "group_id", cluster_id = "cluster_id",
                  FUN = sum, message = FALSE) {

    ## Check arguments
    ## -------------------------------------------------------------------------
    .assertVector(x = TSE, type = "TreeSummarizedExperiment")
    stopifnot(length(assay) == 1)
    if (is.character(assay)) {
        stopifnot(assay %in% SummarizedExperiment::assayNames(TSE))
    } else {
        stopifnot(
            is.numeric(assay) &&
                assay %in% seq_len(length(SummarizedExperiment::assays(TSE))))
    }
    .assertScalar(x = sample_id, type = "character")
    .assertScalar(x = group_id, type = "character")
    .assertScalar(x = cluster_id, type = "character")
    stopifnot(all(c(sample_id, group_id, cluster_id) %in%
                      colnames(SummarizedExperiment::colData(TSE))))
    .assertScalar(x = message, type = "logical")
    .assertVector(x = FUN, type = "function")
    .assertVector(x = TreeSummarizedExperiment::colTree(TSE), type = "phylo")
    stopifnot(all(SummarizedExperiment::colData(TSE)[[cluster_id]] %in%
                      TreeSummarizedExperiment::colTree(TSE)$tip.label))

    ## Cell information
    ## -------------------------------------------------------------------------
    if (message) {
        message("Extracting cell information ...")
    }
    tree <- TreeSummarizedExperiment::colTree(TSE)
    cell_info <- SummarizedExperiment::colData(TSE)[, c(sample_id, group_id,
                                                        cluster_id)]
    colnames(cell_info) <- c("sample_id", "group_id", "cluster_id")
    ## Find the node number for each original cluster and add to cell info
    cell_info$node <- TreeSummarizedExperiment::convertNode(
        tree = tree, node = as.character(cell_info$cluster_id))

    ## Tree information
    ## -------------------------------------------------------------------------
    if (message) {
        message("Extracting tree information ...")
    }
    ## Get node number corresponding to each alias
    node <- TreeSummarizedExperiment::showNode(tree = tree, only.leaf = FALSE,
                                               use.alias = TRUE)
    ## Find descendant leaves (node numbers) for each node
    desd_list <- TreeSummarizedExperiment::findDescendant(
        tree = tree, node = node, only.leaf = TRUE, self.include = TRUE)
    names(desd_list) <- names(node)

    ## Indicators
    ## -------------------------------------------------------------------------
    if (message) {
        message("Grouping cells by samples and nodes ...")
    }
    ## For each node, find the cells belonging to it (included in one of the
    ## corresponding leaf clusters)
    ind_list1 <- lapply(desd_list, FUN = function(x) {
        which(cell_info$node %in% x)
    })

    ## For each sample, find the cells belonging to it
    sp <- unique(cell_info$sample_id)
    ind_list2 <- lapply(sp, FUN = function(x) {
        which(cell_info$sample_id %in% x)
    })
    names(ind_list2) <- sp

    ## For each node, construct a list where each element gives the cells
    ## corresponding to that node and a given sample
    ind_list3 <- lapply(ind_list1, FUN = function(x) {
        xx <- lapply(ind_list2, intersect, y = x)
    })

    if (message) {
        message("Preparing data on each node ...")
    }
    ## Reshape data
    ## -------------------------------------------------------------------------
    asy <- SummarizedExperiment::assays(TSE)[[assay]]
    ## Initialize the list that will hold assays (one for each node)
    dat_list <- vector("list", length(ind_list3))
    names(dat_list) <- names(ind_list3)
    ## Populate each node assay
    for (i in seq_along(dat_list)) {
        ind_i <- ind_list3[[i]]
        sp_i <- names(ind_i)
        asy_i <- lapply(ind_i, FUN = function(x) {
            ## Get assay values for node and sample
            ## Genes in rows, cells in columns
            xx <- asy[, x, drop = FALSE]

            ## Aggregate
            ax <- apply(xx, 1, FUN = FUN)
        })
        ## Assemble aggregated matrix
        asy_i <- do.call(cbind, asy_i)
        colnames(asy_i) <- sp_i
        dat_list[[i]] <- asy_i

        if (message) {
            message(i, " out of ", length(dat_list),
                    " nodes finished", "\r", appendLF = FALSE)
            utils::flush.console()
        }
    }

    ## Sample information
    ## -------------------------------------------------------------------------
    if (message) {
        message("Working on sample information ...")
    }
    sample_df <- SummarizedExperiment::colData(TSE)[, c(sample_id, group_id)] |>
        data.frame() |>
        distinct()
    sample_info <- sample_df[, group_id, drop = FALSE]
    rownames(sample_info) <- sample_df[[sample_id]]
    sample_info <- sample_info[colnames(dat_list[[1]]), , drop = FALSE]

    ## Metadata
    ## -------------------------------------------------------------------------
    if (message) {
        message("Working on metadata ...")
    }
    cell_n <- table(cell_info[["sample_id"]])
    experiment_info <- cell_info |>
        data.frame() |>
        select(sample_id, group_id) |>
        distinct() |>
        mutate(n_cells = cell_n[sample_id])

    agg_pars <- list(assay = assay,
                     by = c(cluster_id, sample_id),
                     fun = deparse(substitute(FUN)),
                     scale = FALSE)

    ## On leaf level
    cell_tab <- table(cell_info[["cluster_id"]],
                      cell_info[["sample_id"]])
    rn <- TreeSummarizedExperiment::convertNode(tree = tree,
                                                node = rownames(cell_tab))

    ## On all levels
    ## 'node' has the node number for each alias
    desd <- TreeSummarizedExperiment::findDescendant(
        tree = tree, node = node, only.leaf = TRUE, self.include = TRUE)
    cell_tab <- as.matrix(cell_tab)
    clus_tab <- lapply(desd, FUN = function(x) {
        ii <- match(x, rn)
        xi <- cell_tab[ii, , drop = FALSE]
        apply(xi, 2, sum)
    })
    clus_tab <- do.call(rbind, clus_tab)
    rownames(clus_tab) <- names(node)
    clus_tab <- clus_tab[, colnames(dat_list[[1]])]

    meta <- list(experiment_info = experiment_info,
                 agg_pars = agg_pars,
                 n_cells = clus_tab)

    ## Output
    if (message) {
        message("Output data ...")
    }
    out <- SummarizedExperiment::SummarizedExperiment(assays = dat_list,
                                                      colData = sample_info,
                                                      metadata = meta)
    return(out)
}
