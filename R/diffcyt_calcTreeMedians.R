#' @export
#' @rdname diffcyt_workflow
#'
#' @importFrom TreeSummarizedExperiment TreeSummarizedExperiment convertNode
#'     aggTSE rowTree
#' @importFrom dplyr select distinct
#' @importFrom methods is
#' @importFrom stats median
#' @importFrom SummarizedExperiment rowData assays colData
#' @importFrom S4Vectors metadata
#'
calcTreeMedians <- function(d_se, tree, message = FALSE) {

    ## Check arguments
    ## -------------------------------------------------------------------------
    .assertScalar(x = message, type = "logical")
    .assertVector(x = d_se, type = "SummarizedExperiment")
    if (!("cluster_id" %in% (colnames(SummarizedExperiment::rowData(d_se))))) {
        stop("Data object does not contain cluster labels.
             Run 'diffcyt::generateClusters' to generate cluster labels.")
    }
    if (!methods::is(tree, "phylo")) {
        stop("tree is not a phylo object.
             Run 'buildTree(d_se)' to generate the tree.")
    }

    ## Create the tse object
    ## -------------------------------------------------------------------------
    rlab <- as.character(SummarizedExperiment::rowData(d_se)$cluster_id)
    d_lse <- TreeSummarizedExperiment::TreeSummarizedExperiment(
        assays = SummarizedExperiment::assays(d_se),
        rowData = SummarizedExperiment::rowData(d_se),
        rowTree = tree,
        rowNodeLab = rlab,
        colData = SummarizedExperiment::colData(d_se),
        metadata = S4Vectors::metadata(d_se))

    ## Retain only type or state markers
    d_lse <- d_lse[, SummarizedExperiment::colData(d_lse)$marker_class %in%
                       c("type", "state")]
    if (ncol(d_lse) == 0) {
        stop("No type or state markers found in the object.")
    }

    ## Calculate required summaries
    ## -------------------------------------------------------------------------
    rd <- SummarizedExperiment::rowData(d_lse)
    sid <- unique(rd$sample_id)
    nodes <- sort(unique(as.vector(tree$edge)))
    labs <- TreeSummarizedExperiment::convertNode(
        tree = tree, node = nodes, use.alias = TRUE)

    lc <- lapply(seq_along(sid), FUN = function(i) {
        if (message) {
            message("Working on ", i, " out of ", length(sid), " samples.")
        }
        x <- sid[i]
        sel <- rd$sample_id == x
        xx <- d_lse[sel, ]
        ax <- TreeSummarizedExperiment::aggTSE(
            x = xx, rowLevel = nodes,
            rowFun = function(x) {
                stats::median(x, na.rm = TRUE)
            }, message = message)

        ## Counts
        cx <- SummarizedExperiment::assays(ax)[[1]]

        ## For missing nodes
        nam <- setdiff(labs, rownames(cx))
        if (length(nam)) {
            mx <- matrix(NA, nrow = length(nam), ncol = ncol(cx))
            rownames(mx) <- nam
            colnames(mx) <- colnames(cx)
        } else {
            mx <- NULL
        }

        cx <- rbind(cx, mx)[labs, ]
        return(cx)
    })

    ## assays: reshape data to split by markers
    rc <- do.call(rbind, lc)
    lcm <- lapply(seq_len(ncol(rc)), FUN = function(x) {
        xx <- matrix(rc[, x], ncol = length(sid), byrow = FALSE)
        colnames(xx) <- sid
        rownames(xx) <- labs
        return(xx)
    })
    names(lcm) <- colnames(d_lse)

    ## Row data
    rd <- cbind.data.frame(cluster_id = factor(labs, levels = labs))
    rownames(rd) <- labs

    ## Column data
    cd <- SummarizedExperiment::rowData(d_se) |>
        data.frame() |>
        select(group_id, sample_id) |>
        distinct()
    cd <- cd[match(cd$sample_id, sid), , drop = FALSE]

    ## Metadata
    tM <- SummarizedExperiment::colData(d_lse)$marker_name[
        SummarizedExperiment::colData(d_lse)$marker_class == "type"]
    sM <- SummarizedExperiment::colData(d_lse)$marker_name[
        SummarizedExperiment::colData(d_lse)$marker_class == "state"]
    md <- list(id_type_markers = colnames(d_lse) %in% tM,
               id_state_markers = colnames(d_lse) %in% sM)

    ## Row node label
    rnl <- TreeSummarizedExperiment::convertNode(
        tree = TreeSummarizedExperiment::rowTree(d_lse), node = nodes,
        use.alias = TRUE)

    ## Assemble output
    ## -------------------------------------------------------------------------
    out <- TreeSummarizedExperiment::TreeSummarizedExperiment(
        assays = lcm, rowData  = rd, colData = cd, metadata = md,
        rowTree = TreeSummarizedExperiment::rowTree(d_lse),
        rowNodeLab = rnl)

    return(out)
}
