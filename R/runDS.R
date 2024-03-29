#' Test for differential state using edgeR
#'
#' Test for differential state of entities using functions from the
#' \code{\link{edgeR}} package. This adapts \code{\link{edgerWrp}} to accept
#' input as a
#' \code{\link[SummarizedExperiment]{SummarizedExperiment}} (SE) object
#' instead of \code{matrix}. Each \code{assay} should correspond to data
#' for one node of the tree. Samples are in columns and
#' features are in rows. The sample information is
#' in \code{colData}. The tree that stores the hierarchical relation between
#' the \code{assays} is provided via the argument \code{tree}.
#'
#' @author Ruizhu Huang
#' @export
#'
#' @param SE A \code{SummarizedExperiment} object, typically generated by
#'     \code{aggDS}.
#' @param tree A \code{phylo} object. Each \strong{assay} of \code{SE}
#'     stores data for one node of the tree.
#' @param option Either \code{"glm"} or \code{"glmQL"}. If \code{"glm"},
#'     \code{\link[edgeR]{glmFit}} and \code{\link[edgeR]{glmLRT}} are used;
#'     otherwise, \code{\link[edgeR]{glmQLFit}} and
#'     \code{\link[edgeR]{glmQLFTest}} are used. Details about the difference
#'     between two options are in the help page of
#'     \code{\link[edgeR]{glmQLFit}}.
#' @param design A numeric design matrix. If \code{NULL}, all columns of
#'     the sample annotation will be used to create the design matrix.
#' @param contrast A numeric vector specifying one contrast of
#'     the linear model coefficients to be tested equal to zero. Its length
#'     must equal to the number of columns of design. If \code{NULL}, the last
#'     coefficient will be tested equal to zero.
#' @param filter_min_count A numeric value, passed to \strong{min.count} of
#'     \code{\link[edgeR]{filterByExpr}}.
#' @param filter_min_total_count A numeric value, passed to
#'     \strong{min.total.count} of \code{\link[edgeR]{filterByExpr}}.
#' @param filter_large_n A numeric value, passed to \strong{large.n} of
#'     \code{\link[edgeR]{filterByExpr}}.
#' @param filter_min_prop A numeric value, passed to \strong{min.prop} of
#'     \code{\link[edgeR]{filterByExpr}}.
#' @param normalize A logical scalar indicating whether to estimate
#'     normalization factors (using \code{\link[edgeR]{calcNormFactors}}).
#' @param normalize_method Normalization method to be used. See
#'     \code{\link[edgeR]{calcNormFactors}} for more details.
#' @param min_cells A numeric scalar specifying the minimal number of cells
#'     in a node required to include a node in the analysis. The information
#'     about the number of cells per node and sample should be available in
#'     \code{metadata(SE)$n_cells}. A node is retained if at least half of the
#'     samples have at least \code{min_cells} cells belonging to the node.
#' @param group_column The name of the column in the sample annotation
#'     providing group labels for samples. This annotation is used for
#'     filtering.
#' @param design_terms The names of columns from the sample annotation
#'     that will be used to generate the design matrix. This is ignored if
#'     \strong{design} is provided.
#' @param message A logical scalar, indicating whether progress messages
#'     should be printed.
#' @param ... More arguments to pass to \code{\link[edgeR]{glmFit}}
#'     (\code{option = "glm"} or \code{\link[edgeR]{glmQLFit}}
#'     (\code{option = "glmQL"}).
#'
#' @returns A list with entries \strong{edgeR_results}, \strong{tree}, and
#' \strong{nodes_drop}.
#' \describe{
#'     \item{edgeR_results}{A list. Each of the elements contains the output of
#'          \code{\link[edgeR]{glmQLFTest}} or
#'          \code{\link[edgeR]{glmLRT}} for one node, depending on the specified
#'          \code{option}.}
#'     \item{tree}{The hierarchical structure of entities that was stored in the
#'          input \code{SE}.}
#'     \item{nodes_drop}{A vector storing the alias node labels of entities
#'          that are filtered before analysis due to low counts.}
#' }
#'
#' @importFrom SummarizedExperiment assays colData assayNames assay
#' @importFrom S4Vectors metadata
#' @importFrom utils flush.console
#'
#' @examples
#' suppressPackageStartupMessages({
#'     library(TreeSummarizedExperiment)
#' })
#' ## Load example data
#' ds_tse <- readRDS(system.file("extdata", "ds_sim_20_500_8de.rds",
#'                               package = "treeclimbR"))
#' ds_se <- aggDS(TSE = ds_tse, assay = "counts", sample_id = "sample_id",
#'                group_id = "group", cluster_id = "cluster_id", FUN = sum)
#' ## Information about the number of cells is provided in the metadata
#' S4Vectors::metadata(ds_se)$n_cells
#'
#' ds_res <- runDS(SE = ds_se, tree = colTree(ds_tse), option = "glmQL",
#'                 group_column = "group", contrast = c(0, 1),
#'                 filter_min_count = 0, filter_min_total_count = 1,
#'                 design = model.matrix(~ group, data = colData(ds_se)),
#'                 filter_min_prop = 0, min_cells = 5, message = FALSE)
#' ## Top differential features (across nodes)
#' nodeResult(ds_res, type = "DS")
#'
runDS <- function(SE, tree, option = c("glm", "glmQL"),
                  design = NULL, contrast = NULL,
                  filter_min_count = 1, filter_min_total_count = 15,
                  filter_large_n = 10, filter_min_prop = 1,
                  min_cells = 10, normalize = TRUE,
                  normalize_method = "TMM", group_column = "group_id",
                  design_terms = "group_id",
                  message = TRUE, ...) {

    ## Check arguments
    ## -------------------------------------------------------------------------
    .assertVector(x = SE, type = "SummarizedExperiment")
    .assertVector(x = tree, type = "phylo")
    .checkEdgeRArgs(design = design, contrast = contrast,
                    filter_min_count = filter_min_count,
                    filter_min_total_count = filter_min_total_count,
                    filter_large_n = filter_large_n,
                    filter_min_prop = filter_min_prop, normalize = normalize,
                    normalize_method = normalize_method,
                    design_terms = design_terms)
    .assertScalar(x = min_cells, type = "numeric")
    .assertVector(x = group_column, type = "character")
    .assertScalar(x = message, type = "logical")

    ## Get node aliases
    ## -------------------------------------------------------------------------
    alias <- SummarizedExperiment::assayNames(SE)

    ## Filter nodes by the number of cells
    ## -------------------------------------------------------------------------
    ncell <- S4Vectors::metadata(SE)$n_cells
    ind_cell <- apply(ncell, 1, FUN = function(x) {
        sum(x >= min_cells) >= 0.5 * length(x)
    })
    message(sum(!ind_cell), " nodes are ignored, as they don't contain ",
            "at least ", min_cells, " cells in at least half of the samples.")
    rem_nodes <- alias[!ind_cell]
    alias <- alias[ind_cell]

    ## Initialize result list (one element per retained node)
    ## -------------------------------------------------------------------------
    res <- vector("list", length(alias))
    names(res) <- alias
    nodes_list <- res

    ## Go through each node and run DS analysis
    ## -------------------------------------------------------------------------
    for (i in seq_along(alias)) {
        if (message) {
            message(i, " out of ", length(alias),
                    " nodes finished", "\r", appendLF = FALSE)
            utils::flush.console()
        }

        ## Extract count matrix, sample information, node information
        count <- SummarizedExperiment::assays(SE)[[alias[i]]]
        sp_info <- SummarizedExperiment::colData(SE)

        ## Check whether sufficient data is available
        sp_keep <- colSums(count, na.rm = TRUE) > 0
        if (sum(!sp_keep)) {
            ## Remove columns without counts
            sp_info <- sp_info[sp_keep, , drop = FALSE]
        }
        ## Get information from the group column
        sp_gr <- unique(sp_info[, group_column])

        if (length(sp_gr) == 1) {
            ## Only one group is available -> return NULL
            res[[alias[i]]] <- NULL
            nodes_list[[alias[i]]] <- alias[i]
        } else {
            ## More than one group is available -> run analysis
            res[[alias[i]]] <- .DS(SE = SE,
                                   assay = alias[i],
                                   option = option,
                                   design = design, contrast = contrast,
                                   filter_min_count = filter_min_count,
                                   filter_min_total_count =
                                       filter_min_total_count,
                                   filter_large_n = filter_large_n,
                                   filter_min_prop = filter_min_prop,
                                   normalize = normalize,
                                   normalize_method = normalize_method,
                                   group_column = group_column,
                                   design_terms = design_terms, ...)
        }
    }
    out <- list(edgeR_results = res,
                tree = tree,
                nodes_drop = c(rem_nodes, unlist(nodes_list)))
    return(out)
}

#' @keywords internal
#' @noRd
#' @importFrom stats as.formula model.matrix
#' @importFrom SummarizedExperiment assay colData
#' @importFrom edgeR filterByExpr
#'
.DS <- function(SE, feature_on_row = TRUE, assay,
                option = c("glm", "glmQL"), design = NULL, contrast = NULL,
                filter_min_count = 10, filter_min_total_count = 15,
                filter_large_n = 10, filter_min_prop = 0.7,
                normalize = TRUE, normalize_method = "TMM",
                group_column = "group", design_terms = "group", ...) {

    ## Get the model to use
    option <- match.arg(option)

    ## Extract count matrix, sample information, node information
    count <- SummarizedExperiment::assay(SE, assay)
    sp_info <- SummarizedExperiment::colData(SE)

    ## Smallest number of samples in a group
    sp_min_0 <- min(table(sp_info[[group_column]]))

    ## Remove samples with zero count
    sp_keep <- colSums(count, na.rm = TRUE) > 0
    if (sum(!sp_keep)) {
        count <- count[, sp_keep]
        sp_info <- sp_info[sp_keep, , drop = FALSE]
    }
    sp_min_1 <- min(table(sp_info[[group_column]]))

    ## Design matrix
    ## -------------------------------------------------------------------------
    if (is.null(design)) {
        formula <- stats::as.formula(paste("~", paste(design_terms,
                                                      collapse = "+")))
        design <- stats::model.matrix(formula, data = data.frame(sp_info))
    }

    ## Filter lowly abundant features
    ## -------------------------------------------------------------------------
    ## Samples with zero library size are removed manually above, so we need to
    ## recalculate the min.prop here
    filter_min_prop <- filter_min_prop * sp_min_0 / sp_min_1
    keep <- edgeR::filterByExpr(count, design = design,
                                min.count = filter_min_count,
                                min.total.count = filter_min_total_count,
                                large.n = filter_large_n,
                                min.prop = filter_min_prop)

    count_keep <- count[keep, , drop = FALSE]
    isLow <- !keep
    feature_drop <- rownames(count)[isLow]

    ## Run edgeR
    ## -------------------------------------------------------------------------
    lrt <- edgerWrp(count = count_keep, lib_size = NULL,
                    option = option,
                    design = design, contrast = contrast,
                    normalize = normalize,
                    normalize_method = normalize_method, ...)

    ## Return result
    ## -------------------------------------------------------------------------
    return(lrt)
}


