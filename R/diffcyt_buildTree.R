#' @export
#' @rdname diffcyt_workflow
#'
#' @importFrom diffcyt calcMediansByClusterMarker
#' @importFrom stats hclust dist
#' @importFrom ape as.phylo
#' @importFrom SummarizedExperiment assay
#' @importFrom TreeSummarizedExperiment addLabel
#' @importFrom S4Vectors metadata
#'
buildTree <- function(d_se, dist_method = "euclidean",
                      hclust_method = "average") {

    ## Check arguments
    ## -------------------------------------------------------------------------
    .assertVector(x = d_se, type = "SummarizedExperiment")
    .assertScalar(x = dist_method, type = "character")
    .assertScalar(x = hclust_method, type = "character")

    ## Build tree
    ## -------------------------------------------------------------------------
    d_medians <- diffcyt::calcMediansByClusterMarker(d_se)
    md <- SummarizedExperiment::assay(d_medians)[, S4Vectors::metadata(
        d_medians)$id_type_markers]
    tree_h <- stats::hclust(d = stats::dist(md, method = dist_method),
                            method = hclust_method)

    tree_p <- ape::as.phylo(tree_h)
    tree_p <- TreeSummarizedExperiment::addLabel(tree_p, on = "internal")

    return(tree_p)
}
