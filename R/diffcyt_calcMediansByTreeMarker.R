#' @export
#' @rdname diffcyt_workflow
#'
#' @importFrom TreeSummarizedExperiment aggTSE showNode
#'     TreeSummarizedExperiment rowLinks
#' @importFrom methods is
#' @importFrom SummarizedExperiment rowData colData assays
#' @importFrom S4Vectors metadata
#' @importFrom stats median
#'
calcMediansByTreeMarker <- function(d_se, tree) {

    ## Check arguments
    ## -------------------------------------------------------------------------
    .assertVector(x = d_se, type = "SummarizedExperiment")
    if (!("cluster_id" %in% (colnames(SummarizedExperiment::rowData(d_se))))) {
        stop("Data object does not contain cluster labels.
             Run 'diffcyt::generateClusters' to generate cluster labels.")
    }
    if (!methods::is(tree, "phylo")) {
        stop("tree is not a phylo object.
             Run 'buildTree(d_se)' to generate the tree.")
    }

    ## Create the TSE object
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

    ## Calculate median values at each node
    ## -------------------------------------------------------------------------
    nodes <- TreeSummarizedExperiment::showNode(
        tree = tree, only.leaf = FALSE, use.alias = FALSE)
    d_tse <- TreeSummarizedExperiment::aggTSE(x = d_lse, rowLevel = nodes,
                                              rowFun = function(x) {
                                                  stats::median(x, na.rm = TRUE)
                                              })
    rownames(d_tse) <- TreeSummarizedExperiment::rowLinks(d_tse)$nodeLab

    return(d_tse)
}
