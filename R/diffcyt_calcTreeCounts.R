#' @export
#' @rdname diffcyt_workflow
#'
#' @importFrom dplyr mutate
#' @importFrom TreeSummarizedExperiment TreeSummarizedExperiment showNode
#'     aggTSE convertNode rowLinks
#' @importFrom methods is
#' @importFrom SummarizedExperiment rowData assays colData
#' @importFrom S4Vectors metadata
#' @importFrom diffcyt calcCounts
#'
calcTreeCounts <- function(d_se, tree) {

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

    ## Calculate counts on the leaf level
    ## -------------------------------------------------------------------------
    d_counts <- diffcyt::calcCounts(d_se)

    ## Build a TreeSummarizedExperiment object
    ## -------------------------------------------------------------------------
    rlab <- as.character(SummarizedExperiment::rowData(d_counts)$cluster_id)
    counts_leaf <- TreeSummarizedExperiment::TreeSummarizedExperiment(
        assays = SummarizedExperiment::assays(d_counts),
        rowData = SummarizedExperiment::rowData(d_counts),
        rowTree = tree,
        rowNodeLab = rlab,
        colData = SummarizedExperiment::colData(d_counts),
        metadata = S4Vectors::metadata(d_counts))

    ## Get counts on all nodes
    ## -------------------------------------------------------------------------
    nodes <- TreeSummarizedExperiment::showNode(tree = tree, only.leaf = FALSE)
    counts_all <- TreeSummarizedExperiment::aggTSE(
        x = counts_leaf, rowLevel = nodes,
        rowFun = function(x) sum(x, na.rm = TRUE))

    ## Construct final TSE object
    ## -------------------------------------------------------------------------
    lab <- TreeSummarizedExperiment::convertNode(
        tree = tree,
        node = TreeSummarizedExperiment::rowLinks(counts_all)$nodeNum,
        use.alias = TRUE)

    SummarizedExperiment::rowData(counts_all) <-
        SummarizedExperiment::rowData(counts_all) |>
        data.frame() |>
        mutate(cluster_id = factor(lab, levels = lab)) |>
        mutate(n_cells = apply(SummarizedExperiment::assays(counts_all)[[1]],
                               1, sum))
    rownames(counts_all) <- lab

    return(counts_all)
}
