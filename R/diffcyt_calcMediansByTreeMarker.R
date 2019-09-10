#' Calculate medians (by nodes (clusters) of a tree and marker)
#' 
#' \code{calcMediansByTreeMarker} calculates medians for each cluster-marker
#' combination. A cluster is represented by a node on the tree. This is a tree
#' version of \code{\link[diffcyt]{calcMediansByClusterMarker}}. The clusters
#' used in \code{\link[diffcyt]{calcMediansByClusterMarker}} are clusters on the
#' leaf level of the tree here. More details about the data could be found in
#' \code{\link[diffcyt]{calcMediansByClusterMarker}}.
#' 
#' 
#' @param d_se Data object from previous steps, in
#'   \code{\link[TreeSummarizedExperiment]{TreeSummarizedExperiment}} format,
#'   containing cluster labels as a column in the row meta-data (from
#'   \code{\link[diffcyt]{generateClusters}}). Column meta-data is assumed to contain a
#'   factor \code{marker_class}.
#' @param tree A phylo object from \code{\link{buildTree}}
#' 
#' 
#' @return a \code{\link[TreeSummarizedExperiment]{TreeSummarizedExperiment}}
#'   object. Clusters (or nodes of a tree) are in rows and markers in columns.
#'   The marker expression values are in \code{assay} and the \code{metadata}
#'   slot contains variables \code{id_type_markers}
#' 
#' 
#' @import TreeSummarizedExperiment
#' @importFrom methods is
#' @importFrom SummarizedExperiment rowData
#' @export
#' 
#' @author Ruizhu Huang
#' 
#' @examples
#' # For a complete workflow example demonstrating each step in the 'diffcyt' pipeline, 
#' # see the package vignette.
#' library(diffcyt)
#' 
#' # Function to create random data (one sample)
#' d_random <- function(n = 20000, mean = 0, sd = 1, ncol = 20, cofactor = 5) {
#'   d <- sinh(matrix(rnorm(n, mean, sd), ncol = ncol)) * cofactor
#'   colnames(d) <- paste0("marker", sprintf("%02d", 1:ncol))
#'   d
#' }
#' 
#' # Create random data (without differential signal)
#' set.seed(123)
#' d_input <- list(
#'   sample1 = d_random(), 
#'   sample2 = d_random(), 
#'   sample3 = d_random(), 
#'   sample4 = d_random()
#' )
#' 
#' experiment_info <- data.frame(
#'   sample_id = factor(paste0("sample", 1:4)), 
#'   group_id = factor(c("group1", "group1", "group2", "group2")), 
#'   stringsAsFactors = FALSE
#' )
#' 
#' marker_info <- data.frame(
#'   channel_name = paste0("channel", sprintf("%03d", 1:20)), 
#'   marker_name = paste0("marker", sprintf("%02d", 1:20)), 
#'   marker_class = factor(c(rep("type", 10), rep("state", 10)), 
#'                         levels = c("type", "state", "none")), 
#'   stringsAsFactors = FALSE
#' )
#' 
#' # Prepare data
#' d_se <- prepareData(d_input, experiment_info, marker_info)
#' 
#' # Transform data
#' d_se <- transformData(d_se)
#' 
#' # Generate clusters
#' d_se <- generateClusters(d_se)
#' 
#' # build a tree
#' tr <- buildTree(d_se)
#' 
#' # Calculate medians (by cluster and marker)
#' d_medians_by_cluster_marker <- calcMediansByTreeMarker(d_se, tr)
#' 

calcMediansByTreeMarker <- function(d_se, tree) {
    
    if (!("cluster_id" %in% (colnames(rowData(d_se))))) {
        stop("Data object does not contain cluster labels. 
             Run 'diffcyt::generateClusters' to generate cluster labels.")
    }
    
    if (!is(tree, "phylo")) {
        stop("tree is not a phylo object.
             Run 'buildTree(d_se)' to generate the tree")
    }
    
    ## create the tse object
    rlab <- as.character(rowData(d_se)$cluster_id)
    d_lse <- TreeSummarizedExperiment(assays = assays(d_se),
                                      rowData = rowData(d_se),
                                      rowTree = tree, 
                                      rowNodeLab = rlab,
                                      colData = colData(d_se),
                                      metadata = metadata(d_se))
    d_lse <- d_lse[, colData(d_lse)$marker_class %in% c("type", "state")]
    
    # calculate median value at each node
    nodes <- showNode(tree = tree, only.leaf = FALSE, use.alias = FALSE)
    d_tse <- aggValue(x = d_lse, rowLevel = nodes, 
                      FUN = function(x){
                          median(x, na.rm = TRUE)
                      })
    rownames(d_tse) <- rowLinks(d_tse)$nodeLab
    return(d_tse)
}