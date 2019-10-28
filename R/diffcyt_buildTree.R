#' build a tree from clusters 
#' \code{treeBuild} creates a tree (\code{phylo} class) 
#' 
#' @param d_se the data object output from
#'   \code{\link[diffcyt]{generateClusters}}
#' @param dist_method the distance measure to be used. This must be one of
#'   "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski".
#'   Any unambiguous substring can be given. Please refer to \code{method} in
#'   \code{\link[stats]{dist}}
#' @param hclust_method the agglomeration method to be used. This should be (an
#'   unambiguous abbreviation of) one of "ward.D", "ward.D2", "single",
#'   "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC)
#'   or "centroid" (= UPGMC). Please refer to \code{method} in
#'   \code{\link[stats]{hclust}}
#' 
#' @importFrom diffcyt calcMediansByClusterMarker
#' @importFrom stats hclust dist
#' @importFrom ape as.phylo
#' @export
#' @return A phylo object
#' @examples 
#' # For a complete workflow example demonstrating each step in the 'diffcyt'
#' # pipeline, please see the diffcyt vignette.
#' library(diffcyt)
#' 
#' # Function to create random data (one sample)
#' d_random <- function(n = 20000, mean = 0, sd = 1, ncol = 20, cofactor = 5) {
#'     d <- sinh(matrix(rnorm(n, mean, sd), ncol = ncol)) * cofactor
#'     colnames(d) <- paste0("marker", sprintf("%02d", 1:ncol))
#'     d
#' }
#' 
#' # Create random data (without differential signal)
#' set.seed(123)
#' d_input <- list(
#'     sample1 = d_random(), 
#'     sample2 = d_random(), 
#'     sample3 = d_random(), 
#'     sample4 = d_random()
#' )
#' 
#' experiment_info <- data.frame(
#'     sample_id = factor(paste0("sample", 1:4)), 
#'     group_id = factor(c("group1", "group1", "group2", "group2")), 
#'     stringsAsFactors = FALSE
#' )
#' 
#' marker_info <- data.frame(
#'     channel_name = paste0("channel", sprintf("%03d", 1:20)), 
#'     marker_name = paste0("marker", sprintf("%02d", 1:20)), 
#'     marker_class = factor(c(rep("type", 10), rep("state", 10)), 
#'                           levels = c("type", "state", "none")), 
#'     stringsAsFactors = FALSE
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
buildTree <- function(d_se, dist_method = "euclidean",
                      hclust_method = "average") {
    d_medians <- calcMediansByClusterMarker(d_se)
    md <- assay(d_medians)[, metadata(d_medians)$id_type_markers]
    tree_h <- hclust(dist(md, method = dist_method), method = hclust_method)
    
    tree_p <- as.phylo(tree_h)
    tree_p <- addLabel(tree_p)
    return(tree_p)
}
