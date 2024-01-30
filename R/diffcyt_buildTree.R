#' Build a tree from diffcyt high-resolution cell clusters
#'
#' Apply hierarchical clustering to build a tree starting from a high-resolution
#' clustering created by the \code{\link[diffcyt]{generateClusters}} function
#' from the \code{diffcyt} package. The function calculates the median
#' abundance for each (ID type) marker and cluster, and uses this data to
#' further aggregate the initial clusters using hierarchical clustering.
#'
#' @author Ruizhu Huang
#' @export
#'
#' @param d_se A \code{SummarizedExperiment} object, with cells as rows and
#'     features as columns. This should be the output
#'     from \code{\link[diffcyt]{generateClusters}}. The \code{colData} is
#'     assumed to contain a factor named \code{marker_class}.
#' @param dist_method The distance measure to be used. This must be one of
#'     "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski".
#'     Any unambiguous substring can be given. Please refer to \code{method} in
#'     \code{\link[stats]{dist}} for more information.
#' @param hclust_method The agglomeration method to be used. This should be (an
#'     unambiguous abbreviation of) one of "ward.D", "ward.D2", "single",
#'     "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC)
#'     or "centroid" (= UPGMC). Please refer to \code{method} in
#'     \code{\link[stats]{hclust}} for more information.
#'
#' @importFrom diffcyt calcMediansByClusterMarker
#' @importFrom stats hclust dist
#' @importFrom ape as.phylo
#' @importFrom SummarizedExperiment assay
#' @importFrom TreeSummarizedExperiment addLabel
#' @importFrom S4Vectors metadata
#'
#' @return A \code{phylo} object representing the hierarchical clustering of
#' the initial high-resolution clusters.
#'
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
#'     sample1 = d_random(), sample2 = d_random(),
#'     sample3 = d_random(), sample4 = d_random()
#' )
#'
#' experiment_info <- data.frame(
#'     sample_id = factor(paste0("sample", seq_len(4))),
#'     group_id = factor(c("group1", "group1", "group2", "group2"))
#' )
#'
#' marker_info <- data.frame(
#'     channel_name = paste0("channel", sprintf("%03d", seq_len(20))),
#'     marker_name = paste0("marker", sprintf("%02d", seq_len(20))),
#'     marker_class = factor(c(rep("type", 10), rep("state", 10)),
#'                           levels = c("type", "state", "none"))
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
#' # Build a tree
#' tr <- buildTree(d_se)
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
    tree_h <- stats::hclust(dist(md, method = dist_method),
                            method = hclust_method)

    tree_p <- ape::as.phylo(tree_h)
    tree_p <- TreeSummarizedExperiment::addLabel(tree_p, on = "internal")

    return(tree_p)
}
