#' calculates number of cells per cluster-sample combination.
#' 
#' \code{calcTreeCounts} calculates number of cells per cluster-sample
#' combination (referred to as cluster cell 'counts', 'abundances', or
#' 'frequencies'). This is a tree version of \code{\link[diffcyt]{calcCounts}}.
#' The clusters used in \code{\link[diffcyt]{calcCounts}} are clusters on the
#' leaf level of the tree here. More details about the data could be found in
#' \code{\link[diffcyt]{calcCounts}}.
#' 
#' 
#' @param d_se Data object from previous steps, in
#'   \code{\link[TreeSummarizedExperiment]{TreeSummarizedExperiment}} format,
#'   containing cluster labels as a column in the row meta-data (from
#'   \code{\link[diffcyt]{generateClusters}}). Column meta-data is assumed to
#'   contain a factor marker_class.
#' @param tree A phylo object from \code{\link{buildTree}}
#' 
#' @importFrom dplyr mutate "%>%" 
#' @import SummarizedExperiment
#' @import TreeSummarizedExperiment
#' @importFrom methods is 
#' @export
#' @author Ruizhu Huang
#' @return A TreeSummarizedExperiment object, where clusters in rows, samples in
#'   columns and counts in \code{assays} 
#' @examples 
#' # For a complete workflow example demonstrating each step, please see the
#' # vignette of 'diffcyt'
#' \dontrun{
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
#' # calculate medians on nodes of a tree
#' d_medians_tree <- calcTreeCounts(d_se = d_se, tree = tr)
#' }

calcTreeCounts <- function(d_se, tree) {
    
    if (!("cluster_id" %in% (colnames(rowData(d_se))))) {
        stop("Data object does not contain cluster labels. 
             Run 'diffcyt::generateClusters' to generate cluster labels.")
    }
    
    if (!is(tree, "phylo")) {
        stop("tree is not a phylo object.
             Run 'buildTree(d_se)' to generate the tree")
    }
    
    # counts on the leaf level
    d_counts <- .calcCounts(d_se)
    
    # build a TreeSummarizedExperiment object
    rlab <- as.character(rowData(d_counts)$cluster_id)
    counts_leaf <- TreeSummarizedExperiment(assays = assays(d_counts),
                                              rowData = rowData(d_counts),
                                              rowTree = tree, 
                                              rowNodeLab = rlab,
                                              colData = colData(d_counts),
                                              metadata = metadata(d_counts))
    
    # counts on all nodes
    nodes <- showNode(tree = tree, only.leaf = FALSE)
    counts_all <- aggTSE(x = counts_leaf, rowLevel = nodes, 
                         rowFun = function(x) {sum(x, na.rm = TRUE)})
    
    lab <- convertNode(tree = tree, 
                     node = rowLinks(counts_all)$nodeNum, use.alias = TRUE)
    rowData(counts_all) <- rowData(counts_all) %>%
        data.frame() %>%
        mutate(cluster_id = factor(lab, levels = lab)) %>%
        mutate(n_cells = apply(assays(counts_all)[[1]], 1, sum))
    rownames(counts_all) <- lab

    return(counts_all)
    }

#' @importFrom reshape2 acast
#' @importFrom dplyr "%>%" group_by 
#' @importFrom tidyr complete
#' @import SummarizedExperiment

.calcCounts <- function(d_se) {
    if (!is(d_se, "SummarizedExperiment")) {
        stop("Data object must be a 'SummarizedExperiment'")
    }
    if (!("cluster_id" %in% (colnames(rowData(d_se))))) {
        stop("Data object does not contain cluster labels. Run 'generateClusters' to generate cluster labels.")
    }
    rowdata_df <- as.data.frame(rowData(d_se))
    n_cells <- rowdata_df %>% group_by(cluster_id, sample_id, 
                                       .drop = FALSE) %>% tally %>% complete(sample_id)
    n_cells <- acast(n_cells, cluster_id ~ sample_id, value.var = "n", 
                     fill = 0)
    if (nrow(n_cells) < nlevels(rowData(d_se)$cluster_id)) {
        ix_missing <- which(!(levels(rowData(d_se)$cluster_id) %in% 
                                  rownames(n_cells)))
        n_cells_tmp <- matrix(0, nrow = length(ix_missing), ncol = ncol(n_cells))
        rownames(n_cells_tmp) <- ix_missing
        n_cells <- rbind(n_cells, n_cells_tmp)
        n_cells <- n_cells[order(as.numeric(rownames(n_cells))), 
                           , drop = FALSE]
    }
    n_cells_total <- rowSums(n_cells)
    row_data <- data.frame(cluster_id = factor(rownames(n_cells), 
                                               levels = levels(rowData(d_se)$cluster_id)), n_cells = n_cells_total, 
                           stringsAsFactors = FALSE)
    col_data <- metadata(d_se)$experiment_info
    n_cells <- n_cells[, match(col_data$sample_id, colnames(n_cells)), 
                       drop = FALSE]
    stopifnot(all(col_data$sample_id == colnames(n_cells)))
    d_counts <- SummarizedExperiment(assays = list(counts = n_cells), 
                                     rowData = row_data, colData = col_data)
    d_counts
}