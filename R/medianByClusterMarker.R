#' calculate median values of markers for each cluster
#' 
#' \code{medianByClusterMarker} calculates median value of each marker in each
#' cluster. The input data should be a
#' \code{\link[SummarizedExperiment]{SummarizedExperiment}} object. The markers
#' could either be in columns or in rows of the \code{assays} slot.
#' 
#' @param SE A \code{\link[SummarizedExperiment]{SummarizedExperiment}} object.
#' @param marker_in_column A logical value, either TRUE or FALSE. The default is
#'   TRUE (markers are in columns)
#' @param column_cluster The name of the column that has information about
#'   cluster id.
#' @param use_marker A logical vector including FALSE and TRUE to specify
#'   markers to be used. The default is NULL and all markers are used.
#' @importFrom dplyr '%>%' mutate summarize_all group_by
#' @importFrom tibble column_to_rownames
#' @importFrom SummarizedExperiment assays rowData colData SummarizedExperiment
#' 
#' @export
#' @return A \code{\link[SummarizedExperiment]{SummarizedExperiment}} object. It
#'   includes clusters in rows and the median value of each marker in each
#'   cluster in columns
#' @author Ruizhu Huang
#' @examples 
#' library(TreeSummarizedExperiment)
#' set.seed(1)
#' count <- matrix(rpois(n = 1000, lambda = 10),
#'                 nrow = 100)
#' colnames(count) <- paste0("mk", 1:10)
#' rowD <- data.frame("cluster" = sample(1:6, 100, replace = TRUE))
#' colD <- data.frame(type_marker = rep(c(FALSE, TRUE), each = 5))
#' d_se <- TreeSummarizedExperiment(assays = list(count),
#'                                  rowData = rowD,
#'                                  colData = colD)
#' d_med <- medianByClusterMarker(SE = d_se, marker_in_column = TRUE,
#'                                column_cluster = "cluster",
#'                                use_marker = colData(d_se)$type_marker)
#' d_med
medianByClusterMarker <- function(SE, marker_in_column = TRUE, 
                                  column_cluster = "cluster_id",
                                  use_marker = NULL) {
    # subset data to get data for markers
    if (marker_in_column) {
        if (is.null(use_marker)) { use_marker <- rep(TRUE, ncol(SE)) }
        SE <- SE[, use_marker]
        # data of type markers
        data_mk <- assays(SE)[[1]]
        # cluster_id
        cluster_id <- rowData(SE)[[column_cluster]]
    } else {
        if (is.null(use_marker)) { use_marker <- rep(TRUE, nrow(SE)) }
        SE <- SE[use_marker, ]
        # data of type marker features
        data_mk <- t(assays(SE)[[1]])
        # cluster_id
        cluster_id <- colData(SE)[[column_cluster]]
    }
   # calculate the median
    med <- data_mk %>%
        data.frame() %>%
        mutate(cluster_id = cluster_id) %>%
        group_by(cluster_id) %>%
        summarize_all(median, na.rm = TRUE) %>%
        column_to_rownames("cluster_id") %>%
        as.matrix()
    
    rowD <- data.frame(column_cluster = rownames(med),
                       stringsAsFactors = FALSE)
    colnames(rowD) <- column_cluster
    # output 
    out <- SummarizedExperiment(
        assays = list(med), 
        rowData = rowD)
    return(out)
    
}
