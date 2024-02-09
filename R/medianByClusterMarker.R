#' Calculate median values of markers for each cluster
#'
#' Calculate median value of each marker in each cluster.
#'
#' @author Ruizhu Huang, Charlotte Soneson
#' @export
#'
#' @param SE A \code{\link[SummarizedExperiment]{SummarizedExperiment}} object.
#' @param assay A numeric index or assay name indicating with assay of \code{SE}
#'     to use to calculate medians.
#' @param marker_in_column A logical scalar, indicating whether markers (genes,
#'     features) are in the columns of \code{SE} or not. The default is
#'     \code{TRUE} (markers are in columns).
#' @param column_cluster The name of the column of \code{colData(SE)} that
#'     contains the cluster assignment of each smaple.
#' @param use_marker A logical or numeric vector such that
#'     \code{SE[use_marker, ]} (if \code{marker_in_column = FALSE}) or
#'     \code{SE[, use_marker]} (if \code{marker_in_column = TRUE}) subsets
#'     \code{SE} to the markers that should be retained. If \code{NULL}
#'     (default), all markers are used.
#'
#' @returns A \code{\link[SummarizedExperiment]{SummarizedExperiment}} object
#'     containing the median value of each marker in each cluster.
#'
#' @importFrom dplyr mutate summarize_all group_by
#' @importFrom tibble column_to_rownames
#' @importFrom SummarizedExperiment assays rowData colData SummarizedExperiment
#'     assayNames
#'
#' @examples
#' library(SummarizedExperiment)
#'
#' ## Simulate data with 100 cells and 10 markers (5 type, 5 state markers)
#' set.seed(1)
#' count <- matrix(rpois(n = 1000, lambda = 10), nrow = 100)
#' colnames(count) <- paste0("mk", 1:10)
#' rowD <- data.frame("cluster" = sample(seq_len(6), 100, replace = TRUE))
#' colD <- data.frame(type_marker = rep(c(FALSE, TRUE), each = 5))
#'
#' ## SE with markers in columns
#' d_se <- SummarizedExperiment(assays = list(counts = count),
#'                              rowData = rowD,
#'                              colData = colD)
#' medianByClusterMarker(SE = d_se, marker_in_column = TRUE,
#'                       column_cluster = "cluster",
#'                       use_marker = colData(d_se)$type_marker)
#'
#' ## SE with markers in rows
#' d_se <- SummarizedExperiment(assays = list(counts = t(count)),
#'                              rowData = colD,
#'                              colData = rowD)
#' medianByClusterMarker(SE = d_se, marker_in_column = FALSE,
#'                       column_cluster = "cluster",
#'                       use_marker = rowData(d_se)$type_marker)
#'
medianByClusterMarker <- function(SE, assay = 1, marker_in_column = TRUE,
                                  column_cluster = "cluster_id",
                                  use_marker = NULL) {

    ## Check arguments
    ## -------------------------------------------------------------------------
    .assertVector(x = SE, type = "SummarizedExperiment")
    stopifnot(length(assay) == 1 && (is.numeric(assay) || is.character(assay)))
    if (is.character(assay)) {
        .assertScalar(x = assay, type = "character",
                      validValues = SummarizedExperiment::assayNames(SE))
    } else if (is.numeric(assay)) {
        .assertScalar(
            x = assay, type = "numeric",
            rngIncl = c(1, length(SummarizedExperiment::assays(SE))))
    }
    .assertScalar(x = marker_in_column, type = "logical")
    .assertScalar(x = column_cluster, type = "character")
    if (!is.null(use_marker)) {
        stopifnot(is.logical(use_marker) || is.numeric(use_marker) ||
                      is.character(use_marker))
    }

    ## Subset data to get data for markers
    ## -------------------------------------------------------------------------
    if (marker_in_column) {
        if (is.null(use_marker)) {
            use_marker <- rep(TRUE, ncol(SE))
        }
        SE <- SE[, use_marker]

        ## Get data (markers in columns)
        data_mk <- SummarizedExperiment::assays(SE)[[assay]]
        ## Get cluster id
        cluster_id <- as.character(
            SummarizedExperiment::rowData(SE)[[column_cluster]])
    } else {
        if (is.null(use_marker)) {
            use_marker <- rep(TRUE, nrow(SE))
        }
        SE <- SE[use_marker, ]

        ## Get data (markers in columns)
        data_mk <- t(SummarizedExperiment::assays(SE)[[assay]])
        ## Get cluster id
        cluster_id <- as.character(
            SummarizedExperiment::colData(SE)[[column_cluster]])
    }

    ## Calculate the median
    ## -------------------------------------------------------------------------
    med <- data_mk |>
        data.frame() |>
        mutate(cluster_id = cluster_id) |>
        group_by(cluster_id) |>
        summarize_all(median, na.rm = TRUE) |>
        column_to_rownames("cluster_id") |>
        as.matrix()

    rowD <- data.frame(column_cluster = rownames(med))
    colnames(rowD) <- column_cluster

    ## Assemble output object
    ## -------------------------------------------------------------------------
    if (marker_in_column) {
        out <- SummarizedExperiment::SummarizedExperiment(
            assays = list(med),
            rowData = rowD)
    } else {
        out <- SummarizedExperiment::SummarizedExperiment(
            assays = list(t(med)),
            colData = rowD)
    }
    if (is.character(assay)) {
        SummarizedExperiment::assayNames(out) <- assay
    } else if (is.numeric(assay)) {
        ## Check if there is an assay name (even if the index was specified)
        anm <- SummarizedExperiment::assayNames(SE)[assay]
        if (!is.null(anm) && anm != "") {
            SummarizedExperiment::assayNames(out) <- anm
        }
    }

    return(out)
}
