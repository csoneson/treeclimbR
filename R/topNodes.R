#' Generate a table of top-ranked entities (nodes)
#'
#' Generate a table of top-ranked nodes from the optimal resolution candidate
#' of entities on a tree.
#'
#' @author Ruizhu Huang, Charlotte Soneson
#' @export
#'
#' @param object An output object from \link{evalCand}.
#' @param n An integer, the maximum number of entities to return.
#' @param sort_by A character string specifying the column of
#'     \code{object$output} to sort by. Set to \code{NULL} to return without
#'     sorting.
#' @param sort_decreasing A logical value indicating whether to sort by
#'     decreasing value of the \code{sort_by} column.
#' @param sort_by_absolute A logical value indicating whether to take the
#'     absolute value of the \code{sort_by} column before sorting.
#' @param p_value A numeric cutoff value for adjusted p-values. Only entities
#'     with adjusted p-values equal or lower than specified are returned.
#'
#' @returns A \code{data.table} with test results. The \strong{node}
#'     column stores the node number for each entity.
#'
#' @importFrom dplyr arrange slice filter desc
#' @importFrom rlang .data
#'
#' @examples
#' library(TreeSummarizedExperiment)
#' library(ggtree)
#'
#' data(tinyTree)
#' ggtree(tinyTree, branch.length = "none") +
#'    geom_text2(aes(label = node)) +
#'    geom_hilight(node = 13, fill = "blue", alpha = 0.5) +
#'    geom_hilight(node = 18, fill = "orange", alpha = 0.5)
#' set.seed(2)
#' pv <- runif(19, 0, 1)
#' pv[c(seq_len(5), 13, 14, 18)] <- runif(8, 0, 0.001)
#'
#' fc <- sample(c(-1, 1), 19, replace = TRUE)
#' fc[c(seq_len(3), 13, 14)] <- 1
#' fc[c(4, 5, 18)] <- -1
#' df <- data.frame(node = seq_len(19),
#'                  pvalue = pv,
#'                  logFoldChange = fc)
#' ll <- getCand(tree = tinyTree, score_data = df,
#'                node_column = "node",
#'                p_column = "pvalue",
#'                sign_column = "logFoldChange")
#' cc <- evalCand(tree = tinyTree, levels = ll$candidate_list,
#'                score_data = df, node_column = "node",
#'                p_column = "pvalue", sign_column = "logFoldChange",
#'                limit_rej = 0.05)
#'
#' ## Unsorted result table
#' topNodes(cc)
#'
#' ## Sort by p-value in increasing order
#' topNodes(cc, sort_by = "pvalue")
#'
topNodes <- function(object, n = 10, sort_by = NULL,
                     sort_decreasing = FALSE,
                     sort_by_absolute = FALSE, p_value = 1) {

    ## Check arguments
    ## -------------------------------------------------------------------------
    .assertVector(x = object, type = "list")
    stopifnot("output" %in% names(object))
    .assertVector(x = object$output, type = "data.frame")
    stopifnot("adj.p" %in% colnames(object$output))
    .assertScalar(x = sort_by, type = "character", allowNULL = TRUE)
    if (!is.null(sort_by)) {
        stopifnot(sort_by %in% colnames(object$output))
    }
    .assertScalar(x = n, type = "numeric")
    .assertScalar(x = sort_decreasing, type = "logical")
    .assertScalar(x = sort_by_absolute, type = "logical")
    .assertScalar(x = p_value, type = "numeric", rngIncl = c(0, 1))

    ## Get the result table from the object, sort and subset
    ## -------------------------------------------------------------------------
    res <- object$output
    n <- min(nrow(res), n)

    if (is.null(sort_by)) {
        res <- res |>
            dplyr::filter(.data$adj.p <= p_value) |>
            dplyr::slice(seq_len(n))
    } else if (!sort_decreasing) {
        if (sort_by_absolute) {
            res <- res |>
                dplyr::arrange(abs(.data[[sort_by]])) |>
                dplyr::filter(.data$adj.p <= p_value) |>
                dplyr::slice(seq_len(n))
        } else {
            res <- res |>
                dplyr::arrange(.data[[sort_by]]) |>
                dplyr::filter(.data$adj.p <= p_value) |>
                dplyr::slice(seq_len(n))
        }
    } else if (sort_decreasing) {
        if (sort_by_absolute) {
            res <- res |>
                dplyr::arrange(dplyr::desc(abs(.data[[sort_by]]))) |>
                dplyr::filter(.data$adj.p <= p_value) |>
                dplyr::slice(seq_len(n))
        } else {
            res <- res |>
                dplyr::arrange(dplyr::desc(.data[[sort_by]])) |>
                dplyr::filter(.data$adj.p <= p_value) |>
                dplyr::slice(seq_len(n))
        }
    }

    return(res)
}
