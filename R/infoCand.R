#' Get information of candidates
#'
#' Extract information about candidates.
#'
#' @author Ruizhu Huang
#' @export
#'
#' @param object An output object from \link{evalCand}.
#'
#' @returns A \code{data.frame} with information about candidates.
#'
#' @importFrom rlang .data
#'
#' @examples
#' suppressPackageStartupMessages({
#'     library(TreeSummarizedExperiment)
#'     library(ggtree)
#' })
#'
#' ## Simulate some data
#' data(tinyTree)
#' ggtree(tinyTree, branch.length = "none") +
#'    geom_text2(aes(label = node)) +
#'    geom_hilight(node = 13, fill = "blue", alpha = 0.3) +
#'    geom_hilight(node = 18, fill = "orange", alpha = 0.3)
#' set.seed(1)
#' pv <- runif(19, 0, 1)
#' pv[c(seq_len(5), 13, 14, 18)] <- runif(8, 0, 0.001)
#'
#' fc <- sample(c(-1, 1), 19, replace = TRUE)
#' fc[c(seq_len(3), 13, 14)] <- 1
#' fc[c(4, 5, 18)] <- -1
#' df <- data.frame(node = seq_len(19),
#'                  pvalue = pv,
#'                  logFoldChange = fc)
#'
#' ## Get candidates
#' ll <- getCand(tree = tinyTree, score_data = df,
#'                node_column = "node",
#'                p_column = "pvalue",
#'                sign_column = "logFoldChange")
#'
#' ## Evaluate candidates
#' cc <- evalCand(tree = tinyTree, levels = ll$candidate_list,
#'                score_data = df, node_column = "node",
#'                p_column = "pvalue", sign_column = "logFoldChange",
#'                limit_rej = 0.05)
#'
#' ## Get summary info about candidates
#' out <- infoCand(object = cc)
#' out
#'
infoCand <- function(object) {

    ## Check arguments
    ## -------------------------------------------------------------------------
    .assertVector(x = object, type = "list")

    ## Extract level_info from list
    ## -------------------------------------------------------------------------
    info <- object$level_info
    if (is.null(info)) {
        stop("Object needs to have a 'level_info' slot - make sure that you ",
             "are using output from 'evalCand()'")
    }

    ## Remove rej_pseudo_leaf and rej_pseudo_node columns if all values are NA
    ## -------------------------------------------------------------------------
    allNA <- all(is.na(info$rej_pseudo_leaf))
    if (allNA) {
        info <- info |>
            dplyr::select(-dplyr::all_of(c("rej_pseudo_leaf",
                                           "rej_pseudo_node")))
    }

    return(info)
}
