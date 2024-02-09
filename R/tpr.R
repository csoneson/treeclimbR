#' Calculate true positive rate (TPR) on a tree structure
#'
#' Calculate the true positive rate on a tree structure, at either leaf or
#' node level.
#'
#' @author Ruizhu Huang
#' @export
#'
#' @param tree A \code{phylo} object.
#' @param truth True signal nodes (e.g., nodes that are truly differentially
#'     abundant between experimental conditions). \strong{Note:} when the
#'     TPR is requested at the leaf level (\code{only.leaf = TRUE}), the
#'     descendant leaves of the given nodes will be found and the TPR will be
#'     estimated on the leaf level.
#' @param found Detected signal nodes (e.g., nodes that have been found to be
#'     differentially abundant via a statistical testing procedure).
#'     \strong{Note:} when the TPR is requested at the leaf level
#'     (\code{only.leaf = TRUE}), the descendant leaves of the given nodes will
#'     be found out and the TPR will be estimated on the leaf level.
#' @param only.leaf A logical scalar. If \code{TRUE}, the false discovery rate
#'     is calculated at the leaf (tip) level; otherwise it is calculated at
#'     the node level.
#'
#' @returns The estimated true positive rate.
#'
#' @importFrom TreeSummarizedExperiment convertNode
#'
#' @examples
#' suppressPackageStartupMessages({
#'     library(ggtree)
#'     library(TreeSummarizedExperiment)
#' })
#'
#' data("tinyTree")
#'
#' ## Two branches are truly differential
#' ggtree(tinyTree) +
#'    geom_text2(aes(label = node)) +
#'    geom_hilight(node = 16, fill = "orange", alpha = 0.3) +
#'    geom_hilight(node = 13, fill = "blue", alpha = 0.3)
#'
#' ## TPR at the leaf level (7/8)
#' tpr(tree = tinyTree, truth = c(16, 13),
#'     found = c(15, 14), only.leaf = TRUE)
#'
#' ## TPR at the node level (12/14)
#' tpr(tree = tinyTree, truth = c(16, 13),
#'     found = c(15, 14), only.leaf = FALSE)
#'
tpr <- function(tree, truth, found, only.leaf = TRUE) {
    .assertVector(x = tree, type = "phylo")
    .assertScalar(x = only.leaf, type = "logical")

    if (is.character(truth)) {
        truth <- TreeSummarizedExperiment::convertNode(
            tree = tree, node = truth, message = FALSE)
    }
    if (is.character(found)) {
        found <- TreeSummarizedExperiment::convertNode(
            tree = tree, node = found, message = FALSE)
    }

    tt <- .tpr0(tree = tree, truth = truth,
                found = found, only.leaf = only.leaf)
    tpr <- tt[1]/tt[2]
    names(tpr) <- "tpr"

    return(tpr)
}


#' Count true positives and all positives
#'
#' Calculate the number of true positives and the total
#' number of positives on a tree structure, at leaf or node level.
#'
#' @author Ruizhu Huang
#'
#' @param tree A \code{phylo} object.
#' @param truth True signal nodes (e.g., nodes that are truly differentially
#'     abundant between experimental conditions).
#' @param found Detected signal nodes (e.g., nodes that have been found to be
#'     differentially abundant via a statistical testing procedure).
#' @param only.leaf A logical scalar. If \code{TRUE}, the false discovery rate
#'     is calculated at the leaf (tip) level; otherwise it is calculated at
#'     the node level.
#'
#' @return A vector with two numbers, the number of true positives and the
#' total number of positives.
#'
#' @keywords internal
#' @noRd
#'
#' @importFrom TreeSummarizedExperiment findDescendant
#'
.tpr0 <- function(tree, truth = NULL, found = NULL, only.leaf = TRUE) {

    ## Check inputs
    ## -------------------------------------------------------------------------
    if (!is.null(truth) && !(is.character(truth) || is.numeric(truth))) {
        stop("'truth' should be either a character vector or a numeric ",
             "vector")
    }

    if (!is.null(found) && !(is.character(found) || is.numeric(found))) {
        stop("'found' should be either a character vector or a numeric ",
             "vector")
    }

    ## Count positives
    ## -------------------------------------------------------------------------
    if (is.null(truth) || length(truth) == 0) {
        ## No truly positive nodes
        c(tp = 1, pos = 1)
    } else {
        ## Some truly positive nodes
        nodeT <- TreeSummarizedExperiment::findDescendant(
            tree = tree, node = truth, only.leaf = only.leaf,
            self.include = TRUE)
        nodeT <- unique(unlist(nodeT))

        if (is.null(found)) {
            ## No significant nodes
            c(tp = 0, pos = length(nodeT))
        } else {
            ## Some significant nodes
            nodeF <- TreeSummarizedExperiment::findDescendant(
                tree = tree, node = found, only.leaf = only.leaf,
                self.include = TRUE)
            nodeF <- unique(unlist(nodeF))

            TP <- intersect(nodeT, nodeF)
            c(tp = length(TP), pos = length(nodeT))
        }
    }
}
