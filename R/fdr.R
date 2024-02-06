#' Calculate false discovery rate (FDR) on a tree structure
#'
#' Calculate the false discovery rate on a tree structure, at either leaf or
#' node level.
#'
#' @author Ruizhu Huang
#' @export
#'
#' @param tree A \code{phylo} object.
#' @param truth True signal nodes (e.g., nodes that are truly differentially
#'     abundant between experimental conditions). \strong{Note:} when the
#'     FDR is requested at the leaf level (\code{only.leaf = TRUE}), the
#'     descendant leaves of the given nodes will be found and the FDR will be
#'     estimated on the leaf level.
#' @param found Detected signal nodes (e.g., nodes that have been found to be
#'     differentially abundant via a statistical testing procedure).
#'     \strong{Note:} when the FDR is requested at the leaf level
#'     (\code{only.leaf = TRUE}), the descendant leaves of the given nodes will
#'     be found out and the FDR will be estimated on the leaf level.
#' @param only.leaf A logical scalar. If \code{TRUE}, the false discovery rate
#'     is calculated at the leaf (tip) level; otherwise it is calculated at
#'     the node level.
#'
#' @return The estimated false discovery rate.
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
#' ## FDR at the leaf level (1/8)
#' fdr(tree = tinyTree, truth = c(16, 13),
#'     found = c(15, 14), only.leaf = TRUE)
#'
#' ## FDR at the node level (2/14)
#' fdr(tree = tinyTree, truth = c(16, 13),
#'     found = c(15, 14), only.leaf = FALSE)
#'
fdr <- function(tree, truth, found, only.leaf = TRUE) {
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

    tt <- .fdr0(tree = tree, truth = truth,
                found = found, only.leaf = only.leaf)
    fdr <- tt[1]/tt[2]
    names(fdr) <- "fdr"

    return(fdr)
}


#' Count false positives and all positives
#'
#' Calculate the number of false discoveries and the total
#' number of discoveries on a tree structure, at leaf or node level.
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
#' @return A vector with two numbers, the number of false discoveries and the
#' total number of discoveries.
#'
#' @keywords internal
#' @noRd
#'
#' @importFrom TreeSummarizedExperiment findDescendant
#'
.fdr0 <- function(tree, truth = NULL, found = NULL, only.leaf = TRUE) {

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

    ## Count discoveries
    ## -------------------------------------------------------------------------
    if (is.null(found) || length(found) == 0) {
        ## Nothing is found
        c(fd = 0, disc = 1)
    } else {
        ## Something is found
        nodeF <- TreeSummarizedExperiment::findDescendant(
            tree = tree, node = found, only.leaf = only.leaf,
            self.include = TRUE)
        nodeF <- unique(unlist(nodeF))

        if (is.null(truth)) {
            ## No truly differential nodes
            c(fd = length(nodeF), disc = length(nodeF))
        } else {
            ## Some truly differential nodes
            nodeT <- TreeSummarizedExperiment::findDescendant(
                tree = tree, node = truth, only.leaf = only.leaf,
                self.include = TRUE)
            nodeT <- unique(unlist(nodeT))

            fd <- setdiff(nodeF, nodeT)
            c(fd = length(fd), disc = length(nodeF))
        }
    }
}


