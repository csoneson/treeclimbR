#' Find branches that are non-overlapping with specified branches in a tree
#'
#' Find all branches whose leaves do not overlap with those of the specified
#' branches.
#'
#' @author Ruizhu Huang
#' @export
#'
#' @param tree A \code{phylo} object.
#' @param node A numeric or character vector specifying the nodes.
#' @param use.alias A logical scalar. If \code{TRUE}, the alias name is
#'     used to name the output vector.
#'
#' @returns A vector of node numbers
#'
#' @importFrom TreeSummarizedExperiment convertNode
#'
#' @examples
#' suppressPackageStartupMessages({
#'     library(ggtree)
#'     library(TreeSummarizedExperiment)
#' })
#'
#' data(tinyTree)
#' ggtree(tinyTree, branch.length = "none") +
#'     geom_text2(aes(label = node)) +
#'     geom_hilight(node = 17, fill = "blue", alpha = 0.3) +
#'     geom_hilight(node = 13, fill = "orange", alpha = 0.3)
#'
#' ## Find branches whose leaves do not overlap with the two colored branches.
#' ## The returned branches are represented at the highest tree level
#' ## possible without including any of the forbidden branches.
#' findExcl(tree = tinyTree, node = c(17, 13))
#'
findExcl <- function(tree, node, use.alias = FALSE) {
    ## Check input arguments
    ## -------------------------------------------------------------------------
    .assertVector(x = tree, type = "phylo")
    .assertScalar(x = use.alias, type = "logical")
    if (!(is.character(node) || is.numeric(node))) {
        stop("'node' should be either a character vector or a numeric ",
             "vector")
    }

    ## Convert node labels to node number
    ## -------------------------------------------------------------------------
    if (is.character(node)) {
        node <- TreeSummarizedExperiment::convertNode(
            tree = tree, node = node, message = FALSE)
    }

    ## Find descendant leaves of nodes
    ## -------------------------------------------------------------------------
    inT <- TreeSummarizedExperiment::findDescendant(
        tree = tree, node = node, only.leaf = TRUE)
    inT <- unlist(inT)

    ## Get all leaves in the tree
    ## -------------------------------------------------------------------------
    allT <- setdiff(tree$edge[, 2], tree$edge[, 1])

    ## Leaves that are not descendants of node
    ## -------------------------------------------------------------------------
    exT <- setdiff(allT, inT)

    ## Replace leaves with their ancestor branch node
    ## -------------------------------------------------------------------------
    out <- TreeSummarizedExperiment::joinNode(tree = tree, node = exT)

    ## Return a vector of the found node (the node number of the node)
    ## Name the vector with the node label
    ## -------------------------------------------------------------------------
    names(out) <- TreeSummarizedExperiment::convertNode(
        tree = tree, node = out, use.alias = use.alias, message = FALSE)

    return(out)
}
