#' Check whether nodes are contained in the same path from a leaf to the root
#' in a tree
#'
#' Perform an elementwise check of whether two vectors of nodes are "connected"
#' in specific ways in a tree. A pair of nodes are considered to be connected
#' if they are part of the same path from a leaf to the root of the tree.
#' They are considered directly connected if one node is the parent of the
#' other, and indirectly connected otherwise.
#'
#' @author Ruizhu Huang, Charlotte Soneson
#' @export
#'
#' @param tree A \code{phylo} object.
#' @param node_a,node_b The two vectors of nodes (either node numbers or
#'     node labels) to check for connections. The vectors should have the same
#'     length (if not, the shorter one will be recycled), as the check for
#'     connectivity is done elementwise.
#' @param connect One of "any", "direct", "indirect", the type of connections
#'     to search for.
#'
#' @returns A logical vector of the same length as \code{node_a} and
#' \code{node_b}, where each element indicates whether the corresponding
#' elements of \code{node_a} and \code{node_b} are connected in the
#' specified way.
#'
#' @importFrom TreeSummarizedExperiment convertNode matTree
#'
#' @examples
#' suppressPackageStartupMessages({
#'     library(ggtree)
#' })
#'
#' data(tinyTree)
#' ggtree(tinyTree, branch.length = "none") +
#'     geom_text2(aes(label = node))
#'
#' node_a <- c(4, 18, 19, 2)
#' node_b <- c(4, 5, 6, 3)
#'
#' isConnect(tree = tinyTree, node_a = node_a,
#'           node_b = node_b, connect = "any")
#'
isConnect <- function(tree, node_a, node_b, connect = "any") {
    ## Check input arguments
    ## -------------------------------------------------------------------------
    .assertVector(x = tree, type = "phylo")
    if (!(is.character(node_a) || is.numeric(node_a))) {
        stop("'node_a' should be either a character vector or a numeric ",
             "vector")
    }
    if (!(is.character(node_b) || is.numeric(node_b))) {
        stop("'node_b' should be either a character vector or a numeric ",
             "vector")
    }
    .assertScalar(x = connect, type = "character",
                  validValues = c("any", "direct", "indirect"))

    ## Convert node labels to node numbers
    ## -------------------------------------------------------------------------
    if (is.character(node_a)) {
        node_a <- TreeSummarizedExperiment::convertNode(
            tree = tree, node = node_a, use.alias = FALSE)
    }
    if (is.character(node_b)) {
        node_b <- TreeSummarizedExperiment::convertNode(
            tree = tree, node = node_b, use.alias = FALSE)
    }

    ## Get paths from each leaf to the root
    ## -------------------------------------------------------------------------
    path <- TreeSummarizedExperiment::matTree(tree = tree)

    ## Find paths that node_a and node_b are in, respectively
    ## -------------------------------------------------------------------------
    num_a <- lapply(node_a, FUN = function(x) {
        which(path == x, arr.ind = TRUE)
    })
    num_b <- lapply(node_b, FUN = function(x) {
        which(path == x, arr.ind = TRUE)
    })

    ## Check for overlap
    ## -------------------------------------------------------------------------
    num_r <- mapply(function(x, y) {
        length(intersect(x[, "row"], y[, "row"])) > 0
    }, num_a, num_b)
    num_neib <- mapply(function(x, y) {
        bdiff <- abs(outer(x[, "row"], y[, "row"], "-"))
        hdiff <- abs(outer(x[, "col"], y[, "col"], "-"))
        any(hdiff == 1 & bdiff == 0)
    }, num_a, num_b)

    if (connect == "direct") {
        out <- num_r & num_neib
    } else if (connect == "any") {
        out <- num_r
    } else if (connect == "indirect") {
        out <- num_r & !num_neib
    }

    return(out)
}
