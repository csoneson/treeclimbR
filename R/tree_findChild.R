#' Find the children of an internal node in a tree
#'
#' Find the direct children of an internal node in a tree.
#'
#' @author Ruizhu Huang, Charlotte Soneson
#' @export
#'
#' @param node Either the node number or node label of an internal node of
#'     \code{tree}.
#' @param tree A \code{phylo} object.
#' @param use.alias A logical scalar. If \code{FALSE} (default), the node label
#'     is used to name the output; otherwise, the alias of the node label is
#'     used. The alias of the node label is created by adding a prefix
#'     \code{"alias_"} to the node number.
#'
#' @return A vector of nodes. The numeric value is the node number, and the
#'     vector name is the corresponding node label. If a node has no label, it
#'     would have NA as name when \code{use.alias = FALSE}, and have the alias
#'     of the node label as name when \code{use.alias = TRUE}.
#'
#' @importFrom TreeSummarizedExperiment matTree convertNode
#'
#' @examples
#' library(ggtree)
#'
#' data(tinyTree)
#' ggtree(tinyTree) +
#'     geom_text2(aes(label = node), color = "darkblue",
#'                hjust = -0.5, vjust = 0.7) +
#'     geom_hilight(node = 17, fill = 'steelblue', alpha = 0.5) +
#'     geom_text2(aes(label = label), color = "darkorange",
#'                hjust = -0.1, vjust = -0.7)
#'
#' ## Specify node numbers
#' findChild(tree = tinyTree, node = c(17, 12))
#'
#' ## Name return values using aliases
#' findChild(tree = tinyTree, node = c(17, 12), use.alias = TRUE)
#'
#' ## Specify node labels
#' findChild(tree = tinyTree, node = c("Node_17", "Node_12"))
#'
#' ## Tips have no children
#' findChild(tree = tinyTree, node = "t4")
#'
findChild <- function(tree, node, use.alias = FALSE) {
    ## Check input arguments
    ## -------------------------------------------------------------------------
    .assertVector(x = tree, type = "phylo")
    .assertScalar(x = use.alias, type = "logical")
    if (!(is.character(node) || is.numeric(node))) {
        stop("'node' should be either a character vector or a numeric ",
             "vector")
    }

    ## Get the edge matrix and the path from each node to the root
    ## -------------------------------------------------------------------------
    mat <- tree$edge
    matN <- TreeSummarizedExperiment::matTree(tree = tree)

    ## Get node number and check that it is part of the tree
    ## -------------------------------------------------------------------------
    if (is.character(node)) {
        numA <- TreeSummarizedExperiment::convertNode(
            tree = tree, node = node, use.alias = TRUE, message = FALSE)
    } else {
        numA <- node
        isOut <- !numA %in% mat
        if (any(isOut)) {
            stop("Node ", numA, " can't be found in ",
                 deparse(substitute(tree)), "\n")
        }
    }

    ## Find positions of each node in matN
    ## -------------------------------------------------------------------------
    loc1 <- lapply(numA, FUN = function(x) {
        which(matN == x, arr.ind = TRUE)
    })

    ## Find positions of children in matN
    ## -------------------------------------------------------------------------
    loc2 <- lapply(loc1, FUN = function(x) {
        x1 <- x[, "row", drop = FALSE]
        x2 <- x[, "col", drop = FALSE] - 1
        cbind(row = x1, col = x2)
    })

    ## Get node numbers for children
    ## -------------------------------------------------------------------------
    chl <- lapply(seq_along(loc2), FUN = function(x) {
        xx <- loc2[[x]]
        y <- matN[xx]
        return(unique(y))
    })

    ## Final output (node number or label)
    ## -------------------------------------------------------------------------
    names(chl) <- TreeSummarizedExperiment::convertNode(
        tree = tree, node = numA, use.alias = use.alias, message = FALSE)

    return(chl)
}
