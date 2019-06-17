#' Find the children
#'
#' \code{findChild} finds children of an internal node.
#'
#' @param node An internal node. It could be the node number or the node
#'   label.
#' @param tree A phylo object.
#' @param use.alias A logical value, TRUE or FALSE. The default is FALSE, and
#'   the node label would be used to name the output; otherwise, the alias of
#'   node label would be used to name the output. The alias of node label is
#'   created by adding a prefix \code{"alia_"} to the node number.
#' @export
#' @return A vector of nodes. The numeric value is the node number, and the
#'   vector name is the corresponding node label. If a node has no label, it
#'   would have NA as name when \code{use.alias = FALSE}, and have the alias of
#'   node label as name when \code{use.alias = TRUE}.
#' @author Ruizhu Huang
#'
#' @examples
#' data(tinyTree)
#'
#' library(ggtree)
#' ggtree(tinyTree) +
#' geom_text2(aes(label = node), color = "darkblue",
#'                hjust = -0.5, vjust = 0.7) +
#' geom_hilight(node = 17, fill = 'steelblue', alpha = 0.5) +
#' geom_text2(aes(label = label), color = "darkorange",
#'            hjust = -0.1, vjust = -0.7)
#'
#' (tips <- findChild(tree = tinyTree, node = 17))
#' 



findChild <- function(tree, node = 11, use.alias = FALSE) {
    
    if (!is(tree, "phylo")) {
        stop("tree requires a phylo object. \n")
    }
    
    if (is.factor(node)) {
        stop("node requires character or numeric.")
    }
    
    if (is.character(node)) {
        node <- transNode(tree = tree, node = node)
    }
    
    mat <- tree$edge
    ind <- which(mat == node, arr.ind = TRUE)
    ind <- ind[ind[,"col"] == 1, , drop = FALSE]
    
    ind[, "col"] <- ind[, "col"] + 1
    child <- mat[ind]
    names(child) <- transNode(tree = tree, node = child,
                              use.alias = use.alias, message = FALSE)
    return(child)
}
