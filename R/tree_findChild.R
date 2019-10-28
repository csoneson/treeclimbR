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
#' (tips <- findChild(tree = tinyTree, 
#'                    node = c(17, 12)))
#' 



findChild <- function(tree, node,
                      use.alias = FALSE) {
    
    if (!inherits(tree, "phylo")) {
        stop("tree: should be a phylo object")
    }
    
    if (!(is.character(node) |
          is.numeric(node) |
          is.integer(node))) {
        stop("The argument (node) should be character or numeric")
    }
    # the edge matrix
    mat <- tree$edge
    matN <- matTree(tree = tree)
    
    if (is.character(node)) {
        numA <- transNode(tree = tree, node = node,
                          use.alias = TRUE,
                          message = FALSE)
    } else {
        numA <- node
        isOut <- !numA %in% mat
        if (any(isOut)) {
            stop("Node ", numA,
                 " can't be found in the ",
                 deparse(substitute(tree)), "\n")
        }
        
    }
    
    # find children
    loc1 <- lapply(numA, FUN = function(x) {
        xi <- which(matN == x, arr.ind = TRUE)
        return(xi)
    })
    
    loc2 <- lapply(loc1, FUN = function(x) {
        x1 <- x[, "row", drop = FALSE]
        x2 <- x[, "col", drop = FALSE] - 1
        cbind(row = x1, col = x2)
        
    })
    
    matNN <- apply(matN, 2, FUN = function(x) {
        xe <- x[!is.na(x)]
        xx <- transNode(tree = tree, node = xe,
                        use.alias = use.alias,
                        message = FALSE)
        x[!is.na(x)] <- xx
        return(x)
    })
    
    
    # descendants: tips or leaves
    tipA <- unique(setdiff(mat[, 2], mat[, 1]))
    
    chl <- lapply(seq_along(loc2), FUN = function(x) {
        xx <- loc2[[x]]
        y <- matN[xx]
        uy <- unique(y)
        return(uy)
    })
    
    # final output (node number or label)
    names(chl) <- transNode(tree = tree, node = numA,
                            use.alias = use.alias,
                            message = FALSE)
    return(chl)
}
