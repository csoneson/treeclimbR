#' test whether nodes are connected
#' 
#' \code{isConnect} tests whether two (vectors of) nodes are connected.
#' 
#' @param tree A phylo object
#' @param node_a A vector of nodes.
#' @param node_b A vector of nodes.
#' @param connect It is choosen among "any", "direct", "indirect" and defines
#'   how two (vectors of ) nodes are connected.
#'   
#' @export
#' @importFrom methods is
#' @author Ruizhu Huang
#' @examples 
#'
#' data("tinyTree")
#' node_a <- c(4, 18, 19, 2)
#' node_b <- c(4, 5, 6, 3)
#' 
#' isConnect(tree = tinyTree, node_a = node_a,
#'           node_b = node_b, connect = "any")
#' 
#' 
isConnect <- function(tree, node_a, node_b, connect = "any") {
    
    if (!is(tree, "phylo")) {
        stop("tree should be a phylo object")
    }
    
    if (is.character(node_a)) {
        node_a <- transNode(tree = tree, node = node_a, use.alias = FALSE)
    }
    
    if (is.character(node_b)) {
        node_b <- transNode(tree = tree, node = node_b, use.alias = FALSE)
    }
    # path
    path <- matTree(tree = tree)
    
    # Which paths are node_a in
    # path_v <- as.vector(path) 
    num_a <- lapply(node_a, FUN = function(x) {
        which(path == x, arr.ind = TRUE)
    })
    
    num_b <- lapply(node_b, FUN = function(x) {
        which(path == x, arr.ind = TRUE)
    })
    
    
    num_r <- mapply(function(x, y) {
        length(intersect(x[, "row"], y[, "row"])) > 0
    }, num_a, num_b)
    num_neib <- mapply(function(x, y) {
        xy <- abs(outer(x[, "col"], y[, "col"], "-"))
        any(xy == 1)
    }, num_a, num_b)
    
    if (connect == "direct") {
        out <- num_r & num_neib
    }
    
    if (connect == "any") {
        out <- num_r
    }
    
    if (connect == "indirect") {
        out <- num_r & !num_neib
    }
    
    return(out)
    
}
