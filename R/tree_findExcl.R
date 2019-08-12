#' find branches that are exclusive from the specified branches
#' 
#' \code{findExcl} is to find all branches that have leaves exclusive from the
#' specified branches. In the same time, these branches group the leaves to have
#' the min number of branches. Examples are given in example section to better
#' clarify this.
#' 
#' @param tree A phylo object.
#' @param node A numeric or charater vector specifying the nodes.
#' @param use.alias A logical value, TRUE or FALSE. If TRUE, the alias name is
#'   used to name the output vector.
#' @export
#' @examples 
#' library(ggtree)
#' data(tinyTree)
#' ggtree(tinyTree, branch.length = "none") +
#' geom_text2(aes(label = node)) +
#' geom_hilight(node = 17, fill = "blue", alpha = 0.3) +
#' geom_hilight(node = 13, fill = "orange", alpha = 0.3)
#' 
#' # To output branches that have leaves not in the two colored branches, the
#' # shortest result here is three. They can't be grouped as less branches
#' # without including the leaves in the specified branches.
#' findExcl(tree = tinyTree, node = c(17, 13))


findExcl <- function(tree, node,
                     use.alias = FALSE){
    
    if (is.character(node)) {
        node <- transNode(tree = tree, node = node,
                           message = FALSE)
    } else {
        node <- node
    }
    
    # find descendant leaves of node
    inT <- findOS(tree = tree, node = node, only.leaf = TRUE)
    inT <- unlist(inT)
    # find all leaves of tree
    allT <- setdiff(tree$edge[,2], tree$edge[, 1])
    
    # Leaves not included in node
    exT <- setdiff(allT, inT)
    
    # replace leaves with their ancestor branch node
    out <- signalNode(tree = tree, node = exT)
    
    # return a vector of the found node (the node number of the node)
    # name the vector with the node label
    names(out) <- transNode(tree = tree, node = out,
                            use.alias = use.alias,
                            message = FALSE)
    return(out)
    
}
