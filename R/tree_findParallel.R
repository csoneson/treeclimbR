findParallel <- function(tree, node,
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
