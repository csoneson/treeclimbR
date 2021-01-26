#' Select branches
#'
#' \code{selNode} selects branches meeting the specified criteria in the number
#' of leaves and the count proportion.
#'
#' @param pr A named numeric vector to provide proportions of entities.
#' @param obj A TreeSummarizedExperiment object.
#' @param tree A phylo object
#' @param data A count table with entities in rows and samples in columns.
#' @param minTip the minimum number of leaves in the selected branch
#' @param maxTip The maximum number of leaves in the selected branch
#' @param minPr The minimum count proportion of the selected branch in a sample.
#'   A value between 0 and 1.
#' @param maxPr The maximum count proportion of the selected branch in a sample.
#'   A value between 0 and 1.
#' @param skip A character vector of node labels. These nodes are not the
#'   descendants or the ancestors of the selected branch.
#' @param all TRUE or FALSE. Default is FALSE. If FALSE, the branch node of a
#'   branch, which meet the requirements and has the minimum count proportion,
#'   is returned; otherwise branch nodes of all branches meeting the
#'   requirements are returned.
#'
#' @export
#' @importFrom methods is
#' @return The node whose descendant branch has the lowest proportion
#' @author Ruizhu Huang
#' @examples
#' set.seed(1)
#' data("tinyTree")
#' toyTable <- matrix(rnbinom(40, size = 1, mu = 10), nrow = 10)
#' colnames(toyTable) <- paste(rep(LETTERS[1:2], each = 2),
#' rep(1:2, 2), sep = "_")
#' rownames(toyTable) <- tinyTree$tip.label
#'
#'
#' dat <- parEstimate(obj = toyTable)
#' (out1 <- selNode(tree = tinyTree, data = dat, all = TRUE))
#' (out2 <- selNode(tree = tinyTree, data = dat, 
#'                  minTip = 4, maxTip = 9,
#'                  minPr = 0, maxPr = 0.8, all = TRUE))
#'
#' ## provide obj
#' lse <- TreeSummarizedExperiment(rowTree = tinyTree,
#'                                 assays = list(toyTable))
#' (out3 <- selNode(obj = lse, all = TRUE))
#' # All nodes except node 1
#' (out4 <- selNode(obj = lse, skip = 1, all = TRUE))
#'

selNode <- function(pr = NULL, obj = NULL, 
                    data = NULL, tree = NULL,
                    minTip = 0, maxTip = Inf,
                    minPr = 0, maxPr = 1,
                    skip = NULL, all = FALSE){
    # if proportions of entities on leaves are not available, a
    # TreeSummarizedExperiment object (obj) or data & tree shouls be provided.
    if (!is.null(pr)) {
        if (!is.numeric(pr)) {stop("pr should be a numeric vector.")}
        if (is.null(tree)) {stop("Tree is required.")}
        if (!is.null(obj) | !is.null(data)) {
            stop("The proportion is already given in 'pr'. 
                 obj and data should be NULL")
        }
        
        pars <- pr
    } else {
    
        # if the obj is provided with a leafSummarizedExperiment object
        # use it; otherwise use tree and data
        if (is(obj, "TreeSummarizedExperiment")) {
            if (any(!is.null(tree), !is.null(data))) {
                stop("Set tree and data as NULL when obj is given. \n")
            }
            obj <- parEstimate(obj = obj)
            tree <- rowTree(obj)
            data <- metadata(obj)$assays.par
            pars <- data$pi
        } else {
            tree <- tree
            data <- data
            pars <- parEstimate(obj = data)$pi
        }
        
    }
    ##------------ descendant tips ------------
    # descendant tips for each internal node
    
    # proportion of internal nodes
    leaf <- setdiff(tree$edge[, 2], tree$edge[, 1])
    nodI <- setdiff(tree$edge[, 1], leaf)
    # nodA <- c(leaf, nodI)
    nodeLab <- convertNode(tree = tree, node = nodI,
                         use.alias = FALSE,
                         message = FALSE)
    nodeLab_alias <- convertNode(tree = tree, node = nodI,
                               use.alias = TRUE,
                               message = FALSE)
    
    # descendant leaves
    tipI <- findDescendant(tree = tree, node = nodI,
                   only.leaf = TRUE, self.include = TRUE)
    names(tipI) <- nodeLab_alias
    
    # number of descendant leaves
    numI <- unlist(lapply(tipI, length))
    #numI[seq_along(leaf)] <- 0
    
    ##------------ proportions on nodes -----------
    # a vector of node numbers with node labels as names
    vnum <- convertNode(tree = tree, node = names(pars),
                      message = FALSE)
    
    # proportion for each node
    propList <- lapply(tipI, FUN = function(x){
        sum(pars[match(x, vnum)])
    })
    nodP <- unlist(propList)
    
    ##---------- sample ---------------
    if (any(duplicated(nodeLab))) {
        tt <- cbind.data.frame(nodeNum = convertNode(tree = tree,
                                                   node = names(nodP),
                                                   message = FALSE),
                               nodeLab = nodeLab,
                               nodeLab_alias = nodeLab_alias,
                               proportion = nodP,
                               numTip = numI,
                               stringsAsFactors =FALSE)
    } else {
        tt <- cbind.data.frame(nodeNum = convertNode(tree = tree,
                                                   node = names(nodP),
                                                   message = FALSE),
                               nodeLab = nodeLab,
                               proportion = nodP,
                               numTip = numI,
                               stringsAsFactors =FALSE)
    }
    
    
    if (maxPr < min(tt$proportion)) {
        stop("maxPr defined is even lower than the minimum value of
             node proportion", signif(min(tt$proportion),2), "\n")
    }
    # only consider nodes with enough tips and
    # desired proportion level
    st <- tt[tt$numTip >= minTip &
                 tt$numTip <= maxTip &
                 tt$proportion >= minPr &
                 tt$proportion <= maxPr,]
    if (nrow(st) == 0) {
        stop("No nodes fullfill the requirements;
             try other settings
             for tip numbers or proportions")
    }
    # remove those overlapped
    if (!is.null(skip)) {
        if (is.character(skip)) {
            skip <- convertNode(tree = tree, node = skip,
                              message = FALSE)
        }
        tipS <- findDescendant(tree = tree, node = skip,
                       only.leaf = TRUE, self.include = TRUE)
        tipS <- unlist(tipS)
        
        rmp <- vapply(st$nodeNum, FUN = function(x){
            tx <- findDescendant(node = x, tree = tree, only.leaf = TRUE,
                         self.include = TRUE)
            ix <- intersect(tipS, unlist(tx))
            length(ix) == 0
        }, FUN.VALUE = TRUE)
        
        new.st <- st[rmp, ]
    } else {
        new.st <- st
    }
    
    if (nrow(new.st) == 0) {
        stop("No nodes fullfill the requirements;
             try other settings
             for tip numbers or proportions")
    }
    
    # return the one has the lowest proportion if all = FALSE
    if (all) {
        final <- new.st
    } else {
        final <- new.st[which.min(new.st$proportion), ]
    }
    
    return(final)
    
}
