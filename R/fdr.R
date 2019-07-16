#' Calculate false discovery rate (FDR) on a tree structure
#'
#' \code{FDR} calculates the false discovery rate on a tree structure at leaf or
#' node level.
#'
#' @param tree A phylo object
#' @param truth Nodes that have signals (eg. differentally abundant at different
#'   experimental conditions.). \strong{Note:} when TPR is required at leaf
#'   level (\code{only.leaf = TRUE}), the descendant leaves of the given nodes
#'   will be found out and the TPR is calculated on the leaf level;
#' @param found Nodes that have been found to have signal (eg. differentally
#'   abundant at different experimental conditions). \strong{Note:} when TPR is
#'   required at leaf level (\code{only.leaf = TRUE}), the descendant leaves of
#'   the given nodes will be found out and the TPR is calculated on the leaf
#'   level;
#' @param only.leaf A logical value, TRUE or FALSE. If TRUE, true positive rate
#'   is calculated at the leaf (tip) level; otherwise it is calculated at node
#'   level. The default is TRUE
#'   
#' @export
#' @return A true positive rate
#' @author Ruizhu Huang
#' @examples
#' library(ggtree)
#' data("tinyTree")
#' ggtree(tinyTree) + 
#'    geom_text2(aes(label = node)) + 
#'    geom_hilight(node = 16, fill = "orange", alpha = 0.3) +
#'    geom_hilight(node = 13, fill = "blue", alpha = 0.3)
#'
#' # Truth: two branches have differential
#' # abundance under different conditions.
#'
#' fdr1 <- fdr(tree = tinyTree, truth = c(16, 13),
#'             found = c(17, 14), only.leaf = TRUE)
#'  # fdr at the node level
#' fdr2 <- fdr(tree = tinyTree, truth = c(16, 13),
#'             found = c(15, 14), only.leaf = FALSE)
#'
#'


fdr <- function(tree, truth, found,
                only.leaf = TRUE) {
    
    if (!inherits(tree, "phylo")) {
        stop("tree: should be a phylo object")
    }
    
    if (is.character(truth)) {
        truth <- transNode(tree = tree, input = truth,
                           message = FALSE)
    }
    if (is.character(found)) {
        found <- transNode(tree = tree, input = found,
                           message = FALSE)
    }
    
    tt <- .fdr0(tree = tree, truth = truth,
                found = found, only.leaf = only.leaf)
    fdr <- tt[1]/tt[2]
    
    
    # return final results
    names(fdr) <- "fdr"
    return(fdr)
}


#' Calculate true positives and positives
#'
#' \code{.fdr0} calculates the number of false discovery and the total number of
#' discovery on a tree structure at leaf or node level.
#'
#' @param tree A phylo object
#' @param truth Nodes that have signals (eg. differentally abundant at different
#'  experimental conditions.).
#' @param found Nodes that have been found to have signal
#' @param only.leaf A logical value, TRUE or FALSE. If TRUE, true positive rate
#'   is calculated at the leaf (tip) level; otherwise it is calculated at node
#'   level. The default is TRUE
#' @return a vector
#' @keywords internal


.fdr0 <- function(tree,
                  truth = NULL,
                  found = NULL,
                  only.leaf = TRUE) {
    
    # ================= check inputs =======================
    if (!is.null(truth)) {
        # if diff exists, check whether the truth has correct input format
        if (!(is.character(truth) |
              is.numeric(truth) |
              is.integer(truth))) {
            stop("truth should include character or numeric")
        }
    }
    
    if (!is.null(found)) {
        # check whether 'found' has correct input
        if (!(is.character(found) |
              is.numeric(found) |
              is.integer(found))) {
            stop("found should include character or numeric")
        } 
    }
    
    # ================= without discovery =======================
    if (is.null(found)) {
        c(fd = 0,  disc = 0)
    } else {
        
        # =================   with discovery  =======================
        nodeF <- findOS(tree = tree, node = found,
                        only.leaf = only.leaf, self.include = TRUE)
        nodeF <- unlist(nodeF)
        
        # if no diff
        if (is.null(truth)) {
            c(fd = length(nodeF), disc = length(nodeF)) 
        } else {
           # if has diff
           nodeT <- findOS(tree = tree, node = truth,
                           only.leaf = only.leaf, self.include = TRUE)
           nodeT <- unlist(nodeT)
           # false discovery & discovery
           fd <- setdiff(nodeF, nodeT)
           c(fd = length(fd), disc = length(nodeF))
        }
    }
}


