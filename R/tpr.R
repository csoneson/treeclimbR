#' Calculate true positive rate (TPR) on a tree structure
#'
#' \code{tpr} calculates the true positive rate (TPR) on a tree structure at
#' leaf or node level.
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
#' library(TreeSummarizedExperiment)
#' data("tinyTree")
#' ggtree(tinyTree) + 
#'    geom_text2(aes(label = node)) + 
#'    geom_hilight(node = 16, fill = "orange", alpha = 0.3) +
#'    geom_hilight(node = 13, fill = "blue", alpha = 0.3)
#'
#' # Truth: two branches have differential
#' # abundance under different conditions.
#'
#' # Found: branches with node 17 and node 14
#'  # TPR at the tip level
#' tpr1 <- tpr(tree = tinyTree, truth = c(16, 13),
#'             found = c(17, 14), only.leaf = TRUE)
#'  # TPR at the node level
#' tpr2 <- tpr(tree = tinyTree, truth = c(16, 13),
#'             found = c(15, 14), only.leaf = FALSE)
#'
#'


tpr <- function(tree, truth, found,
                only.leaf = TRUE) {
    
    if (!inherits(tree, "phylo")) {
        stop("tree: should be a phylo object")
    }
    
   if (is.character(truth)) {
        truth <- transNode(tree = tree, node = truth,
                           message = FALSE)
    }
    if (is.character(found)) {
        found <- transNode(tree = tree, node = found,
                           message = FALSE)
    }
    
    tt <- .tpr0(tree = tree, truth = truth,
                found = found, only.leaf = only.leaf)
    tpr <- tt[1]/tt[2]
    
    
    # return final results
    names(tpr) <- "tpr"
    return(tpr)
}


#' Calculate true positives and positives
#'
#' \code{.tpr0} calculates the number of true positives and the total number of
#' positives on a tree structure at leaf or node level.
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


.tpr0 <- function(tree,
                  truth = NULL,
                  found = NULL,
                  only.leaf = TRUE) {
    
    # ================= check inputs ===========================
    if (!is.null(truth)) {
        if (!(is.character(truth) |
              is.numeric(truth) |
              is.integer(truth))) {
            stop("truth should include character or numeric")
        }
    }
    
    if (!is.null(found)) {
        if (!(is.character(found) |
              is.numeric(found) |
              is.integer(found))) {
            stop("found should include character or numeric")
        } 
    }
    
    # ================= without diff ===========================
    if (is.null(truth) | length(truth) == 0) {
        c(tp = 1, pos = 1)
    } else {
        # ================= with diff ===========================
        nodeT <- findOS(tree = tree, node = truth,
                        only.leaf = only.leaf, self.include = TRUE)
        nodeT <- unique(unlist(nodeT))
        
        # no discovery
        if (is.null(found)) {
            c(tp = 0, pos = length(nodeT)) 
        } else {
            # found 
            nodeF <- findOS(tree = tree, node = found,
                            only.leaf = only.leaf, self.include = TRUE)
            nodeF <- unique(unlist(nodeF))
            # true positive & positive
            TP <- intersect(nodeT, nodeF)
            c(tp = length(TP), pos = length(nodeT))
        }
    }
}
