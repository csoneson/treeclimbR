#' subset the TreeSummarizedExperiment object 
#'
#' \code{subsetTse} is to subset the TreeSummarizedExperiment object. The
#' corresponding tree could be optionally updated. If the tree is allowed to be
#' updated, branches with no data in the \code{assays} table will be removed to
#' create a new tree structure. Users might find leaf nodes that have no
#' corresponding rows (or columns) in the \code{assays} table are kept in the
#' \code{rowTree} (or \code{colTree}). That's because at least one of its
#' ancestor nodes (nodes in the path connecting the leaf and the root) could be
#' mapped to rows (or columns) in the \code{assays} table.
#'
#' @param x A TreeSummarizedExperiment object
#' @param subset_rows A logical expression indicating rows to keep.
#' @param subset_columns A logical expression indicating columns to keep
#' @param update_rowTree A logical value, TRUE or FALSE. If TRUE, the row tree
#'   is updated correspondingly.
#' @param update_colTree A logical value, TRUE or FALSE. If TRUE, the column
#'   tree is updated correspondingly.
#' @param trim_internal_rowTree A logical value, TRUE or FALSE. If TRUE,
#'   internal nodes with chldren removed are dropped from the
#'   \code{rowTree}.
#' @param trim_internal_colTree A logical value, TRUE or FALSE. If TRUE,
#'   internal nodes with chldren removed are dropped from the
#'   \code{colTree}.
#' @param collapse_signles_rowTree A logical value, TRUE or FALSE. If TRUE,
#'   internal nodes with only one direct child node is removed from the
#'   \code{rowTree}.
#' @param collapse_signles_colTree A logical value, TRUE or FALSE. If TRUE,
#'   internal nodes with only one direct child node is removed from the
#'   \code{colTree}.
#' @import TreeSummarizedExperiment
#' @importFrom ape drop.tip
#' @importFrom methods as
#' @export
#' @return A TreeSummarizedExperiment object
#' @author Ruizhu Huang
#' @examples 
#' 
#' data("tinyTree")
#' set.seed(1)
#' # the count table
#' count <- matrix(rpois(100, 50), nrow = 10)
#' rownames(count) <- c(tinyTree$tip.label)
#' colnames(count) <- paste("C_", 1:10, sep = "_")
#'
#' # The sample information
#' sampC <- data.frame(condition = rep(c("control", "trt"), each = 5),
#'                     gender = sample(x = 1:2, size = 10, replace = TRUE))
#' rownames(sampC) <- colnames(count)
#'
#' # build a TreeSummarizedExperiment object
#' tse <- TreeSummarizedExperiment(assays = list(count),
#'                                 colData = sampC,
#'                                 rowTree = tinyTree)
#' rowLinks(tse)
#'                        
#' sel <- sample(sample(c(TRUE, FALSE), size = 10, replace = TRUE))         
#' tse1 <- subsetTse(x = tse, subset_rows = sel, update_rowTree = TRUE )    
#' rowLinks(tse1) 
subsetTse <- function(x,
                     subset_rows,
                     subset_columns,
                     update_rowTree = FALSE,
                     update_colTree = FALSE,
                     trim_internal_rowTree = TRUE,
                     trim_internal_colTree = TRUE,
                     collapse_singles_rowTree = FALSE,
                     collapse_singles_colTree = FALSE,
                     message = FALSE) {
    # ==================== check inputs ==========================
    # The input data should be a TreeSummarizedExperiment object
    if (!is(x, "TreeSummarizedExperiment")) {
        stop("x should be a TreeSummarizedExperiment.")
    }
    
    
    if (!missing(subset_rows)) {
        if (!is.logical(subset_rows)) {
            stop("subset_rows should be a logical vector")
        }
        if (length(subset_rows) != nrow(x)) {
            stop("subset_rows should have length equal to ", nrow(x))
        }
    }   
    
    if (!missing(subset_columns)) {
        if (!is.logical(subset_columns)) {
            stop("subset_columns should be a logical vector")
        }
        if (length(subset_columns) != ncol(x)) {
            stop("subset_columns should have length equal to ", ncol(x))
        }
    }
    
    # ==================== subset TSE ===========================
    if (missing(subset_rows)) {
        subset_rows <- rep(TRUE, nrow(x))
    }
    if (missing(subset_columns)) {
        subset_columns <- rep(TRUE, ncol(x))
    }
    nx <- x[subset_rows, subset_columns]
    
    # ==================== update trees ==========================
    if (update_rowTree) {
        if (message) {
            message("working on the row tree...")
        }
        rowLN <- DataFrame(rowLinks(nx))
        treeR <- rowTree(nx)
        
        if (is.null(treeR)) {
            stop("rowTree is not available.")
        }
        # leaves to be removed
        allLeaf <- printNode(tree = treeR, type = "leaf")$nodeNum
        existNode <- rowLN$nodeNum
        existLeaf <- unionLeaf(tree = treeR, node = existNode)
        rnode <- setdiff(allLeaf, unlist(existLeaf))
        
        
        track <- trackNode(tree = treeR)
        
        if (length(rnode)) {
        # track <- pruneTree(tree = track, rmLeaf = rnode,
        #                    mergeSingle = merge_rowTree)
        track <- drop.tip(phy = track, tip = rnode, 
                          trim.internal = trim_internal_rowTree,
                          collapse.singles = collapse_singles_rowTree)
        treeR <- drop.tip(phy = treeR, tip = rnode, 
                          trim.internal = trim_internal_rowTree,
                          collapse.singles = collapse_singles_rowTree)
        
        }
        
        if (!is.null(track)) {
        # the new row link data
        rowLN$nodeNum <- transNode(tree = track, node = rowLN$nodeLab_alias)
        rowLN$nodeLab <- transNode(tree = treeR, node = rowLN$nodeNum, 
                                   use.alias = FALSE)
        rowLN$nodeLab_alias <- transNode(tree = treeR, node = rowLN$nodeNum,
                                         use.alias = TRUE)
        rowLN <- as(rowLN, "LinkDataFrame")
        } else {
            rowLN <- NULL
        }
       
        # update the row tree and likns
        nx <- BiocGenerics:::replaceSlots(nx, rowLinks = rowLN, 
                                         rowTree = list(phylo = treeR))
    }
    
    if (update_colTree) {
        if (message) {
            message("working on the row tree...")
        }
        colLN <- DataFrame(colLinks(nx))
        treeC <- colTree(nx)
        if (is.null(treeC)) {
            stop("colTree is not available.")
        }
        # leaves to be removed
        allLeaf <- printNode(tree = treeC, type = "leaf")$nodeNum
        existNode <- colLN$nodeNum
        existLeaf <- unionLeaf(tree = treeC, node = existNode)
        rnode <- setdiff(allLeaf, unlist(existLeaf))
        
        
        track <- trackNode(tree = treeC)
        
        if (length(rnode)) {
            track <- drop.tip(phy = track, tip = rnode, 
                              trim.internal = trim_internal_colTree,
                              collapse.singles = collapse_singles_colTree)
            treeC <- drop.tip(phy = treeC, tip = rnode, 
                              trim.internal = trim_internal_colTree,
                              collapse.singles = collapse_singles_colTree)
            
        }
        
        # the new col link data
        colLN$nodeNum <- transNode(tree = track, node = colLN$nodeLab_alias)
        colLN$nodeLab <- transNode(tree = treeC, node = colLN$nodeNum, 
                                   use.alias = FALSE)
        colLN$nodeLab_alias <- transNode(tree = treeC, node = colLN$nodeNum,
                                         use.alias = TRUE)
        colLN <- as(colLN, "LinkDataFrame")
        
        # update the col tree and likns
        nx <- BiocGenerics:::replaceSlots(nx, colLinks = colLN, 
                                          colTree = list(phylo = treeC))
    }
    
    return(nx)
}
