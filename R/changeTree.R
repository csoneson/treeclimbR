#' change the tree in TreeSummarizedExperiment
#' 
#' \code{changeTree} changes a row or column tree to a new one by matching the
#' node labels.
#' 
#' @param tse A TreeSummarizedExperiment object.
#' @param new_colTree A phylo object. A new tree to replace the current column
#'   tree.
#' @param new_rowTree A phylo object. A new tree to replace the current row
#'   tree.
#' @param lab_colTree A data frame. A data frame with two columns (oldLab and
#'   newLab), one having the node labels of the current column tree, and the
#'   other having those of the new tree.
#' @param lab_rowTree A data frame. A data frame with two columns (oldLab and
#'   newLab), one having the node labels of the current row tree, and the other
#'   having those of the new tree.
#'   
#' @importFrom methods is
#' @export
#' @return A TreeSummarizedExperiment object

changeTree <- function(tse, 
                       new_colTree = NULL,
                       new_rowTree = NULL,
                       lab_colTree = NULL, 
                       lab_rowTree = NULL) {
    # check TreeSummarizedExperiment object
    if (!is(tse, "TreeSummarizedExperiment")) {
        stop("tse should be a TreeSummarizedExperiment object")
    }
    
    # check trees
    if (!is.null(new_colTree) & !is(new_colTree, "phylo")) {
        stop("new_colTree should be a phylo object")
    }
    
    if (!is.null(new_rowTree) & !is(new_rowTree, "phylo")) {
        stop("new_rowTree should be a phylo object")
    }
    
    
    # check labels
    if (!is.null(lab_colTree) & !is(lab_colTree, "data.frame")) {
        stop("lab_colTree should be a data frame object")
    }
    
    if (!is.null(lab_rowTree) & !is(lab_rowTree, "data.frame")) {
        stop("lab_rowTree should be a data frame object")
    }
    
    if (!is.null(lab_colTree)) {
        cn <- colnames(lab_colTree)
        isLab <- all(c("oldLab", "newLab") %in% cn)
        if (!isLab) {
            stop("lab_colTree should include two columns: oldLab and newLab.")
        }
        if (anyDuplicated(lab_colTree$oldLab)) {
            stop("duplicated values are not allowed in oldLab of lab_colTree")
        }
        if (anyDuplicated(lab_colTree$newLab)) {
            stop("duplicated values are not allowed in newLab of lab_colTree")
        }
    }
    
    if (!is.null(lab_rowTree)) {
        cn <- colnames(lab_rowTree)
        isLab <- all(c("oldLab", "newLab") %in% cn)
        if (!isLab) {
            stop("lab_rowTree should include two columns: oldLab and newLab.")
        }
        if (anyDuplicated(lab_rowTree$oldLab)) {
            stop("duplicated values are not allowed in oldLab of lab_rowTree")
        }
        if (anyDuplicated(lab_rowTree$newLab)) {
            stop("duplicated values are not allowed in newLab of lab_rowTree")
        }
    }
    
    if (!is.null(new_colTree)) {
        # update the column link data
        cLink <- colLinks(x = tse)
        ii <- match(cLink$nodeLab, lab_colTree$oldLab)
        nLab <- lab_colTree$newLab[ii]
        nNum <- transNode(tree = new_colTree, node = nLab)
        nAlias <- transNode(tree = new_colTree, node = nNum)
        isLf <- isLeaf(tree = new_colTree, node = nNum)
        ncLink <- LinkDataFrame(nodeLab = nLab,
                                nodeLab_alias = nAlias,
                                nodeNum = nNum,
                                isLeaf = isLf)
        tse <- BiocGenerics:::replaceSlots(tse, 
                                           colLinks = ncLink)
    }
    
    if (!is.null(new_rowTree)) {
        # update the rowumn link data
        rLink <- rowLinks(x = tse)
        ii <- match(rLink$nodeLab, lab_rowTree$oldLab)
        nLab <- lab_rowTree$newLab[ii]
        nNum <- transNode(tree = new_rowTree, node = nLab)
        nAlias <- transNode(tree = new_rowTree, node = nNum)
        isLf <- isLeaf(tree = new_rowTree, node = nNum)
        nrLink <- LinkDataFrame(nodeLab = nLab,
                                nodeLab_alias = nAlias,
                                nodeNum = nNum,
                                isLeaf = isLf)
        tse <- BiocGenerics:::replaceSlots(tse, 
                                           rowLinks = nrLink)
    }
    
    return(tse)
}
