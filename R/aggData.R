#' Perform data aggregations on the tree structure
#'
#' \code{aggData} aggregates data on each node of the tree. Users could decide
#' on which dimension (row or column) and how should the aggregation be
#' performed.
#'
#' @param x A TreeSummarizedExperiment object.
#' @param onRow A logical value, TRUE or FALSE.
#' @param FUN A function for aggregation.
#' @param message A logical value, TRUE or FALSE.
#'
#' @import TreeSummarizedExperiment
#' @importFrom methods is
#' @export
#'
#' @return A TreeSummarizedExperiment object
#' @author Ruizhu HUANG
#'
#' @examples
#' #' # assays data
#' set.seed(1)
#' toyTable <- matrix(rnbinom(20, size = 1, mu = 10), nrow = 5)
#' colnames(toyTable) <- paste(rep(LETTERS[1:2], each = 2),
#'                             rep(1:2, 2), sep = "_")
#' rownames(toyTable) <- paste("entity", seq_len(5), sep = "")
#'
#' toyTable
#'
#' # the column data
#' colInf <- DataFrame(gg = c(1, 2, 3, 3),
#'                     group = rep(LETTERS[1:2], each = 2),
#'                     row.names = colnames(toyTable))
#' colInf
#'
#' # the toy tree
#' library(ape)
#' set.seed(4)
#' treeC <- rtree(4)
#' treeC$node.label <- c("All", "GroupA", "GroupB")
#'
#' library(ggtree)
#' ggtree(treeC, size = 2) +
#'     geom_text2(aes(label = node), color = "darkblue",
#'                hjust = -0.5, vjust = 0.7, size = 6) +
#'     geom_text2(aes(label = label), color = "darkorange",
#'                hjust = -0.1, vjust = -0.7, size = 6)
#'
#'
#' tse <- TreeSummarizedExperiment(assays = list(toyTable),
#'                                 colData = colInf,
#'                                 colTree = treeC,
#'                                 colNodeLab = treeC$tip.label)
#'
#' aggCol <- aggData(x = tse, onRow = FALSE, FUN = sum)
#'

aggData <- function(x, onRow = TRUE, FUN = sum,
                    message = FALSE) {
    # The input data should be a TreeSummarizedExperiment object
    # It should includes only data on the leaf nodes
    if (!is(x, "TreeSummarizedExperiment")) {
        stop("x should be a TreeSummarizedExperiment.")
    }

    # The input data should be on the leaf level
    if (onRow) {
        linkD <- rowLinks(x)
        treeD <- rowTree(x)
    } else {
        linkD <- colLinks(x)
        treeD <- colTree(x)
    }
    isTip <- all(linkD$isLeaf)
    if (!isTip) {
        stop("x should include data only on the leaf level.")
    }

    # All leaves should be available
    leaf <- printNode(tree = treeD, type = "leaf")$nodeNum
    isAll <- setequal(leaf, linkD$nodeNum[linkD$isLeaf])
    if (!isAll) {
        nn <- length(setdiff(leaf, linkD$nodeNum[linkD$isLeaf]))
        stop(nn, " leaves are missing")
    }

    # The aggregation is to each node of the tree
    level <- printNode(tree = treeD, type = "all")$nodeNum

    if (onRow) {
        out <- aggValue(x = x, rowLevel = level)
    } else {
        out <- aggValue(x = x, colLevel = level)
    }

    return(out)

}
