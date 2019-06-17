
#' Generate average score in a family
#' 
#' \code{avgScore} takes the tree structure into account when calculating the
#' score for an internal node. If an internal node A has two children B and C,
#' \code{(A->B, A->C)}, the new score at node A would be calculated as the
#' mean of the scores in the whole family (A, B and C). 
#' 
#' @param tree A phylo object
#' @param scoreData A data frame that includes at least two columns. One column
#'   stores the number of the node, and the other stores the original score of
#'   the corresponding node.
#' @param nodeCol The name of the column that stores the number of the node.
#' @param scoreCol The name of the column that stores the original score of the
#'   node.
#' 
#' @importFrom TreeSummarizedExperiment printNode matTree
#' @export
#' 
#' @author Ruizhu HUANG
#' 
#' @examples 
#' 
#' library(TreeSummarizedExperiment)
#' library(ggtree)
#' library(dplyr)
#' 
#' # tree
#' data("tinyTree")
#' ggtree(tinyTree) + geom_text2(aes(label = node))
#' 
#' # score
#' exScore <- data.frame(nodeNum = 1:19, score = (1:19)/10)
#'
#' # update the score based on the family provided by the tree
#' newScore <- avgScore(tree = tinyTree, scoreData = exScore,
#'                       nodeCol = "nodeNum", scoreCol = "score")
#'
#' # visualize the result
#' df <- newScore %>% 
#'       rename(node = nodeNum) %>%    
#'       mutate(score = round(score, 3)) %>%
#'       mutate(avgScore = round(avgScore, 3))
#' 
#' # the original scores are in black texts and the new ones in blue
#' ggtree(tinyTree) %<+% df + 
#' geom_text2(aes(label = score), hjust = -0.05) +
#' geom_text2(aes(label = avgScore, hjust = -1.2), color = "blue") +
#' xlim(0, 5)  



avgScore <- function(tree, scoreData, nodeCol, 
                       scoreCol) {
    
    df <- printNode(tree, type = "all")
    nodeA <- sort(df$nodeNum)
    tip <- sort(df$nodeNum[df$isLeaf])
    nodeI <- sort(setdiff(nodeA, tip))
    
    if (!all.equal(nodeA, c(tip, nodeI))) {
        warnings("The nodeNum should start from leaves with number 1.")
    }
    
    
    desd <- lapply(nodeI, FUN = function(x) {
        xx <- findChild(tree = tree, node = x)
        out <- c(x, xx)
        return(out)
    })
    
    desdA <- c(as.list(tip), desd)
    score <- lapply(desdA, FUN = function(x) {
        ix <- match(x, scoreData[, nodeCol])
        mean(scoreData[ix, scoreCol])
    })
    outData <- scoreData
    outData$avgScore <- unlist(score)[match(c(tip, nodeI), 
                                            scoreData[, nodeCol])]
    return(outData)
    
}
