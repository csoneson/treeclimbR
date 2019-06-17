#' Generate weighted score accounting for the family effect
#' 
#' \code{wScoreTree} takes the tree structure into account when calculating the
#' score for an internal node. If an internal node A has two children B and C,
#' \code{(A->B, A->C)}, the new score at node A would be calculated as the
#' weighted mean of the scores in the whole family (A, B and C). The weights are
#' based on the number of descendant leaves. For example, if the node B has 2
#' descendant leaves, and C has 3 descendant leaves, then A would have 5. The
#' calculation would be \eqn{(Score\_A * 5 + Score\_B*2 +Score\_C*3)/10}. The
#' generation starts from the leaves and the new generated scores are used to
#' update those in higher level of the tree until the root is reached.
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
#' newScore <- wScoreTree(tree = tinyTree, scoreData = exScore,
#'                       nodeCol = "nodeNum", scoreCol = "score")
#'
#' # visualize the result
#' df <- newScore %>% 
#'       rename(node = nodeNum) %>%    
#'       mutate(score = round(score, 3)) %>%
#'       mutate(wScore = round(wScore, 3))
#' 
#' # the original scores are in black texts and the new ones in blue
#' ggtree(tinyTree) %<+% df + 
#' geom_text2(aes(label = score), hjust = -0.05) +
#' geom_text2(aes(label = wScore, hjust = -1.2), color = "blue") +
#' xlim(0, 5)  



wScoreTree <- function(tree, scoreData, nodeCol, 
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
    
    # the number of descendant leaves for each node
    leafA <- findOS(tree = tree, node = nodeA, 
                    only.leaf = TRUE, self.include = TRUE)
    leafDF <- data.frame(node = nodeA, num = unlist(lapply(leafA, length)))
    leafA <- lapply(desdA, FUN = function(x) {
        leafDF$num[x]
    })
    # leafA <- lapply(desdA[1:4], FUN = function(x) {
    #     xx <- findOS(tree = tree, node = x, only.leaf = TRUE, 
    #                  self.include = TRUE)
    #     xl <- lapply(xx, length)
    #     ul <- unlist(xl)
    #     return(ul)
    # })
    
    # nodes at different levels
    mat <- matTree(tree = tree)
    tempData <- scoreData
    tempData$wScore <- tempData[, scoreCol]
    
    # rule out leaves
    tempMat <- mat[, -1]
    for (i in seq_len(ncol(tempMat))) {
        mat.i <- tempMat[, i]
        mat.o <- tempMat[, -seq_len(i)]
        node.i <- setdiff(as.vector(mat.i), as.vector(mat.o))
        node.i <- node.i[!is.na(node.i)]
        
        # take the weighted average score of the whole family
        # (a-b; a-c) => scoreNA = wmean(scoreA, scoreB, scoreC) 
        score.i <- lapply(node.i, FUN = function(x){
            desd.i <- desdA[[x]]
            nLeaf.i <- leafA[[x]]
            
            rr.i <- match(desd.i, tempData[, nodeCol])
            sc.i <- tempData[rr.i, "wScore"]
            isNa <- is.na(sc.i)
            sc.i <- sc.i[!isNa]
            lf.i <- nLeaf.i[!isNa]
            ot.i <- sum(sc.i*lf.i)/sum(lf.i)
            return(ot.i)
        })
        
        row.i <- match(node.i, tempData[, nodeCol])
        tempData[row.i, "wScore"] <- unlist(score.i)
        
    }
    
    return(tempData)
    
}
