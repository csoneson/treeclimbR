#' Generate weighted score accounting for the family effect
#' 
#' \code{treeScore} takes the tree structure into account when calculating the
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
#' @param score_data A data frame that includes at least two columns. One column
#'   stores the number of the node, and the other stores the original score of
#'   the corresponding node.
#' @param node_column The name of the column that stores the number of the node.
#' @param score_column The name of the column that stores the original score of the
#'   node.
#' @param new_score The name of the column that stores the generated score.
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
#' newScore <- treeScore(tree = tinyTree, score_data = exScore,
#'                       node_column = "nodeNum", 
#'                       score_column = "score",
#'                       new_score = "wScore")
#'
#' # visualize the result
#' df <- newScore %>% 
#'       rename(node = nodeNum) %>%    
#'       mutate(score = round(score, 3)) %>%
#'       mutate(wScore = round(wScore, 3))
#' 
#' # the original scores are in black texts and the new ones in blue
#' ggtree(tinyTree) %<+% 
#'    df + 
#'    geom_text2(aes(label = score), hjust = -0.05) +
#'    geom_text2(aes(label = wScore, hjust = -1.2), 
#'    color = "blue") +
#'    xlim(0, 5)  



treeScore <- function(tree, score_data, node_column, 
                      score_column, new_score) {
    
    df <- printNode(tree, type = "all")
    nodeA <- sort(df$nodeNum)
    tip <- sort(df$nodeNum[df$isLeaf])
    nodeI <- sort(setdiff(nodeA, tip))
    
    if (!all.equal(nodeA, c(tip, nodeI))) {
        warnings("The nodeNum should start from leaves with number 1.")
    }
    
    # a list: each element (a node and its children)
    chl <- findChild(tree = tree, node = nodeI)
    desd <- lapply(seq_along(nodeI), FUN = function(x) {
        c(nodeI[x], chl[[x]])
    })
    desdA <- c(as.list(tip), desd)
    
    # the number of descendant leaves for each node
    leafA <- findOS(tree = tree, node = nodeA, 
                    only.leaf = TRUE, self.include = TRUE)
    leafDF <- data.frame(node = nodeA, 
                         num = unlist(lapply(leafA, length)))
    leafA <- lapply(desdA, FUN = function(x) {
        leafDF$num[x]
    })
    
    
    # nodes at different levels
    mat <- matTree(tree = tree)
    tempData <- score_data
    tempData[new_score] <- tempData[, score_column]
    
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
            
            rr.i <- match(desd.i, tempData[, node_column])
            sc.i <- tempData[rr.i, new_score]
            isNa <- is.na(sc.i)
            sc.i <- sc.i[!isNa]
            lf.i <- nLeaf.i[!isNa]
            ot.i <- sum(sc.i*lf.i)/sum(lf.i)
            return(ot.i)
        })
        
        # remove nodes without scores
        na.i <- unlist(lapply(score.i, is.na))
        score.ii <- score.i[!na.i]
        node.ii <- node.i[!na.i]
        
        row.i <- match(node.ii, tempData[, node_column])
        # if (any(is.na(row.i))) {
        #     warnings("score_data doesn't provide all nodes required.")
        # }
        tempData[row.i, new_score] <- unlist(score.ii)
        
    }
    
    return(tempData)
    
}
