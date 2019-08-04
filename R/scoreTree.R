#' Generate score accounting for the family effect
#' 
#' \code{scoreTree} takes the tree structure into account when calculating the
#' score for an internal node. If an internal node A has two children B and C,
#' \code{(A->B, A->C)}, the new score at node A would be calculated as the mean
#' (sum or other user-specified function) of the whole family (A, B and C). The
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
#' @param new_score The name of the column to store the new generated score.
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
#' newScore <- scoreTree(tree = tinyTree, 
#'                       score_data = exScore,
#'                       node_column = "nodeNum", 
#'                       score_column = "score",
#'                       new_score = "treeScore")
#'
#' # visualize the result
#' df <- newScore %>% 
#'       rename(node = nodeNum) %>%    
#'       mutate(score = round(score, 3)) %>%
#'       mutate(treeScore = round(treeScore, 3))
#' 
#' # the original scores are in black texts and the new ones in blue
#' ggtree(tinyTree) %<+% df + 
#' geom_text2(aes(label = score)) +
#' geom_text2(aes(label = treeScore, hjust = -1), 
#' color = "blue")   



scoreTree <- function(tree, score_data, node_column,
                      score_column, new_score) {
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
    
    # nodes at different levels
    mat <- matTree(tree = tree)
    tempData <- score_data
    tempData[[new_score]] <- tempData[, score_column]
    
    # rule out leaves
    tempMat <- mat[, -1]
    for (i in seq_len(ncol(tempMat))) {
        mat.i <- tempMat[, i]
        mat.o <- tempMat[, -seq_len(i)]
        node.i <- setdiff(as.vector(mat.i), as.vector(mat.o))
        node.i <- node.i[!is.na(node.i)]
        
        # take the average score of the whole family
        # (a-b; a-c) => scoreNA = mean(scoreA, scoreB, scoreC) 
        score.i <- lapply(node.i, FUN = function(x){
            desd.i <- desdA[[x]]
            rr.i <- match(desd.i, tempData[, node_column])
            temp.i <- tempData[rr.i, new_score]
            temp.i <- cbind(temp.i[!is.na(temp.i)])
            auc <- apply(temp.i, 2, FUN = function(x) {
                mean(x, na.rm = TRUE)})
            return(auc)
        })
        # aucSE.i <- lapply(node.i, FUN = function(x){
        #     desd.i <- desdA[[x]]
        #     rr.i <- match(desd.i, tempData[, "nodeNum"])
        #     aucSE <- tempData[rr.i, "aucSE"]
        #     newSE <- sqrt(sum(aucSE^2))
        #     return(newSE)
        # })
        # 
        row.i <- match(node.i, tempData[, node_column])
        tempData[row.i, new_score] <- unlist(score.i)
        #tempData[row.i, "aucSE"] <- unlist(aucSE.i)
    }
    
    return(tempData)
}
# function(tree, score_data, node_column, 
#                       score_column, FUN) {
#     
#     df <- printNode(tree, type = "all")
#     nodeA <- sort(df$nodeNum)
#     tip <- sort(df$nodeNum[df$isLeaf])
#     nodeI <- sort(setdiff(nodeA, tip))
#     
#     if (!all.equal(nodeA, c(tip, nodeI))) {
#         warnings("The nodeNum should start from leaves with number 1.")
#     }
#     
#     
#     desd <- lapply(nodeI, FUN = function(x) {
#         xx <- findChild(tree = tree, node = x)
#         out <- c(x, xx)
#         return(out)
#     })
#     
#     desdA <- c(as.list(tip), desd)
#     
#     # nodes at different levels
#     mat <- matTree(tree = tree)
#     tempData <- score_data
#     tempData$treeScore <- tempData[, score_column]
#     
#     # rule out leaves
#     tempMat <- mat[, -1]
#     for (i in seq_len(ncol(tempMat))) {
#         mat.i <- tempMat[, i]
#         mat.o <- tempMat[, -seq_len(i)]
#         node.i <- setdiff(as.vector(mat.i), as.vector(mat.o))
#         node.i <- node.i[!is.na(node.i)]
#         
#         # take the average score of the whole family
#         # (a-b; a-c) => scoreNA = mean(scoreA, scoreB, scoreC) 
#         score.i <- lapply(node.i, FUN = function(x){
#             desd.i <- desdA[[x]]
#             rr.i <- match(desd.i, tempData[, node_column])
#             temp.i <- tempData[rr.i, "treeScore"]
#             temp.i <- cbind(temp.i[!is.na(temp.i)])
#             auc <- apply(temp.i, 2, FUN = FUN)
#             return(auc)
#         })
#         # aucSE.i <- lapply(node.i, FUN = function(x){
#         #     desd.i <- desdA[[x]]
#         #     rr.i <- match(desd.i, tempData[, "nodeNum"])
#         #     aucSE <- tempData[rr.i, "aucSE"]
#         #     newSE <- sqrt(sum(aucSE^2))
#         #     return(newSE)
#         # })
#         # 
#         row.i <- match(node.i, tempData[, node_column])
#         tempData[row.i, "treeScore"] <- unlist(score.i)
#         #tempData[row.i, "aucSE"] <- unlist(aucSE.i)
#     }
#     
#     return(tempData)
#     
# }
