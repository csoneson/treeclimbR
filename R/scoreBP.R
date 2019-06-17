#' Calculate the score for each node of the tree
#' 
#' \code{scoreBP} calculates score values for nodes of the tree using the
#' Bootstrap method.
#' 
#' @param tse a \code{TreeSummarizedExperiment} object (see
#'   \code{\link[TreeSummarizedExperiment]{TreeSummarizedExperiment}}.
#' @param onRow a logical value, TRUE or FALSE. If true, the score value is
#'   calculated for nodes of the row tree; otherwise, nodes of the column tree.
#' @param groupCol a column name. The column provides information about the
#'   experimental groups that samples belong to.
#' @param normLibSize a logical value, TRUE or FALSE. If true, samples are
#'   scaled to have the same library size that is the maximum library size of
#'   the provided samples.
#' @param n the number of bootstrap.
#' @param assayNum a numeric value to decide which \code{assays} table is used
#'   to calculate the score.
#' 
#' @importFrom SummarizedExperiment assays rowData colData
#' @importFrom TreeSummarizedExperiment rowLinks colLinks rowTree colTree
#' @importFrom stats sd
#' @export
#' @author Ruizhu HUANG
#' @examples 
#' library(TreeSummarizedExperiment)
#' library(ggtree)
#'  
#' # tree
#' data("tinyTree")
#' # We simulate a count table below to have abundance changed only on two
#' # colored branches
#' ggtree(tinyTree, size = 2)+
#'    geom_text2(aes(label = node), color = "darkblue",
#'                hjust = -0.5, vjust = 0.7, size = 6) +
#'    geom_hilight(node = 18, fill = "blue", alpha = 0.3) +
#'    geom_hilight(node = 13, fill = "orange", alpha = 0.3) +
#'    xlim(c(0, 3))
#'     
#'     
#' # count table
#' set.seed(1)
#' p1 <- rep(0.1, 10)
#' p2 <- c(rep(0.04, 3), rep(0.19, 2), rep(0.1, 5))
#' n <- 50   # the number of replicates
#' count <- cbind(rmultinom(n = n, size = 100, prob = p1),
#'               rmultinom(n = n, size = 100, prob = p2))
#' #rowname(count) <- tinyTree$tip.label
#' 
#' rowD <- data.frame(Truth = rep(c(TRUE, FALSE), each = 5),
#'                    stringsAsFactors = FALSE)
#' colD <- data.frame(group = rep(LETTERS[1:2], each = n))
#'  
#'  # build a TreeSummarizedExperiment object to store data on the leaf level
#'  tse1 <- TreeSummarizedExperiment(assays = list(count), 
#'                                   rowData = rowD,
#'                                   colData = colD,
#'                                   rowTree = tinyTree, 
#'                                   rowNodeLab = tinyTree$tip.label)
#'                                   
#' # aggregate data to have data availabel for all nodes
#' # data on an internal node is the sum of counts from its descendants 
#' tse2 <- aggData(x = tse1, onRow = TRUE, FUN = sum, message = TRUE)
#' 
#' #' # output score
#' score <- scoreBP(tse = tse2, onRow = TRUE, 
#'                  groupCol = "group", normLibSize = TRUE,
#'                  n = 100, assayNum = 1)

scoreBP <- function(tse, onRow = TRUE, groupCol, 
                     normLibSize, n = 100, assayNum = 1) {
    
    if (!is(tse, "TreeSummarizedExperiment")) {
        stop("tse requires a TreeSummarizedExperiment object.")
    }
    
    if (onRow) {
        dat <- assays(tse)[[assayNum]]
        linkD <- rowLinks(tse)
        treeD <- rowTree(tse)
        annD <- colData(tse)
    } else {
        dat <- assays(tse)[[assayNum]]
        dat <- t(dat)
        linkD <- colLinks(tse)
        treeD <- colTree(tse)
        annD <- rowData(tse)
    }
    
    if (is.null(treeD)) {
        stop("Tree is not available.")
    }
    # the number of nodes in the tree
    nodeA <- unique(as.vector(treeD$edge))
    nodeA <- sort(nodeA)
    nodeOut <- setdiff(nodeA, linkD$nodeNum)
    if (length(nodeOut)) {
        stop (length(nodeOut), " nodes have no data available")
    }
    
    if (!groupCol %in% colnames(annD)) {
        stop(groupCol, " doesn't exist.")
    }
    
    group <- annD[[groupCol]]
    ug <- unique(group)
    lg <- length(ug)
    if (lg > 2) {
        stop(lg, " groups are found; Only two-group comparison is allowed.")
    }
    rownames(dat) <- linkD$nodeLab_alias
    
    # normalisation on library size
    libSize <- apply(dat, 2, sum)
    if (normLibSize) {
        rr <- max(libSize)/libSize
    } else {
        rr <- rep(1, ncol(dat))
    }
    dat <- dat %*% diag(rr)
    
    # positive group
    scoreP <- vector("list", length = length(nodeA))
    gr <- lapply(ug, FUN = function(x) {
        which(group == x)
    })
    names(gr) <- ug
    selP <- cbind(sample(gr[[1]], n*n, replace = TRUE),
                  sample(gr[[2]], n*n, replace = TRUE))
    colnames(selP) <- ug
    
    # negative group
    scoreN <- scoreP
    selN <- selP
    rg <- sample(seq_along(gr), n*n, replace = TRUE)
    selN[rg == 1, ] <- matrix(sample(gr[[1]], 2*sum(rg == 1), replace = TRUE),
                              ncol = 2)
    selN[rg == 2, ] <- matrix(sample(gr[[2]], 2*sum(rg == 2), replace = TRUE),
                              ncol = 2)
    
    # output
    out <- vector("list", length = length(nodeA))
    for (i in seq_along(nodeA)) {
        # positive
        sp1 <- selP[, 1]
        sp2 <- selP[, 2]
        abP <- matrix(dat[i, sp1] > dat[i, sp2], nrow = n)
        abP <- apply(abP, 2, FUN = function(x){sum(x)/n})
        blP <- matrix(dat[i, sp1] < dat[i, sp2], nrow = n)
        blP <- apply(blP, 2, FUN = function(x){sum(x)/n})
        scp.i <- cbind(a = abP, b = blP, c = 1-abP-blP)
        scoreP[[i]] <- scp.i
        
        # negative
        sn1 <- selN[, 1]
        sn2 <- selN[, 2]
        abN <- matrix(dat[i, sn1] > dat[i, sn2], nrow = n)
        abN <- apply(abN, 2, FUN = function(x){sum(x)/n})
        blN <- matrix(dat[i, sn1] < dat[i, sn2], nrow = n)
        blN <- apply(blN, 2, FUN = function(x){sum(x)/n})
        scn.i <- cbind(a = abN, b = blN, c = 1-abN-blN)
        scoreN[[i]] <- scn.i
        
        # the probability to get a higher score in positive group than in
        # negative
        hP <- apply(scp.i[, "a", drop = FALSE], 1,
                    FUN = function(x){
                        sum(x > scn.i[, "a"], na.rm = TRUE)})
        
        hN <- apply(scp.i[, "b", drop = FALSE], 1,
                    FUN = function(x){
                        sum(x > scn.i[, "b"], na.rm = TRUE)})
        
        # out[[i]] <- c(aucH = sum(hP)/(n*n),
        #               aucL = - sum(hN)/(n*n),
        #               scorePH = .gmean(scp.i[, "a"]),
        #               scorePL = -.gmean(scp.i[, "b"]),
        #               scoreNH = .gmean(scn.i[, "a"]),
        #               scoreNL = -.gmean(scn.i[, "b"]))
        
        shP <- sum(hP)/(n*n)
        shL <- - sum(hN)/(n*n)
        
        if (abs(shP) > abs(shL)) {
            auc <- shP
            aucSE <- sd(hP/n)/sqrt(length(hP))
            score <- .gmean(scp.i[, "a"])
            scoreE <- .gmean(scp.i[, "c"])
        } else {
            auc <- shL
            aucSE <- sd(hN/n)/sqrt(length(hN))
            score <- -.gmean(scp.i[, "b"])
            scoreE <- .gmean(scp.i[, "c"])
        }
        node.i <- linkD$nodeNum[i]
        
        # if all observed values are 0, then set score and auc as NA
        xi <- table(dat[i, ])
        if (all(dat[i, ] == 0) | ncol(dat) < 5 | 
            any(xi/sum(xi) > 0.5)) {
            auc <- NA 
            aucSE <- NA
            score <- NA
            scoreE <- NA
        }
        out[[i]] <- c(nodeNum = node.i,
                      auc = auc, aucSE = aucSE, 
                      score = score, scoreEq = scoreE)
        
        
    }
    
    out <- do.call(rbind, out)
    colnames(out) <- c("nodeNum", "auc", "aucSE", "score", "scoreEq")
    # outList <- list(out = out, scoreP = scoreP, scoreN = scoreP)
    out <- data.frame(out)
    return(out)
    
}

# calculate the geometric mean
.gmean <-  function(x, na.rm=TRUE){
    exp(sum(log(x[x >= 0]), na.rm=na.rm) / length(x))
}
# scoreBP <- function(tse, onRow = TRUE, groupCol, normLibSize = TRUE,
#                   n = 100, assayNum = 1) {
#     score <- .scoreBP(tse = tse, onRow = onRow, groupCol = groupCol,
#                       normLibSize = normLibSize, n = n, 
#                       assayNum = assayNum)
#     # tree
#     if (onRow) {tree <- rowTree(tse)} else {tree <- colTree(tse)}
#     
#     # average score based on the family
#     out <- .scoreAv(tree = tree, scoreData = score)
#     
#     return(out)
# }
# 

# # take average over the family and update the score from the leaf level
# .scoreAv <- function(tree, scoreData) {
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
#         })
#         
#     desdA <- c(as.list(tip), desd)
#     
#     # nodes at different levels
#     mat <- matTree(tree = tree)
#     tempData <- scoreData
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
#         auc.i <- lapply(node.i, FUN = function(x){
#             desd.i <- desdA[[x]]
#             rr.i <- match(desd.i, tempData[, "nodeNum"])
#             auc <- mean(tempData[rr.i, "auc"])
#             return(auc)
#             })
#         aucSE.i <- lapply(node.i, FUN = function(x){
#             desd.i <- desdA[[x]]
#             rr.i <- match(desd.i, tempData[, "nodeNum"])
#             aucSE <- tempData[rr.i, "aucSE"]
#             newSE <- sqrt(sum(aucSE^2))
#             return(newSE)
#         })
#         
#         row.i <- match(node.i, tempData[, "nodeNum"])
#         tempData[row.i, "auc"] <- unlist(auc.i)
#         tempData[row.i, "aucSE"] <- unlist(aucSE.i)
#     }
#     
#     return(tempData)
#     
# }

# Calculate score at each node

