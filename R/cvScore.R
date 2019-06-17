#' calculte cross validation error
#' 
#' \code{cvEr} calculates the cross validation error.
#' 
#' @param tse A TreeSummarizedExperiment object.
#' @param kFold A numeric value. K-fold cross validation.
#' @param node A numeric or character vector. 
#' @param onRow A logical value, TRUE or FALSE.
#' @param groupCol The name of the column that stores the group information of
#'   samples.
#' @param normLibSize A logical value, TRUE or FALSE.
#' @param n A numeric value. The number of pairs.
#' @param assayNum A numeric value.
#' 
#' @importFrom TreeSummarizedExperiment rowLinks rowTree colLinks colTree
#' @importFrom SummarizedExperiment colData rowData assays
#' @export
#' @return a numeric value
#' @author Ruizhu HUANG
#' @examples 
#' library(TreeSummarizedExperiment)
#' library(ape)
#' set.seed(1)
#' exTree <- rtree(1000)
#' 
#' p <- runif(1000)
#' p <- p/sum(p)
#' y <- matrix(rmultinom(50, size =5000, prob = p),nrow=1000)
#' 
#' colnames(y) <- paste("S", 1:50, sep = "")
#' rownames(y) <- exTree$tip.label
#'
#'
#' toy_lse <- TreeSummarizedExperiment(assays = list(y), 
#'                                     rowTree = exTree)
#' 
#'
#' toy_lse <- parEstimate(obj = toy_lse)
#' dat1 <- simData(obj = toy_lse, ratio = 4, scenario = "S1",
#'                 minTip.A = 20, maxTip.A = 40,
#'                 minTip.B = 20)
#' # view simulation scenario
#' viewSim(dat1)
#' (br <- metadata(dat1)$branch)
#' 
#' # Aggregation
#' dat2 <- aggData(x = dat1, onRow = TRUE, FUN = sum)
#' # output score
#' score <- scoreBP(tse = dat2, onRow = TRUE, 
#'                  groupCol = "group", normLibSize = TRUE,
#'                  n = 100, assayNum = 1)
#' score <- round(score, digits = 4)
#' # search level
#' final <- searchLevel(tree = rowTree(dat2),
#'                      scoreData = data.frame(score),
#'                      scoreCol = "auc",
#'                      searchMax = TRUE,
#'                      parentPrior = TRUE, 
#'                      message = FALSE)
#'                      
#' # the nodes selected
#' nodeS <- final$nodeNum[final$keep]                  
#'                      
#' # cross validation error
#' cv <- lapply(seq_along(nodeS), FUN = function(x) {
#'                      cat(x, "out of ", length(nodeS), " has finished \n")
#'                      cvEr(tse = dat2, kFold = 10, node = nodeS[x])} )
#' out <- data.frame(nodeNum = nodeS, err = unlist(cv))
#' out <- out[order(out$err, decreasing = FALSE), ]


cvEr <- function(tse, kFold = 1, node, 
               onRow = TRUE, groupCol = "group", 
               normLibSize = TRUE,
               n = 100, assayNum = 1) {
    # check the input object is TreeSummarizedExperiment
    if (!is(tse, "TreeSummarizedExperiment")) {
        stop("tse requires a TreeSummarizedExperiment object.")
    }

    # decide the dimension: row or column
    if (onRow) {
        dat <- assays(tse)[[assayNum]]
        linkD <- rowLinks(tse)
        treeD <- rowTree(tse)
        annD <- colData(tse)
        sp <- ncol(tse)
    } else {
        dat <- assays(tse)[[assayNum]]
        dat <- t(dat)
        linkD <- colLinks(tse)
        treeD <- colTree(tse)
        annD <- rowData(tse)
        sp <- nrow(tse)
    }

    # groups
    group <- annD[[groupCol]]
    nr <- which(linkD$nodeNum == node)
    
    # check the tree is available
    if (is.null(treeD)) {
        stop("Tree is not available.")
    }

    # split data into k folds
    ind <- .createFolds(m = sp, k = kFold)
    
    # test and training
    er <- rep(NA, kFold)
    for (i in seq_len(kFold)) {
        ts <- unlist(ind[i])
        tr <- unlist(ind[-i])
        
        # The group that the (test & training) samples belong to
        gts <- group[ts]
        gtr <- group[tr]
        
        # samples in two different groups in the training data
        gr <- split(seq_along(gtr), gtr)
        
        # test & train data
        test <- dat[, ts]
        train <- dat[, tr]
        
        predG <- rep(NA, ncol(test))
        # run samples in test data
        for (j in seq_len(ncol(test))) {
            sel1 <- sample(gr[[1]], n, replace = TRUE)
            sel2 <- sample(gr[[2]], n, replace = TRUE)
            
            score <- cbind(abs(train[nr, sel1] - test[nr, j]),
                           abs(train[nr, sel2] - test[nr, j]))
            nn <- matrix(NA, nrow = n, ncol =2)
            for (k in seq_len(n)){
                nn[k, 1] <- sum(score[k, 1] > score[, 2], na.rm = TRUE)
                nn[k, 2] <- sum(score[k, 2] > score[, 1], na.rm = TRUE)
            }
            sn <- c(sum(nn[, 1])/(n*n), sum(nn[, 2])/(n*n))
            
            if (sn[1] < sn[2]) {
                predG[[j]] <- names(gr)[1]
            }
            
            if (sn[1] > sn[2]) {
                predG[[j]] <- names(gr)[2]
            }
        }
        
        er[i] <- sum(predG != gts)/length(gts)
        
        
    }
    
    return(mean(er))
    

}


# df <- data.frame(value = c(train[nr, sel1], train[nr, sel2]),
#                  group = rep(names(gr), c(length(sel1), length(sel2))))
# df <- data.frame(value = c(dat[nr, 1:50], dat[nr, 1:50]),
#                  group = rep(names(gr), each = 50))
# library(ggplot2)
# ggplot(data = df) +
#     geom_boxplot(aes(group, value, color = group))


.createFolds <- function(m, k) {
    # shuffle samples
    mm <- sample(seq_len(m), m, replace = FALSE)
    
    # k folds
    ks <- rep(seq_len(k), len = m)
    
    out <- split(x = mm, f = ks)
    
    return(out)
    
}
