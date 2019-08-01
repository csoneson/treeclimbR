#' search and output results at the optimal level
#' 
#' \code{searchOptimal} search the optimal level on the tree to test the NULL
#' hypotheses and output the result.
#' 
#' @param tree A phylo object.
#' @param score_data A data frame includes at least one column about the nodes,
#'   one column about the p value and one column about the difference direction.
#' @param node_column The name of the column that gives the node information.
#' @param p_column The name of the column that gives p values of nodes.
#' @param sign_column The name of the column that gives the direction of the
#'   difference.
#' @param threshold A sequence of values with the range between 0 and 1. The
#'   different threshold values that the search is performed at.
#' @param score_column The name of the column that will store the generated scores
#' @param method The multiple testing correction method. Please refer to the
#'   argument \code{method} in \code{\link[stats]{p.adjust}}
#' @param output_all A logical value, TRUE or FALSE. Default is TRUE. If TRUE,
#'   intermediate results are saved.
#' @param message A logical value, TRUE or FALSE. Default is FALSE. If TRUE, the
#'   message about running process is printed out.
#'   
#' @importFrom utils flush.console
#' @importFrom methods is
#' @importFrom stats p.adjust
#' @importFrom dplyr mutate "%>%" filter
#' @export
#' @return a list
#' \item{results}{
#'    \itemize{
#'    \item
#'    }}
#' \item{method}{the method used to do multiple testing correction.}
#' \item{threshold}{the threshold values used to search the levels}
#' \item{best_threshold}{which threshold value is used to determine the optimal level}
#' @author Ruizhu Huang
#' @examples 
#' library(TreeSummarizedExperiment)
#' data(tinyTree)
#' set.seed(2)
#' df <- data.frame(node = 1:19, pvalue = runif(19), 
#'                  foldChange = sample(c(1, -1), 19,
#'                   replace = TRUE))
#'                  
#' out <- searchOptimal(tree = tinyTree, score_data = df,
#'                      node_column = "node", p_column = "pvalue",
#'                      sign_column = "foldChange",
#'                      threshold = seq(0, 1, 0.1))
#'                       

searchOptimal <- function(tree, score_data, node_column,
                          p_column, sign_column, threshold,
                          score_column = "S",
                          method = "BH", limit_rej = 0.05,
                          output_all = TRUE, message = FALSE) {
    # ======================= check inputs ===============================
    if (!is(tree, "phylo")) {
        stop("tree should be a phylo object.")
    }
    
    if (anyDuplicated(score_data)) {
        stop("Duplicated rows are detected in the score_data")
    }
    
    if (anyDuplicated(score_data[node_column])) {
        stop("More than one score is detected for a same node")
    }
    
    
    # matrix to store score (S), level, and result
    mm <- matrix(NA, nrow = nrow(score_data), 
                 ncol = length(threshold))
    level_mat <- result_mat <- score_mat <- mm
    colnames(score_mat) <- paste(score_column, 
                                 seq_along(threshold), sep = "_")
    colnames(level_mat) <- paste("level", 
                                 seq_along(threshold), sep = "_")
    colnames(result_mat) <- paste("result", 
                                  seq_along(threshold), sep = "_")
    
    # the number of rejections: on the leaf node
    nleaf <- rep(NA, length(threshold))
    # the number of rejections: on the test level
    nlevel <- rep(NA, length(threshold))
    # the average size of branch with alternative hypotheses
    brSize <- rep(NA, length(threshold))

    
    for (i in seq_along(threshold)) {
        # message
        if (message) {
            message("working on ", i , " out of ",
                    length(threshold), "\r", appendLF = FALSE)
            flush.console()
        }
        
        # the threshold
        t <- threshold[i]
        
        # the score data
        dat_i <- score_data
        
        # S: transform p value to S score
        s <- 1 - dat_i[[p_column]]
        s[!dat_i[p_column] > t] <- 1
        s <- ifelse(dat_i[[sign_column]] > 0 , s, -s)
        dat_i[score_column] <- s
        
        # U: transform S to U
        dat_iu <- scoreTree(tree = tree,
                            score_data = dat_i,
                            node_column = node_column,
                            score_column = score_column,
                            new_score = "U")
        
        # absU: transform U to absU
        score_mat[, i] <- dat_iu[["U"]]
        dat_iu <- dat_iu %>%
            mutate(absU = abs(U))
        
        # the old version
        # score_mat[, i] <- round(dat_iu[["U"]], 2)
        # dat_iu <- dat_iu %>%
        #     mutate(absU = round(abs(U), 2))
        
        # the level selection
        lev <- getLevel(tree = tree, 
                        score_data = dat_iu,
                        score_column = "absU", 
                        node_column = "node",
                        get_max = TRUE, 
                        parent_first = TRUE)
        level_mat[, i] <- lev$keep
        
        # multiple testing correction
        sel <- lev$keep
        pv <- lev[[p_column]][sel]
        adp <- p.adjust(p = pv, method = method)
        rej <- adp <= limit_rej
        nodeSel <- lev[[node_column]][sel]
        nodeRej <- nodeSel[rej]
        result_mat[, i] <- lev[[node_column]] %in% nodeRej
        
        # the rejection on the test level 
        nlevel[i] <- length(nodeRej) 
        
        # the rejection on the leaf level 
        leaf <- findOS(tree = tree, node = nodeRej, 
                       only.leaf = TRUE, self.include = TRUE)
        leaf <- unlist(leaf)
        nleaf[i] <- length(leaf)
        
        # the average branch size
        brSize[i] <- 2*limit_rej*(length(leaf)/length(nodeRej) - 1)
    }
    
    
    
    # The threshold slot
    dt <- cbind.data.frame(T = threshold, nleaf = nleaf,
                           r = brSize, nlevel = nlevel)
    
    # The intermediate results
    result_data <- cbind.data.frame(score_data, score_mat, 
                                    level_mat, result_mat)
    
    # output
    ##  check the validation
    isValid <- brSize >= limit_rej & threshold <= brSize
    if (any(isValid)){
    dts <- dt %>%
        filter(T <= r & r >= limit_rej) %>%
        filter(nleaf == max(nleaf)) %>%
        filter(nlevel == min(nlevel))
    best_t <- dts$T[1]
    best_i <- which(threshold == best_t)
    out <- cbind.data.frame(score_data[[node_column]], 
                            score_data[[p_column]],
                            score_data[[sign_column]],
                            result_mat[, best_i])
    colnames(out) <- c(node_column, p_column, sign_column, "reject")
    
    } else {
       outNode <- setdiff(tree$edge[, 2], tree$edge[, 1])
       outNode <- unique(outNode)
       selDat <- score_data[score_data[[node_column]] %in% outNode, ]
       pv <- selDat[[p_column]]
       adp <- p.adjust(p = pv, method = method)
       rejNode <- selDat[[node_column]][adp <= limit_rej]
       
       out <- cbind.data.frame(score_data[[node_column]], 
                               score_data[[p_column]],
                               score_data[[sign_column]],
                               score_data[[node_column]] %in% rejNode)
       colnames(out) <- c(node_column, p_column, sign_column, "reject")
       
       best_t <- NULL
    }
    
    outList <- list(results = result_data, 
                    method = method, 
                    threshold = dt,
                    best_threshold = best_t,
                    limit = limit_rej,
                    output = out)
    
    if (output_all) {
        return(outList)
    } else {
        return(outList$out) 
    }
    
    
}

# searchOptimal <- function(tree, score_data, node_column,
#                           p_column, sign_column, threshold,
#                           score_column = "S",
#                           method = "BH", limit_rej = 0.05,
#                           message = FALSE) {
#     # ======================= check inputs ===============================
#     if (!is(tree, "phylo")) {
#         stop("tree should be a phylo object.")
#     }
#     
#     if (anyDuplicated(score_data)) {
#         stop("Duplicated rows are detected in the score_data")
#     }
#     
#     if (anyDuplicated(score_data[node_column])) {
#         stop("More than one score is detected for a same node")
#     }
#     
#     
#     score_mat <- matrix(NA, nrow = nrow(score_data), 
#                         ncol = length(threshold))
#     colnames(score_mat) <- paste(score_column, 
#                                  seq_along(threshold), sep = "_")
#     
#     level_mat <- result_mat <- score_mat
#     colnames(level_mat) <- paste("level", 
#                                  seq_along(threshold), sep = "_")
#     colnames(result_mat) <- paste("result", 
#                                   seq_along(threshold), sep = "_")
#     
#     # the number of leaf nodes
#     nleaf <- rep(NA, length(threshold))
#     avSize <- rep(NA, length(threshold))
#     for (i in seq_along(threshold)) {
#         if (message) {
#             message("working on ", i , " out of ",
#                     length(threshold), "\r", appendLF = FALSE)
#             flush.console()
#         }
#         # the threshold
#         t <- threshold[i]
#         
#         # the score data
#         dat_i <- score_data
#         
#         # score value: transformed p value
#         si <- 1 - dat_i[[p_column]]
#         si[!dat_i[p_column] > t] <- 1
#         si <- ifelse(dat_i[[sign_column]] > 0 , si, -si)
#         dat_i[score_column] <- si
#         
#         # calculate node scores
#         # dat_iu <- treeScore(tree = tree,
#         #                     score_data = dat_i,
#         #                     node_column = node_column,
#         #                     score_column = score_column,
#         #                     new_score = "U")
#         dat_iu <- scoreTree(tree = tree,
#                             score_data = dat_i,
#                             node_column = node_column,
#                             score_column = score_column,
#                             new_score = "U")
#         
#         # the old version
#         # score_mat[, i] <- round(dat_iu[["wScore"]], 2)
#         # dat_iu <- dat_iu %>%
#         #     mutate(wScore = round(abs(wScore), 2))
#         
#         # try this new version
#         score_mat[, i] <- dat_iu[["U"]]
#         dat_iu <- dat_iu %>%
#             mutate(absU = abs(U))
#         
#         # search the level using node scores
#         lev <- getLevel(tree = tree, 
#                         score_data = dat_iu,
#                         score_column = "absU", 
#                         node_column = "node",
#                         get_max = TRUE, 
#                         parent_first = TRUE)
#         level_mat[, i] <- lev$keep
#         
#         # run multiple testing correction
#         lev2 <- lev[lev$keep, , drop = FALSE]
#         level <- lev2[[node_column]]
#         pv <- lev2[[p_column]]
#         adp <- p.adjust(p = pv, method = method)
#         node_sel <- unique(level[adp <= limit_rej])
#         result_mat[, i] <- lev[[node_column]] %in% node_sel
#         leaf_sel <- findOS(tree = tree, node = node_sel, 
#                            only.leaf = TRUE, self.include = TRUE)
#         len_leaf <- length(unlist(leaf_sel))
#         # the number of leaf nodes covered by the selected nodes
#         nleaf[i] <- lf
#         rr[i] <- 0.5*limit_rej*(lf/length(node_sel) - 1)
#     }
#     
#     nt <- which(nleaf == max(nleaf))[1]
#     
#     result_data <- cbind.data.frame(score_data, score_mat, 
#                                     level_mat, result_mat)
#     
#     outList <- list(results = result_data, 
#                     method = method, 
#                     threshold = cbind(threshold, rr),
#                     best_threshold = threshold[nt],
#                     limit = limit_rej)
#     
#     return(outList)
#     
# }
# 
