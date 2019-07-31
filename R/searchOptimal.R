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
#' @param message A logical value, TRUE or FALSE. Default is FALSE. If TRUE, the
#'   message about running process is printed out.
#'   
#' @importFrom utils flush.console
#' @importFrom methods is
#' @importFrom stats p.adjust
#' @importFrom dplyr mutate "%>%"
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
#' set.seed(1)
#' df <- data.frame(node = 1:19, pvalue = runif(19), 
#'                  foldChange = sample(c(1, -1), 19,
#'                   replace = TRUE))
#'                  
#' out <- searchOptimal(tree = tinyTree, score_data = df,
#'                      node_column = "node", p_column = "pvalue",
#'                      sign_column = "foldChange",
#'                       threshold = seq(0, 1, 0.1))
#'                       
searchOptimal <- function(tree, score_data, node_column,
                          p_column, sign_column, threshold,
                          score_column = "treeP",
                          method = "BH", limit_rej = 0.05,
                          message = FALSE) {
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
    
    
    score_mat <- matrix(NA, nrow = nrow(score_data), 
                        ncol = length(threshold))
    colnames(score_mat) <- paste(score_column, 
                                 seq_along(threshold), sep = "_")
    
    level_mat <- result_mat <- score_mat
    colnames(level_mat) <- paste("level", 
                                 seq_along(threshold), sep = "_")
    colnames(result_mat) <- paste("result", 
                                  seq_along(threshold), sep = "_")
    
    # the number of leaf nodes
    nleaf <- rep(NA, length(threshold))
    
    for (i in seq_along(threshold)) {
        if (message) {
            message("working on ", i , " out of ",
                    length(threshold), "\r", appendLF = FALSE)
            flush.console()
        }
        # the threshold
        t <- threshold[i]
        
        # the score data
        dat_i <- score_data
        
        # score value: transformed p value
        si <- 1 - dat_i[[p_column]]
        si[!dat_i[p_column] > t] <- 1
        si <- ifelse(dat_i[[sign_column]] > 0 , si, -si)
        dat_i[score_column] <- si
        
        # calculate node scores
        # dat_iu <- treeScore(tree = tree,
        #                     score_data = dat_i,
        #                     node_column = node_column,
        #                     score_column = score_column,
        #                     new_score = "wScore")
        dat_iu <- scoreTree(tree = tree,
                            score_data = dat_i,
                            node_column = node_column,
                            score_column = score_column,
                            new_score = "wScore")
        
        # the old version
        # score_mat[, i] <- round(dat_iu[["wScore"]], 2)
        # dat_iu <- dat_iu %>%
        #     mutate(wScore = round(abs(wScore), 2))
        
        # try this new version
        score_mat[, i] <- dat_iu[["wScore"]]
        dat_iu <- dat_iu %>%
            mutate(wScore = abs(wScore))
        
        # search the level using node scores
        out <- getLevel(tree = tree, 
                        score_data = dat_iu,
                        score_column = "wScore", 
                        node_column = "node",
                        get_max = TRUE, 
                        parent_first = TRUE)
        level_mat[, i] <- out$keep
        
        # run multiple testing correction
        out_keep <- out[out$keep, , drop = FALSE]
        level <- out_keep[[node_column]]
        pv <- out_keep[[p_column]]
        adp <- p.adjust(p = pv, method = method)
        node_sel <- level[adp <= limit_rej]
        result_mat[, i] <- out[[node_column]] %in% node_sel
        leaf_sel <- findOS(tree = tree, node = node_sel, 
                           only.leaf = TRUE, self.include = TRUE)
        
        # the number of leaf nodes covered by the selected nodes
        nleaf[i] <- length(unlist(leaf_sel))
    }
    
    nt <- which(nleaf == max(nleaf))[1]
    
    result_data <- cbind.data.frame(score_data, score_mat, 
                                    level_mat, result_mat)
    outList <- list(results = result_data, method = method, 
                    threshold = threshold,
                    best_threshold = threshold[nt],
                    limit = limit_rej)
    
    return(outList)
    
}

