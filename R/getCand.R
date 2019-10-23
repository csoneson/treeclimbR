#' search candidates under different thresholds
#'
#' \code{getCand} search candidates under different thresholds
#'
#' @param tree A phylo object.
#' @param t A sequence of values with the range between 0 and 1.
#'   Thresholds used to search candidates. The default is to use a sequence from
#'   0 to 1 with step 0.05.
#' @param score_data A data frame includes at least one column about the nodes,
#'   one column about the p value (\code{p_column}) and one column about the
#'   direction of change (\code{sign_column}).
#' @param node_column The name of the column that gives the node information.
#' @param p_column The name of the column that gives p values of nodes.
#' @param sign_column The name of the column that gives the direction of the
#'   difference.
#' @param message A logical value, TRUE or FALSE. Default is FALSE. If TRUE, the
#'   message about running process is printed out.
#'
#' @importFrom utils flush.console
#' @importFrom methods is
#' @importFrom stats p.adjust
#' @importFrom dplyr mutate "%>%"
#' @export
#' @return A list includes \code{candidate_list} and \code{score_data}. 
#' \item{condidate_list}{A list of candidates under different thresholds.}
#' \item{score_data}{a data frame that includes columns from the input
#'     \code{score_data} and additional columns to store score U and s under
#'     different thresholds.}
#' @author Ruizhu Huang
#' @examples
#' library(TreeSummarizedExperiment)
#' library(ggtree)
#'
#' data(tinyTree)
#' ggtree(tinyTree, branch.length = "none") +
#'    geom_text2(aes(label = node)) +
#'    geom_hilight(node = 13, fill = "blue", alpha = 0.5) +
#'    geom_hilight(node = 18, fill = "orange", alpha = 0.5)
#' set.seed(2)
#' pv <- runif(19, 0, 1)
#' pv[c(1:5, 13, 14, 18)] <- runif(8, 0, 0.001)
#'
#' fc <- sample(c(-1, 1), 19, replace = TRUE)
#' fc[c(1:3, 13, 14)] <- 1
#' fc[c(4, 5, 18)] <- -1
#' df <- data.frame(node = 1:19,
#'                  pvalue = pv,
#'                  foldChange = fc)
#'
#' ll <- getCand(tree = tinyTree, score_data = df,
#'               node_column = "node", 
#'               p_column = "pvalue", 
#'               sign_column = "foldChange")

getCand <- function(tree, t = NULL,
                    score_data, node_column,
                    p_column, sign_column,
                    message = FALSE) {
    
    if (!is(tree, "phylo")) {
        stop("tree should be a phylo object.")
    }
    
    if (is.null(t)) {
        t <- seq(0, 1, by = 0.05)
    }
    # a list to store levels under different ts
    level_list <- vector("list", length(t) + 1)
    names(level_list) <- c(t, "leaf")
    
    p_col <- score_data[[p_column]]
    sign_col <- score_data[[sign_column]]
    
    for (i in seq_along(t)) {
        if (message) {
            message("Calculating U at t = ", x, " ...")
        }
        
        # S
        name_S <- paste0("S_", t[i])
        score_data[[name_S]] <- ifelse(p_col > t[i], 1-p_col, 1) * sign(sign_col)
        
        # U
        name_U <- paste0("U_", t[i])
        score_data <- treeScore(tree = tree,
                                score_data = score_data,
                                node_column = node_column,
                                score_column = name_S,
                                new_score = name_U)
        
        # U: transform u to U (U = abs(u))
        score_data[[name_U]] <- abs(score_data[[name_U]])
        
        # get levels
        if (message) {
            message("Searching the candidate level at t = ", t, " ...")
        }
        lev <- getLevel(tree = tree,
                        score_data = score_data,
                        score_column = name_U,
                        node_column = node_column,
                        get_max = TRUE,
                        parent_first = TRUE,
                        message = FALSE)
        level_list[[i]] <- lev[[node_column]][lev$keep]
        
    }
    
    # the leaf level
    leaf <- showNode(tree = tree, only.leaf = TRUE)
    level_list$leaf <- leaf
    
    out <- list(candidate_list = level_list,
                score_data = score_data)
    return(out)
    
}
