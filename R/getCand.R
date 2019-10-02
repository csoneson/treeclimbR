#' search candidates under different thresholds
#'
#' \code{getCand} search candidates under different thresholds
#'
#' @param tree A phylo object.
#' @param threshold A sequence of values with the range between 0 and 1.
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

getCand <- function(tree, threshold = NULL,
                    score_data, node_column,
                    p_column, sign_column,
                    message = FALSE) {

    if (!is(tree, "phylo")) {
        stop("tree should be a phylo object.")
    }

    if (anyDuplicated(score_data)) {
        stop("Duplicated rows are detected in the score_data")
    }

    if (anyDuplicated(score_data[node_column])) {
        stop("More than one score is detected for a same node")
    }

    if (is.null(threshold)) {
        threshold <- seq(0, 1, by = 0.05)
    }
    # a list to store levels under different thresholds
    level_list <- vector("list", length(threshold) + 1)
    names(level_list) <- c(paste0("level_", threshold), "level_leaf")
    
    dat_U <- score_data
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
        dat_i$S <- s

        # u: transform S to u
        dat_iu <- treeScore(tree = tree,
                            score_data = dat_i,
                            node_column = node_column,
                            score_column = "S",
                            new_score = "u")

        # U: transform u to U (U = abs(u))
        dat_iu <- dat_iu %>%
            mutate(U = abs(u))
        U_i <- paste0("U_", t)
        s_i <- paste0("s_", t)
        dat_U[[U_i]] <- dat_iu$U
        dat_U[[s_i]] <- dat_iu$S
            
       
        
        # get levels
        lev <- getLevel(tree = tree,
                        score_data = dat_iu,
                        score_column = "U",
                        node_column = node_column,
                        get_max = TRUE,
                        parent_first = TRUE)
        level_list[[i]] <- lev[[node_column]][lev$keep]

    }
    
    # the leaf level
    leaf <- showNode(tree = tree, only.leaf = TRUE)
    level_list$level_leaf <- leaf
    
    out <- list(candidate_list = level_list,
                score_data = dat_U)
    return(out)

}
