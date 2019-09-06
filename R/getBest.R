#' run levels to find candidates that are expected to control FDR on leaf level
#'
#' \code{getBest} run levels to find candidates that are expected to control
#' FDR on leaf level
#'
#' @param tree A phylo object.
#' @param levels A list of levels that is the output of
#'   \code{\link{searchLevel}}.
#' @param score_data A data frame includes at least one column about the nodes,
#'   one column about the p value (\code{p_column}) and one column about the
#'   direction of change (\code{sign_column}).
#' @param node_column The name of the column that gives the node information.
#' @param p_column The name of the column that gives p values of nodes.
#' @param sign_column The name of the column that gives the direction of the
#'   difference.
#' @param method method The multiple testing correction method. Please refer to the
#'   argument \code{method} in \code{\link[stats]{p.adjust}}. Default is "BH".
#' @param limit_rej The FDR level. Default is 0.05.
#' @param message A logical value, TRUE or FALSE. Default is FALSE. If TRUE, the
#'   message about running process is printed out.
#'
#' @importFrom utils flush.console
#' @importFrom methods is
#' @importFrom stats p.adjust
#' @importFrom dplyr select
#' @export
#' @return a list.
#'   \describe{
#'   \item{\code{candidate_best}}{the best candidate level}
#'   \item{\code{output}}{the result of best candidate level}
#'   \item{\code{candidate_list}}{a list of candidates}
#'   \item{\code{level_info}}{the information of all candidates}
#'   \item{FDR}{the specified FDR level}
#'   \item{method}{the method to perform multiple test correction.}
#'   }
#'  More details about columns in \code{level_info}.
#'  \itemize{
#'  \item T the thresholds
#'  \item r the upper limit of T to control FDR on the leaf level
#'  \item is_valid whether the threshold is in the range to control leaf FDR
#'  \item \code{limit_rej} the specified FDR
#'  \item \code{level_name} the name of the candidate level
#'  \item \code{rej_leaf} the number of rejection on the leaf level
#'  \item \code{rej_node} the number of rejection on the tested candidate level
#'  }
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
#' ll <- getCand(tree = tinyTree, score_data = df,
#'                   node_column = "node", 
#'                   p_column = "pvalue", 
#'                   sign_column = "foldChange")
#' cc <- getBest(tree = tinyTree, 
#'               levels = ll$candidate_list, 
#'               score_data = df,
#'               node_column = "node", 
#'               p_column = "pvalue", 
#'               sign_column = "foldChange",
#'               limit_rej = 0.01)
#'

getBest <- function(tree, levels,
                    score_data, node_column,
                    p_column, sign_column,
                    method = "BH", limit_rej = 0.05,
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

    # ---------------- info about the agg level ---------------------------
    # candidates in the aggregated level
    # a data frame: T, br_size, candidate, method, limit_rej, level_name
    t <- gsub(pattern = "level_", replacement = "", x = names(levels))
    t[t == "leaf"] <- NA
    t <- as.numeric(t)
    
    level_info <- data.frame(T = t, r = NA, is_valid = FALSE,
                           method = method, limit_rej = limit_rej,
                           level_name = names(levels),
                           rej_leaf = NA, rej_node = NA)
    
    for (i in seq_along(t)) {
        # message
        if (message) {
            message("working on ", i , " out of ",
                    length(t), "\r", appendLF = FALSE)
            flush.console()
        }

        # estimate the signal branch size 
        name_i <- as.character(level_info$level_name[i])
        level_i <- levels[[name_i]]
        sel_i <- match(level_i, score_data[[node_column]])
        
        pv <- score_data[[p_column]][sel_i]
        adp <- p.adjust(p = pv, method = method)
        
        rej <- adp <= limit_rej
        node_s <- score_data[[node_column]][sel_i]
        node_r <- node_s[rej]
        
        leaf_r <- findOS(tree = tree, node = node_r, 
                         only.leaf = TRUE, self.include = TRUE)
        leaf_r <- unlist(leaf_r)
        leaf_r <- unique(leaf_r)
        level_info$rej_node[i] <- length(node_r)
        level_info$rej_leaf[i] <- length(leaf_r)
        
        r_i <- 2*limit_rej*(length(leaf_r)/max(c(length(node_r)-1, 1)))
        level_info$r[i] <- r_i
        
        # the leaf is always valid; 
        # the aggregated level should be in a specific range to be valid
        if (is.na(t[i]) & name_i == "level_leaf") {
            level_info$is_valid[i] <- TRUE
        } else {
            level_info$is_valid[i] <- r_i >= limit_rej & t[i] <= r_i
        }
        
    }
    
    # candidates: levels that fullfil the requirement to control FDR on the leaf
    # level when multiple hypothesis correction is performed on it
    level_c <- levels
    isB <- level_info %>%
        filter(is_valid) %>%
        filter(rej_leaf == max(rej_leaf)) %>%
        filter(rej_node == min(rej_node)) %>%
        select(level_name) %>% 
        unlist() %>%
        as.character()
    level_b <- levels[[isB[1]]]
    
    # output the results on the best level
    sel_b <- score_data[[node_column]] %in% level_b
    dat_b <- score_data[sel_b, ]
    pv <- dat_b[[p_column]]
    adp <- p.adjust(p = pv, method = method)
    adp_v <- rep(NA, nrow(score_data))
    adp_v[sel_b] <- adp
    rej_b <- score_data[[node_column]][sel_b][adp <= limit_rej]
    n_leaf <- findOS(tree = tree, node = score_data[[node_column]], 
                     only.leaf = TRUE, self.include = TRUE)
    n_leaf <- unlist(lapply(n_leaf, length))
    
    outB <- cbind.data.frame(score_data[[node_column]], 
                             score_data[[p_column]],
                             score_data[[sign_column]],
                             adp_v,
                             score_data[[node_column]] %in% rej_b,
                             n_leaf)
    colnames(outB) <- c(node_column, p_column, sign_column, 
                        "adj.p", "signal.node", "n_leaf")
    out <- list(candidate_best = level_b, output = outB,
                candidate_list = level_c,  
                level_info = level_info, 
                FDR = limit_rej, method = method)
    return(out)

}
