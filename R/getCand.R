#' search candidates under different thresholds
#'
#' \code{getCand} search candidates under different thresholds
#'
#' @param tree A phylo object.
#' @param t A sequence of values with the range between 0 and 1. Thresholds used
#'   to search candidates. The default is to use a sequence \code{c(seq(0, 0.04,
#'   by = 0.01), seq(0.05, 1, by = 0.05))}
#' @param score_data A data frame includes at least one column about the nodes,
#'   one column about the p value (\code{p_column}) and one column about the
#'   direction of change (\code{sign_column}).
#' @param node_column The name of the column that gives the node information.
#' @param p_column The name of the column that gives p values of nodes.
#' @param sign_column The name of the column that gives the direction of the
#'   difference.
#' @param threshold Default is 0.05. A value to hinder an internal to be picked
#'   due to randomness.
#' @param pct_na A numeric value. Default is 0.5. Internal nodes that are 
#'   finally selected should have more than half of direct children (> 0.5) 
#'   without missing p values. 
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
                    threshold = 0.05,
                    pct_na = 0.5,
                    message = FALSE) {
    
    if (!is(tree, "phylo")) {
        stop("tree should be a phylo object.")
    }
    
    # t values
    if (is.null(t)) {
        t <- c(0,
               seq(0.01, 0.04, by = 0.01), 
               seq(0.05, 1, by = 0.05))
    }
    
    # a list to store levels under different t 
    level_list <- vector("list", length(t))
    names(level_list) <- c(t)
    
    # columns: p value, sign, node
    p_col <- score_data[[p_column]]
    sign_col <- score_data[[sign_column]]
    node_col <- score_data[[node_column]]
    
    # paths
    path <- matTree(tree = tree)
    
    # for each t
    for (i in seq_along(t)) {
        if (message) {
            message("Searching candidates on t =  ", t[i], " ...")
        }
        # q
        name_q <- paste0("q_", t[i])
        score_data[[name_q]] <- ifelse(p_col > t[i], 0,
                                       1) * sign(sign_col)
        
        # add a column to store the result of search: keep
        if (any(colnames(score_data) == "keep")) {
            stop("The result will be output in the 'keep' column;
             Please use other name for the current 'keep' column.")
        }
        
        keep <- !is.na(p_col)
        node_keep <- score_data[[node_column]][keep]
        node_all <- printNode(tree = tree, type = "all")
        node_in <- node_all[["nodeNum"]][!node_all[["isLeaf"]]]
        node_in <- intersect(node_keep, node_in)
        
        # nodes with p value available
        leaf <- .pseudoLeaf(tree = tree, score_data = score_data, 
                            node_column = node_column,
                            p_column = p_column)
        
        
        # For an internal node, if more than half of its direct child nodes has
        # NA score, it would not be picked.
        chl_I <- findChild(tree = tree, node = node_in)
        sel_0 <- lapply(chl_I, FUN = function(x){
            xx <- match(x, node_col)
            qx <- score_data[[name_q]][xx]
            sum(!is.na(qx))/length(qx) > pct_na
        })
        sel_0 <- unlist(sel_0)
        node_0 <- node_in[sel_0]
        
        
        
        # For an internal nodes, if itself and all its descendant have q score
        # equals 1 or -1, pick the node
        br_I <- findDescendant(tree = tree, node = node_0, only.leaf = FALSE,
                         self.include = TRUE, use.alias = TRUE)
        
        sel_1 <- lapply(br_I, FUN = function(x){
            xx <- match(x, node_col)
            qx <- score_data[[name_q]][xx]
            abs(mean(qx, na.rm = TRUE)) == 1
        })
        sel_1 <- unlist(sel_1)
        
        # filter by threshold
        sel_2 <- p_col[match(node_0, node_col)] <= threshold
        
        node_1 <- node_0[sel_1 & sel_2]
        
        # remove nodes whose ancestor is selected
        ind0 <- apply(path, 2, FUN = function(x) {
            x %in% node_1
        })
        ind0 <- which(ind0, arr.ind = TRUE)
        rs <- split(seq_len(nrow(ind0)), ind0[, "row"])
        rl <- unlist(lapply(rs, length)) > 1
        rs <- rs[rl]
        node_rm <- lapply(rs, FUN = function(x) {
            xx <- ind0[x, , drop = FALSE]
            mx <- xx[, "col"] < max(xx[, "col"])
            fx <- xx[mx, , drop = FALSE]
            path[fx]
        })
        node_rm <- unlist(node_rm)
        node_2 <- setdiff(node_1, node_rm)
        desd_2 <- findDescendant(tree = tree, node = node_2, 
                         only.leaf = FALSE, self.include = TRUE)
        desd_2 <- unlist(desd_2)
        
        level_list[[i]] <- c(setdiff(leaf, desd_2), node_2)
        
    }
    
    out <- list(candidate_list = level_list,
                score_data = score_data)
    return(out)
    
}


