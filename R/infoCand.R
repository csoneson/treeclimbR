#' Get information of candidates
#' 
#' \code{infoCand} extract information about candidates 
#' 
#' @param object the output of the function \link{evalCand}.
#' @export
#' @return a data frame
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
#'               #t = seq(0, 1, by = 0.05),
#'                node_column = "node",
#'                p_column = "pvalue",
#'                sign_column = "foldChange")
#' cc <- evalCand(tree = tinyTree, levels = ll$candidate_list,
#'                score_data = df, node_column = "node",
#'                p_column = "pvalue", sign_column = "foldChange",
#'                limit_rej = 0.05)
#'                
#' out <- infoCand(object = cc)

infoCand <- function(object){
    info <- object$level_info
    
    if (is.null(info)) {
        stop("obj should be the output from 'evalCand'")
    } 
    
    isNA <- all(is.na(info$rej_pseudo_leaf))
    if (isNA) {
        info <- subset(info, select = -c(rej_pseudo_leaf, rej_pseudo_node))
    }
  
    return(info)
}