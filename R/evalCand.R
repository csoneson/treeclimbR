#' evaluate candidate levels and select the best one
#'
#' \code{evalCand} evaluate all candidate levels and select the one with best
#' performance
#'
#' @param tree A phylo object.
#' @param type "single" or "multiple".
#' @param levels A list of candidate levels that are selected by
#'   \code{\link{getCand}}. If \code{type = "DA"}, elements in the list are
#'   candidate levels, and are named by values of tuning parameter that are
#'   used. If \code{type = "DS"}, a nested list is required and the list should
#'   be named by the feature (e.g., genes or antibodies). Each element is a list
#'   of candidate levels for a feature (e.g. gene or antibody) that are selected
#'   by \code{\link{getCand}}.
#' @param score_data A data frame (\code{type = "DA"}) or a of data frame
#'   (\code{type = "DS"}). Each data frame includes at least one column about
#'   the nodes (\code{node_column}), one column about the p value
#'   (\code{p_column}), one column about the direction of change
#'   (\code{sign_column}) and one optional column about the feature
#'   (\code{feature_column}, this is to distinct the results from different
#'   features for \code{type = "DS"} in the final output.)
#' @param node_column The name of the column that gives the node information.
#' @param p_column The name of the column that gives p values of nodes.
#' @param sign_column The name of the column that gives the direction of the
#'  (estimated) change.
#' @param feature_column The name of the column that gives information about the
#'   feature.
#' @param method method The multiple testing correction method. Please refer to
#'   the argument \code{method} in \code{\link[stats]{p.adjust}}. Default is
#'   "BH".
#' @param limit_rej The FDR level. Default is 0.05.
#' @param control_fdr_on "leaf" or "pseudo-leaf". If "leaf", the FDR is
#'   calculated on the leaf level of the tree; If "pseudo-leaf", the FDR is
#'   calculated on the pseudo leaf level. The pseudo-leaf level is the level on
#'   which entities have sufficient data to run analysis and the level that is
#'   closest to the leaf level.
#' @param message A logical value, TRUE or FALSE. Default is FALSE. If TRUE, the
#'   message about running process is printed out.
#'
#' @importFrom utils flush.console
#' @importFrom methods is
#' @importFrom stats p.adjust
#' @importFrom dplyr select
#' @importFrom data.table rbindlist
#' @importFrom TreeSummarizedExperiment findOS
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
#'  \item t the thresholds
#'  \item r the upper limit of T to control FDR on the leaf level
#'  \item is_valid whether the threshold is in the range to control leaf FDR
#'  \item \code{limit_rej} the specified FDR
#'  \item \code{level_name} the name of the candidate level
#'  \item \code{rej_leaf} the number of rejection on the leaf level
#'  \item \code{rej_pseudo_leaf} the number of rejected pseudo leaf nodes.
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
#'               #t = seq(0, 1, by = 0.05),
#'                node_column = "node",
#'                p_column = "pvalue",
#'                sign_column = "foldChange")
#' cc <- evalCand(tree = tinyTree, levels = ll$candidate_list,
#'                score_data = df, node_column = "node",
#'                p_column = "pvalue", sign_column = "foldChange",
#'                limit_rej = 0.05 )
#' cc$output
evalCand <- function(tree, 
                     type = c("single", "multiple"),
                     levels = cand_list,
                     score_data = NULL, 
                     node_column, p_column,
                     sign_column = sign_column,
                     feature_column = NULL,
                     method = "BH", 
                     limit_rej = 0.05,
                     # control_fdr_on = c("leaf", "pseudo-leaf"),
                     message = FALSE) {
    
    if (!is(tree, "phylo")) {
        stop("tree should be a phylo object.")
    }
    
    type <- match.arg(type)
    if (type == "single") {
        score_data <- list(score_data)
        levels <- list(levels)
    }
    
    
    if (type == "multiple" & is.null(feature_column)) {
        warning("To distinct results from different features, 
                feature_column is required")
    }
    
   
    # ------------------------- the pseudo leaf level -------------------------
    # some nodes might not be included in the analysis step because they have no
    # enough data. In such case, an internal node would become a pseudo leaf 
    # node if its descendant nodes are filtered due to lack of sufficient data.
    
    if (message) {
        message("collecting the pseudo leaf level for all features ...")
    }
    
    pseudo_leaf <- lapply(seq_along(score_data), FUN = function(x) {
        if (message) {
            message(x, " out of ", length(score_data),
                    " features finished", "\r", appendLF = FALSE)
            flush.console()}
        
        .pseudoLeaf(tree = tree, score_data = score_data[[x]], 
                    node_column = node_column, p_column = p_column)
    })
    names(pseudo_leaf) <- names(score_data)
    
    node_list <- lapply(score_data, FUN = function(x) {
        x[[node_column]]
    })
    
    if (message) {
        message("Calculating the number of pseudo-leaves of each node 
                for all features ...")
    }
    info_nleaf <- lapply(seq_along(node_list), FUN = function(x) {
        if (message) {
            message(x, " out of ", length(node_list),
                    " features finished", "\r", appendLF = FALSE)
            flush.console()}

        xx <- node_list[[x]]
        ps.x <- pseudo_leaf[[x]]

        desd.x <- findOS(tree = tree, node = xx,
                         only.leaf = FALSE, self.include = TRUE)
        leaf.x <- findOS(tree = tree, node = xx,
                         only.leaf = TRUE, self.include = TRUE)
        psLeaf.x <- lapply(desd.x, FUN = function(x) {
            intersect(x, ps.x)})
        info <- cbind(n_leaf = unlist(lapply(leaf.x, length)),
                      n_pseudo_leaf = unlist(lapply(psLeaf.x, length)))

        return(info)
    })
    names(info_nleaf) <- names(score_data)
    
    # add two columns in score_data
    # ---------------- info about the candidate level --------------------------
    # candidates in the candidate level
    # a data frame: t, br_size, candidate, method, limit_rej, level_name, 
    # rej_leaf, rej_node, rej_pseudo_leaf
    
    if (message) {
        message("Evaluating candidates ... ")
    }
    tlist <- lapply(levels, names)
    t <- tlist[!duplicated(tlist)]
    if (length(t) > 1) {
        stop("the names of elements in 'levels' are different")
    }
    t <- unlist(t)
    t[t == "leaf"] <- NA
    t <- as.numeric(t)
    
    level_info <- data.frame(t = t, lower_t = NA, upper_t = NA, 
                             is_valid = FALSE,
                             method = method, limit_rej = limit_rej,
                             level_name = tlist[[1]],
                             rej_leaf = NA, rej_node = NA,
                             rej_pseudo_leaf = NA,
                             rej_pseudo_node = NA)
    
    sel <- vector("list", length(t))
    names(sel) <- tlist[[1]]
    for (i in seq_along(t)) {
        message(i)
        # message
        if (message) {
            message("working on ", i , " out of ",
                    length(t), " candidates \r", appendLF = FALSE)
            flush.console()
        }
        
        # get the candidate level at t[i]
        name_i <- as.character(level_info$level_name[i])
        level_i <- lapply(levels, FUN = function(x) {
            x[[i]]
        })
        sel_i <- mapply(function(x, y) {
            ii <- match(x, y[[node_column]])
        }, level_i, score_data, SIMPLIFY = FALSE)
        len_i <- lapply(sel_i, length)
        # ind_i <- mapply(FUN = function(x, y) { rep(x, y)}, 
        #                 seq_along(sel_i), len_i, SIMPLIFY = FALSE)
        
        
        # adjust p-values
        p_i <- mapply(FUN = function(x, y) {
            x[y, p_column]
        }, score_data, sel_i, SIMPLIFY = FALSE)
        adp_i <- p.adjust(p = unlist(p_i), method = method)
        rej_i <- adp_i <= limit_rej
        
        # the low bound of t
        low_i <- min(unlist(p_i)[rej_i])
        
        # calculate r
        # rejected nodes
        nd_i <- mapply(FUN = function(x, y) {
            x[y, node_column]
        }, score_data, sel_i, SIMPLIFY = FALSE)
        sign_i <- mapply(FUN = function(x, y) {
            sign(x[y, sign_column])
        }, score_data, sel_i, SIMPLIFY = FALSE)
        rej_sign <- unlist(sign_i)[rej_i %in% TRUE]
        rej_N <- unlist(nd_i)[rej_i %in% TRUE]
        rej_N <- split(rej_N, rej_sign)
        
        # rejected branches
        path <- matTree(tree = tree)
        is_L <- lapply(rej_N, FUN = function(x) {
            isLeaf(tree = tree, node = x)})
        rej_L <- mapply(FUN = function(x, y) {
            unique(x[y])}, rej_N, is_L)
        rej_I <- mapply(FUN = function(x, y) {
            unique(x[!y]) }, rej_N, is_L)
        rej_L2 <- lapply(rej_L, FUN = function(x) {
            unique(path[path[, "L1"] %in% x, "L2"])}) 
        n_C <- length(rej_I) + length(unlist(rej_L2))
        
        
        rej_m1 <- mapply(FUN = function(x, y) {
            x[y, "n_pseudo_leaf"]
        }, info_nleaf, sel_i, SIMPLIFY = FALSE)
        n_m1 <- sum(unlist(rej_m1)[rej_i %in% TRUE])
        up_i <- 2 * limit_rej * (n_m1/max(n_C, 1) - 1)
        
        
        rej_m2 <- mapply(FUN = function(x, y) {
            x[y, "n_leaf"]
        }, info_nleaf, sel_i, SIMPLIFY = FALSE)
        n_m2 <- sum(unlist(rej_m2)[rej_i %in% TRUE])
        
        level_info$lower_t[i] <- low_i
        level_info$upper_t[i] <- up_i
        level_info$rej_leaf[i] <- n_m2
        level_info$rej_pseudo_leaf[i] <- n_m1
        level_info$rej_node[i] <- length(unlist(rej_N))
        level_info$rej_pseudo_node[i] <- n_C
        sel[[i]] <- sel_i
        
        # the leaf is always valid; 
        # the aggregated level should be in a specific range to be valid
        if (is.na(t[i]) & name_i == "leaf") {
            level_info$is_valid[i] <- TRUE
        } else {
            level_info$is_valid[i] <- up_i >= t[i] & t[i] >= low_i
        }
    }
    
    # candidates: levels that fullfil the requirement to control FDR on the 
    # (pseudo) leaf level when multiple hypothesis correction is performed on it
    isB <- level_info %>%
        filter(is_valid) %>%
        filter(rej_pseudo_leaf == max(rej_pseudo_leaf)) %>%
        filter(rej_node == min(rej_node)) %>%
        select(level_name) %>% 
        unlist() %>%
        as.character()
    
    level_b <- lapply(levels, FUN = function(x) {x[[isB[1]]]})
    
    # output the result on the best level
    if (message) {
        message("mulitple-hypothesis correction on the best candidate ...")
    }
    sel_b <- sel[[isB[1]]]
    
    outB <- lapply(seq_along(score_data), FUN = function(i) {
        si <- sel_b[[i]]
        score_data[[i]][si, , drop = FALSE]
    })
    
    outB <- rbindlist(outB)
    pv <- outB[[p_column]]
    apv <- p.adjust(pv, method = method)
    outB$adj.p <- apv
    outB$signal.node <- apv <= limit_rej
    
    if (message) {
        message("output the results ...")
    }
    out <- list(candidate_best = level_b, output = outB,
                candidate_list = levels,  
                level_info = level_info, 
                FDR = limit_rej, method = method)
    return(out)
    
}

#' @importFrom TreeSummarizedExperiment matTree
#' @keywords internal
.pseudoLeaf <- function(tree, score_data, node_column, p_column) {
    mat <- matTree(tree = tree)
    nd <- score_data[[node_column]][!is.na(score_data[[p_column]])]
    exist_mat <- apply(mat, 2, FUN = function(x) {x %in% nd})
    
    ww <- which(exist_mat, arr.ind = TRUE)
    ww <- ww[order(ww[, 1]), , drop = FALSE]
    loc_leaf <- ww[!duplicated(ww[, 1]), ]
    leaf_0 <- unique(mat[loc_leaf])
    ind_0 <- lapply(leaf_0, FUN = function(x) {
        xx <- which(mat == x, arr.ind = TRUE)
        ux <- xx[!duplicated(xx), , drop = FALSE]
        y0 <- nrow(ux) == 1
        if (nrow(ux) > 1) {
            ux[, "col"] <- ux[, "col"] - 1
            y1 <- all(is.na(exist_mat[ux]))
            y0 <- y0 | y1
        }
        return(y0)
    })
    leaf_1 <- leaf_0[unlist(ind_0)]
    
    
    return(leaf_1)
}

