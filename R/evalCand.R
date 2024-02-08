#' Evaluate candidate levels and select the optimal one
#'
#' Evaluate all candidate levels proposed by \code{\link{getCand}} and select
#' the one with best performance. For more details about how the scoring is
#' done, see Huang et al (2021): https://doi.org/10.1186/s13059-021-02368-1.
#'
#' @author Ruizhu Huang
#' @export
#'
#' @param tree A \code{phylo} object.
#' @param type A character scalar indicating whether the evaluation is for a
#'     DA-type workflow (set \code{type="single"}) or a DS-type workflow
#'     (set \code{type="multiple"}).
#' @param levels A list of candidate levels that are returned by
#'     \code{\link{getCand}}. If \code{type = "single"}, elements in the list
#'     are candidate levels, and are named by the value of the tuning parameter.
#'     If \code{type = "multiple"}, a nested list is required and
#'     the list should be named by the feature (e.g., genes). In that case,
#'     each element is a list of candidate levels for that feature.
#' @param score_data A \code{data.frame} (\code{type = "single"}) or a list of
#'     \code{data.frame}s (\code{type = "multiple"}). Each \code{data.frame}
#'     must have at least one column containing the node IDs
#'     (defined by \code{node_column}), one column with p-values
#'     (defined by \code{p_column}), one column with the direction of change
#'     (defined by \code{sign_column}) and one optional column with the feature
#'     (\code{feature_column}, for \code{type="multiple"}).
#' @param node_column The name of the column that contains the node information.
#' @param p_column The name of the column that contains p-values of nodes.
#' @param sign_column The name of the column that contains the direction of the
#'     (estimated) change.
#' @param feature_column The name of the column that contains information about
#'     the feature ID.
#' @param method method The multiple testing correction method. Please refer to
#'     the argument \code{method} in \code{\link[stats]{p.adjust}}. Default is
#'     "BH".
#' @param limit_rej The desired false discovery rate threshold.
#' @param use_pseudo_leaf A logical scalar. If \code{FALSE}, the FDR is
#'     calculated on the leaf level of the tree; If \code{TRUE}, the FDR is
#'     calculated on the pseudo leaf level. The pseudo-leaf level is the level
#'     on which entities have sufficient data to run analysis and the that is
#'     closest to the leaf level.
#' @param message A logical scalar, indicating whether progress messages should
#'     be printed.
#'
#' @returns A list with the following components:
#' \describe{
#'     \item{\code{candidate_best}}{The best candidate level}
#'     \item{\code{output}}{Node-level information for best candidate level}
#'     \item{\code{candidate_list}}{A list of candidates}
#'     \item{\code{level_info}}{Summary information of all candidates}
#'     \item{\code{FDR}}{The specified FDR level}
#'     \item{\code{method}}{The method to perform multiple test correction.}
#'     \item{\code{column_info}}{A list with the specified node, p-value, sign
#'     and feature column names}
#' }
#' More details about the columns in \code{level_info}:
#' \itemize{
#'     \item t The thresholds.
#'     \item r The upper limit of t to control FDR on the leaf level.
#'     \item is_valid Whether the threshold is in the range to control leaf FDR.
#'     \item \code{limit_rej} The specified FDR.
#'     \item \code{level_name} The name of the candidate level.
#'     \item \code{rej_leaf} The number of rejections on the leaf level.
#'     \item \code{rej_pseudo_leaf} The number of rejected pseudo leaf nodes.
#'     \item \code{rej_node} The number of rejections on the tested candidate
#'     level (leaves or internal nodes).
#' }
#'
#' @importFrom utils flush.console
#' @importFrom stats p.adjust
#' @importFrom dplyr select mutate filter
#' @importFrom data.table rbindlist
#' @importFrom TreeSummarizedExperiment findDescendant showNode matTree
#'
#' @examples
#' library(TreeSummarizedExperiment)
#' library(ggtree)
#'
#' data(tinyTree)
#' ggtree(tinyTree, branch.length = "none") +
#'    geom_text2(aes(label = node)) +
#'    geom_hilight(node = 13, fill = "blue", alpha = 0.5) +
#'    geom_hilight(node = 18, fill = "orange", alpha = 0.5)
#' set.seed(1)
#' pv <- runif(19, 0, 1)
#' pv[c(seq_len(5), 13, 14, 18)] <- runif(8, 0, 0.001)
#'
#' fc <- sample(c(-1, 1), 19, replace = TRUE)
#' fc[c(seq_len(3), 13, 14)] <- 1
#' fc[c(4, 5, 18)] <- -1
#' df <- data.frame(node = seq_len(19),
#'                  pvalue = pv,
#'                  foldChange = fc)
#' ll <- getCand(tree = tinyTree, score_data = df,
#'                node_column = "node",
#'                p_column = "pvalue",
#'                sign_column = "foldChange")
#' cc <- evalCand(tree = tinyTree, levels = ll$candidate_list,
#'                score_data = ll$score_data, node_column = "node",
#'                p_column = "pvalue", sign_column = "foldChange",
#'                limit_rej = 0.05)
#'
#' ## Best candidate
#' cc$candidate_best
#'
#' ## Details for best candidate
#' cc$output
#'
evalCand <- function(tree, type = c("single", "multiple"),
                     levels, score_data = NULL,
                     node_column, p_column, sign_column,
                     feature_column = NULL, method = "BH",
                     limit_rej = 0.05, use_pseudo_leaf = FALSE,
                     message = FALSE) {

    ## Check arguments
    ## -------------------------------------------------------------------------
    .assertVector(x = tree, type = "phylo")

    type <- match.arg(type)
    if (type == "single") {
        .assertVector(x = score_data, type = "data.frame")
        .assertVector(x = levels, type = "list")
        score_data <- list(data.frame(score_data))
        levels <- list(levels)
    } else {
        .assertVector(x = score_data, type = "list")
        .assertVector(x = levels, type = "list")
        score_data <- lapply(score_data, data.frame)
    }

    .assertScalar(x = node_column, type = "character")
    .assertScalar(x = p_column, type = "character")
    .assertScalar(x = sign_column, type = "character")
    for (i in seq_along(score_data)) {
        stopifnot(all(c(node_column, p_column, sign_column) %in%
                      colnames(score_data[[i]])))
    }
    .assertScalar(x = feature_column, type = "character", allowNULL = TRUE)
    if (type == "multiple") {
        if (is.null(feature_column)) {
            warning("To distinguish results from different features, ",
                    "the feature_column is required.")
        } else {
            for (i in seq_along(score_data)) {
                stopifnot(feature_column %in% colnames(score_data[[i]]))
            }
        }
    }
    .assertScalar(x = method, type = "character")
    .assertScalar(x = limit_rej, type = "numeric", rngIncl = c(0, 1))
    .assertScalar(x = use_pseudo_leaf, type = "logical")
    .assertScalar(x = message, type = "logical")

    ## Get nodes for each candidate
    ## -------------------------------------------------------------------------
    node_list <- lapply(score_data, FUN = function(x) {
        x[[node_column]]
    })

    ## Find pseudo-leaves (if requested)
    ## -------------------------------------------------------------------------
    ## Some nodes might not be included in the analysis step because they don't
    ## have enough data (invalid p-values). In such case, an internal node
    ## would become a pseudo leaf node if its descendant nodes are filtered
    ## due to lack of sufficient data.
    if (use_pseudo_leaf) {
        ## Find pseudo leaves for each element in score_data
        if (message) {
            message("Finding the pseudo leaf level for all features ...")
        }
        pseudo_leaf <- lapply(seq_along(score_data), FUN = function(x) {
            if (message) {
                message(x, " out of ", length(score_data),
                        " features finished", "\r", appendLF = FALSE)
                utils::flush.console()}
            .pseudoLeaf(tree = tree, score_data = score_data[[x]],
                        node_column = node_column, p_column = p_column)
        })
        names(pseudo_leaf) <- names(score_data)

        ## Count the number of leaves and pseudo leaves for each node
        if (message) {
            message("Calculating the number of pseudo-leaves of each node",
                    "for all features ...")
        }
        info_nleaf <- lapply(seq_along(node_list), FUN = function(x) {
            if (message) {
                message(x, " out of ", length(node_list),
                        " features finished", "\r", appendLF = FALSE)
                utils::flush.console()
            }

            xx <- node_list[[x]]
            ps.x <- pseudo_leaf[[x]]

            desd.x <- TreeSummarizedExperiment::findDescendant(
                tree = tree, node = xx, only.leaf = FALSE, self.include = TRUE)
            leaf.x <- TreeSummarizedExperiment::findDescendant(
                tree = tree, node = xx, only.leaf = TRUE, self.include = TRUE)
            psLeaf.x <- lapply(desd.x, FUN = function(x) {
                intersect(x, ps.x)
            })
            info <- cbind(n_leaf = unlist(lapply(leaf.x, length)),
                          n_pseudo_leaf = unlist(lapply(psLeaf.x, length)))

            return(info)
        })

        names(info_nleaf) <- names(score_data)
    } else {
        ## Count number of leaves for each node
        node_all <- TreeSummarizedExperiment::showNode(
            tree = tree, only.leaf = FALSE)
        desc_all <- TreeSummarizedExperiment::findDescendant(
            tree = tree, node = node_all, only.leaf = TRUE, self.include = TRUE)
        info_nleaf <- data.frame(
            node = node_all,
            n_leaf = unlist(lapply(desc_all, length)))
    }

    ## Evaluate candidates
    ## -------------------------------------------------------------------------
    if (message) {
        message("Evaluating candidates ... ")
    }

    ## Get vector of t-values (should be the same for all features)
    ## -------------------------------------------------------------------------
    tlist <- lapply(levels, names)
    t <- tlist[!duplicated(tlist)]
    if (length(t) > 1) {
        stop("The names of elements in 'levels' are different")
    }
    t <- unlist(t)
    t <- as.numeric(t)

    ## Initialize level_info data.frame
    ## -------------------------------------------------------------------------
    level_info <- data.frame(t = t, upper_t = NA,
                             is_valid = FALSE,
                             method = method,
                             limit_rej = limit_rej,
                             level_name = tlist[[1]],
                             best = FALSE,
                             rej_leaf = NA,
                             rej_node = NA,
                             rej_pseudo_leaf = NA,
                             rej_pseudo_node = NA)

    ## Go through each t-value and collect information
    ## -------------------------------------------------------------------------
    sel <- vector("list", length(t))
    names(sel) <- tlist[[1]]
    for (i in seq_along(t)) {
        if (message) {
            message("Working on ", i , " out of ",
                    length(t), " candidates \r", appendLF = FALSE)
            utils::flush.console()
        }

        ## Get the candidate level at t[i]
        ## ---------------------------------------------------------------------
        name_i <- as.character(level_info$level_name[i])
        level_i <- lapply(levels, FUN = function(x) {
            x[[i]]
        })

        ## Get the nodes in the candidate
        ## ---------------------------------------------------------------------
        sel_i <- mapply(function(x, y) {
            ii <- match(x, y[[node_column]])
        }, level_i, score_data, SIMPLIFY = FALSE)
        len_i <- lapply(sel_i, length)

        ## Adjust p-values
        ## ---------------------------------------------------------------------
        p_i <- mapply(FUN = function(x, y) {
            x[y, p_column]
        }, score_data, sel_i, SIMPLIFY = FALSE)
        adp_i <- stats::p.adjust(p = unlist(p_i), method = method)
        rej_i <- adp_i <= limit_rej

        ## The largest p value that is rejected
        ## ---------------------------------------------------------------------
        maxp_i <- max(c(-1, unlist(p_i)[rej_i]))

        ## Number of rejected branches
        ## ---------------------------------------------------------------------
        path <- TreeSummarizedExperiment::matTree(tree = tree)
        n_C <- mapply(FUN = function(x, y) {
            ## Nodes rejected in each feature
            xx <- x[y, c(node_column, sign_column, p_column)]
            xs <- xx[xx[[p_column]] <= maxp_i, ]

            ## Split nodes by sign
            sn <- split(xs[[node_column]], sign(xs[[sign_column]]))

            is_L <- lapply(sn, FUN = function(x) {
                isLeaf(tree = tree, node = x)})
            rej_L <- mapply(FUN = function(x, y) {
                unique(x[y])}, sn, is_L)
            rej_I <- mapply(FUN = function(x, y) {
                unique(x[!y]) }, sn, is_L)
            rej_L2 <- lapply(rej_L, FUN = function(x) {
                unique(path[path[, "L1"] %in% x, "L2"])})
            length(unlist(rej_I)) + length(unlist(rej_L2))
        }, score_data, sel_i, SIMPLIFY = FALSE)
        n_C <- sum(unlist(n_C))

        ## Number of rejected leaves
        ## ---------------------------------------------------------------------
        if (use_pseudo_leaf) {
            rej_m1 <- mapply(FUN = function(x, y) {
                x[y, "n_pseudo_leaf"]
            }, info_nleaf, sel_i, SIMPLIFY = FALSE)
            n_m <- sum(unlist(rej_m1)[rej_i %in% TRUE])
            av_size <- n_m/max(n_C, 1)
        } else {
            node_i <- mapply(FUN = function(x, y) {
                x[y, node_column]
            }, score_data, sel_i, SIMPLIFY = FALSE)
            node_r <- unlist(node_i)[rej_i %in% TRUE]
            ind_r <- match(node_r, info_nleaf[["node"]])
            n_m <- sum(info_nleaf[ind_r, "n_leaf"])
            av_size <- n_m/max(n_C, 1)
        }

        ## Estimate the maximal t that still controls the FDR and add
        ## columns to level_info
        ## ---------------------------------------------------------------------
        up_i <- min(2 * limit_rej * (max(av_size, 1) - 1), 1)
        up_i <- round(up_i, 10)

        level_info$upper_t[i] <- up_i
        level_info$rej_leaf[i] <- n_m
        level_info$rej_node[i] <- sum(rej_i)

        if (use_pseudo_leaf) {
            level_info$rej_pseudo_leaf[i] <- n_m
            level_info$rej_pseudo_node[i] <- n_C
        }
        sel[[i]] <- sel_i

        level_info$is_valid[i] <- up_i > t[i] | t[i] == 0
    }

    ## Compare candidates and find the best one
    ## -------------------------------------------------------------------------
    ## candidates: levels that fulfill the requirement to control FDR on the
    ## (pseudo) leaf level when multiple hypothesis correction is performed
    isB <- level_info |>
        dplyr::filter(is_valid) |>
        dplyr::filter(rej_leaf == max(rej_leaf)) |>
        dplyr::filter(rej_node == min(rej_node)) |>
        dplyr::select(level_name) |>
        unlist() |>
        as.character()
    level_info <- level_info |>
        dplyr::mutate(best = level_name %in% isB)
    level_b <- lapply(levels, FUN = function(x) {
        x[[isB[1]]]
    })

    ## Output the result on the best level
    ## -------------------------------------------------------------------------
    if (message) {
        message("Multiple-hypothesis correction on the best candidate ...")
    }
    sel_b <- sel[[isB[1]]]

    outB <- lapply(seq_along(score_data), FUN = function(i) {
        si <- sel_b[[i]]
        score_data[[i]][si, , drop = FALSE]
    })

    outB <- data.table::rbindlist(outB)
    pv <- outB[[p_column]]
    apv <- stats::p.adjust(pv, method = method)
    outB$adj.p <- apv
    outB$signal.node <- apv <= limit_rej

    ## Assemble final output
    ## -------------------------------------------------------------------------
    if (message) {
        message("output the results ...")
    }
    out <- list(candidate_best = level_b,
                output = outB,
                candidate_list = levels,
                level_info = level_info,
                FDR = limit_rej,
                method = method,
                column_info = list("node_column" = node_column,
                                   "p_column" = p_column,
                                   "sign_column" = sign_column,
                                   "feature_column" = feature_column))
    return(out)

}

#' @author Ruizhu  Huang
#' @keywords internal
#'
#' @returns A vector of node numbers representing the pseudo leaf level
#'
#' @importFrom TreeSummarizedExperiment matTree
#'
.pseudoLeaf <- function(tree, score_data, node_column, p_column) {
    ## Create matrix with paths from leaves to root
    ## -------------------------------------------------------------------------
    mat <- TreeSummarizedExperiment::matTree(tree = tree)

    ## Check which nodes in each path have valid p-values
    ## -------------------------------------------------------------------------
    nd <- score_data[[node_column]][!is.na(score_data[[p_column]])]
    exist_mat <- apply(mat, 2, FUN = function(x) {
        x %in% nd
    })

    ## Find leaves with valid values
    ## -------------------------------------------------------------------------
    ww <- which(exist_mat, arr.ind = TRUE)
    ww <- ww[order(ww[, 1]), , drop = FALSE]
    loc_leaf <- ww[!duplicated(ww[, 1]), ]
    leaf_0 <- unique(mat[loc_leaf])

    ## Find lowest nodes with valid values
    ## -------------------------------------------------------------------------
    ind_0 <- lapply(leaf_0, FUN = function(x) {
        xx <- which(mat == x, arr.ind = TRUE)
        ux <- xx[!duplicated(xx), , drop = FALSE]
        y0 <- nrow(ux) == 1
        if (nrow(ux) > 1) {
            ux[, "col"] <- ux[, "col"] - 1
            y1 <- all(!exist_mat[ux])
            y0 <- y0 | y1
        }
        return(y0)
    })
    leaf_1 <- leaf_0[unlist(ind_0)]

    return(leaf_1)
}

