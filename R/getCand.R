#' Generate candidates for different thresholds
#'
#' Generate candidates for different thresholds (t). A candidate consists of
#' a disjoint collection of leaves and internal branches, that collectively
#' cover all leaves in the tree, and represents a specific aggregation pattern
#' along the tree.
#'
#' @author Ruizhu Huang
#' @export
#'
#' @param tree A \code{phylo} object.
#' @param t A vector of threshold values used to search for candidates,
#'     in the range [0, 1]. The default (\code{NULL}) uses a sequence
#'     \code{c(seq(0, 0.04, by = 0.01), seq(0.05, 1, by = 0.05))}
#' @param score_data A \code{data.frame} including at least one column with
#'     node IDs (specified with the \code{node_column} argument),
#'     one column with p-values (specified with the \code{p_column} argument)
#'     and one column with directions of change (specified with the
#'     \code{sign_column} argument).
#' @param node_column The name of the column of \code{score_data} that
#'     contains the node information.
#' @param p_column The name of the column of \code{score_data} that
#'     contains p-values for nodes.
#' @param sign_column The name of the column of \code{score_data} that
#'     contains the direction of change (e.g., the log-fold change). Only
#'     the sign of this column will be used.
#' @param threshold Numeric scalar; any internal node where the value of
#'     the p-value column is above this value will not be returned. The default
#'     is 0.05. The aim of this threshold is to avoid arbitrarily picking up
#'     internal nodes without true signal.
#' @param pct_na Numeric scalar. In order for an internal node to be eligible
#'     for selection, more than \code{pct_na} of its direct child nodes must
#'     have a valid (i.e., non-missing) value in the \code{p_column} column.
#'     Hence, increasing this number implies a more strict selection (in terms
#'     of presence of explicit values).
#' @param message A logical scalar, indicating whether progress messages
#'     should be printed to the console.
#'
#' @returns A list with two elements: \code{candidate_list} and
#' \code{score_data}. \code{condidate_list} is a list of candidates obtained
#' for the different thresholds. \code{score_data} is a \code{data.frame}
#' that includes columns from the input \code{score_data} and additional
#' columns with q-scores for different thresholds.
#'
#' @importFrom TreeSummarizedExperiment matTree printNode findDescendant
#' @importFrom dplyr bind_cols
#'
#' @examples
#' suppressPackageStartupMessages({
#'     library(TreeSummarizedExperiment)
#'     library(ggtree)
#' })
#'
#' data(tinyTree)
#' ggtree(tinyTree, branch.length = "none") +
#'    geom_text2(aes(label = node)) +
#'    geom_hilight(node = 13, fill = "blue", alpha = 0.3) +
#'    geom_hilight(node = 18, fill = "orange", alpha = 0.3)
#'
#' ## Simulate p-values and directions of change for nodes
#' ## (Nodes 1, 2, 3, 4, 5, 13, 14, 18 have a true signal)
#' set.seed(1)
#' pv <- runif(19, 0, 1)
#' pv[c(seq_len(5), 13, 14, 18)] <- runif(8, 0, 0.001)
#'
#' fc <- sample(c(-1, 1), 19, replace = TRUE)
#' fc[c(seq_len(3), 13, 14)] <- 1
#' fc[c(4, 5, 18)] <- -1
#' df <- data.frame(node = seq_len(19),
#'                  pvalue = pv,
#'                  logFoldChange = fc)
#'
#' ll <- getCand(tree = tinyTree, score_data = df,
#'               t = c(0.01, 0.05, 0.1, 0.25, 0.75),
#'               node_column = "node", p_column = "pvalue",
#'               sign_column = "logFoldChange")
#'
#' ## Candidates
#' ll$candidate_list
#'
#' ## Score table
#' ll$score_data
#'
getCand <- function(tree, t = NULL, score_data, node_column, p_column,
                    sign_column, threshold = 0.05, pct_na = 0.5,
                    message = FALSE) {

    ## Check arguments and initialize variables
    ## -------------------------------------------------------------------------
    .assertVector(x = tree, type = "phylo")
    .assertVector(x = t, type = "numeric", rngIncl = c(0, 1),
                  allowNULL = TRUE)
    .assertVector(x = score_data, type = "data.frame")
    .assertScalar(x = node_column, type = "character")
    .assertScalar(x = p_column, type = "character")
    .assertScalar(x = sign_column, type = "character")
    stopifnot(all(c(node_column, p_column, sign_column) %in%
                      colnames(score_data)))
    .assertScalar(x = threshold, type = "numeric", rngIncl = c(0, 1))
    .assertScalar(x = pct_na, type = "numeric", rngIncl = c(0, 1))
    .assertScalar(x = message, type = "logical")

    ## Set t values if missing
    ## -------------------------------------------------------------------------
    if (is.null(t)) {
        t <- c(0, seq(0.01, 0.04, by = 0.01), seq(0.05, 1, by = 0.05))
    }

    ## Initiate a list to store levels under different t
    ## -------------------------------------------------------------------------
    level_list <- vector("list", length(t))
    names(level_list) <- c(t)

    ## Get columns: p-value, sign, node
    ## -------------------------------------------------------------------------
    p_col <- score_data[[p_column]]
    sign_col <- score_data[[sign_column]]
    node_col <- score_data[[node_column]]

    ## Transform the tree into a matrix. Matrix entries are node numbers, and
    ## each row represents a path connecting a leaf node and the root
    ## -------------------------------------------------------------------------
    path <- TreeSummarizedExperiment::matTree(tree = tree)

    ## Generate candidates
    ## -------------------------------------------------------------------------
    ## All nodes
    node_all <- TreeSummarizedExperiment::printNode(tree = tree, type = "all")

    ## Initialize q columns (one for each t)
    qcols <- do.call(dplyr::bind_cols,
                     lapply(structure(seq_along(t),
                                      names = paste0("q_", t)), function(i) {
        rep(NA_real_, nrow(score_data))
    }))
    score_data <- dplyr::bind_cols(score_data, qcols)

    for (i in seq_along(t)) {
        if (message) {
            message("Searching candidates on t = ", t[i], " ...")
        }

        ## Add a q column (does the node have pvalue <= t[i], add sign)
        ## ---------------------------------------------------------------------
        name_q <- paste0("q_", t[i])
        score_data[[name_q]] <- ifelse(p_col > t[i], 0, 1) * sign(sign_col)

        ## Get vector of node numbers with valid p-value
        ## ---------------------------------------------------------------------
        keep <- !is.na(p_col)
        node_keep <- score_data[[node_column]][keep]

        ## Internal nodes with valid p-values
        ## ---------------------------------------------------------------------
        node_in <- node_all[["nodeNum"]][!node_all[["isLeaf"]]]
        node_in <- intersect(node_keep, node_in)

        ## Get leaves (with valid p-values)
        ## ---------------------------------------------------------------------
        leaf <- .pseudoLeaf(tree = tree, score_data = score_data,
                            node_column = node_column, p_column = p_column)

        ## For an internal node, if more than pct_na of its direct child nodes
        ## have a valid score, it will be retained. If less than pct_na of its
        ## direct child nodes have a valid score, it will be excluded.
        ## ---------------------------------------------------------------------
        chl_I <- findChild(tree = tree, node = node_in)
        sel_0 <- lapply(chl_I, FUN = function(x) {
            xx <- match(x, node_col)
            qx <- score_data[[name_q]][xx]
            sum(!is.na(qx))/length(qx) > pct_na
        })
        sel_0 <- unlist(sel_0)
        node_0 <- node_in[sel_0]

        ## For an internal node, if itself and all its descendants have q score
        ## equal to 1 or -1, pick the node
        ## ---------------------------------------------------------------------
        br_I <- TreeSummarizedExperiment::findDescendant(
            tree = tree, node = node_0, only.leaf = FALSE,
            self.include = TRUE, use.alias = TRUE)

        ## Generate a logical vector (with names being internal node aliases)
        ## indicating whether to pick the node
        sel_1 <- vapply(br_I, FUN = function(x) {
            xx <- match(x, node_col)
            qx <- score_data[[name_q]][xx]
            abs(mean(qx, na.rm = TRUE)) == 1
        }, FUN.VALUE = FALSE)

        ## Filter by threshold
        ## ---------------------------------------------------------------------
        sel_2 <- p_col[match(node_0, node_col)] <= threshold
        node_1 <- node_0[sel_1 & sel_2]

        ## Remove nodes whose ancestor is selected
        ## ---------------------------------------------------------------------
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

        ## Leaves that are not descendants of any selected internal node will
        ## be reported individually
        ## ---------------------------------------------------------------------
        desd_2 <- TreeSummarizedExperiment::findDescendant(
            tree = tree, node = node_2, only.leaf = FALSE, self.include = TRUE)
        desd_2 <- unlist(desd_2)

        level_list[[i]] <- c(setdiff(leaf, desd_2), node_2)
    }

    ## Return results
    ## -------------------------------------------------------------------------
    out <- list(candidate_list = level_list,
                score_data = score_data)

    return(out)
}


