#' Search for a target level on the tree via a specified score
#'
#' Search for the target level of the tree via a specified score. The score
#' value needs to be provided for each node of the tree.
#'
#' @author Ruizhu Huang
#' @export
#'
#' @param tree A \code{phylo} object.
#' @param score_data A \code{data.frame} providing scores for all nodes in the
#'     tree. The data frame should have at least 2 columns, one with information
#'     about nodes (the node number) and the other with score for each node.
#' @param drop A logical expression indicating elements or rows to keep.
#'     Missing values are taken as FALSE.
#' @param score_column The name of the column of \code{score_data} that
#'     contains the original scores of the nodes.
#' @param node_column The name of the column of \code{score_data} that contains
#'     the numbers of the nodes.
#' @param get_max A logical scalar. If \code{TRUE}, search for nodes
#'     that has higher score value than its descendants; otherwise, search for
#'     nodes that has lower score value than its descendants.
#' @param parent_first A logical scalar. If \code{TRUE}, the parent
#'     node is selected if tied values occur on the parent node and some of the
#'     children nodes.
#' @param message A logical scalar indicating whether progress messages should
#'     be printed.
#'
#' @returns A \code{data.frame} similar to \code{score_data} but with an
#' additional column named \code{keep} indicating which nodes to retain.
#'
#' @importFrom methods is
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom TreeSummarizedExperiment findDescendant printNode
#'
#' @examples
#' library(TreeSummarizedExperiment)
#' library(ggtree)
#' data(tinyTree)
#' ggtree(tinyTree, size = 2, branch.length = "none") +
#'     geom_text2(aes(label = node), color = "darkblue",
#'            hjust = -0.5, vjust = 0.7, size = 4) +
#'     geom_text2(aes(label = label), color = "darkorange",
#'                hjust = -0.1, vjust = -0.7, size = 4) +
#'     geom_hilight(node = 13, fill = "blue", alpha = 0.4) +
#'     geom_hilight(node = 16, fill = "orange", alpha = 0.4)
#'
#' ## Generate score for each node
#' pv <- rep(0.1, 19)
#' pv[c(16, 13, 17)] <- c(0.01, 0.05, 0.005)
#' out <- data.frame(node = 1:19, pvalue = pv)
#'
#' ## search nodes
#' final <- getLevel(tree = tinyTree,
#'                   score_data = out,
#'                   drop =  pvalue > 0.05,
#'                   score_column = "pvalue",
#'                   node_column = "node",
#'                   get_max = FALSE,
#'                   parent_first = TRUE,
#'                   message = FALSE)
#'
#' Nodes to keep
#' final$node[final$keep]
#'
getLevel <- function(tree, score_data, drop, score_column,
                      node_column, get_max, parent_first = TRUE,
                      message = FALSE) {

    ## Check arguments
    ## -------------------------------------------------------------------------
    .assertVector(x = tree, type = "phylo")
    if (!methods::is(score_data, "data.frame")) {
        if (methods::is(score_data, "DataFrame")) {
            score_data <- data.frame(score_data)
        } else {
            stop("score_data should be a data.frame")
        }
    }
    .assertScalar(x = score_column, type = "character",
                  validValues = colnames(score_data))
    .assertScalar(x = node_column, type = "character",
                  validValues = colnames(score_data))
    .assertScalar(x = get_max, type = "logical")
    .assertScalar(x = parent_first, type = "logical")
    .assertScalar(x = message, type = "logical")

    ## Check that there is not already a column named 'keep' (will be added
    ## to hold the result of the search)
    if (any(colnames(score_data) == "keep")) {
        stop("The result will be output in the 'keep' column;
             Please use other name for the current 'keep' column.")
    }

    ## Filter nodes
    ## -------------------------------------------------------------------------
    if (message) {
        message("Preparing data ... ")
    }
    score_data$keep <- TRUE

    ## Drop nodes using the specified standard
    if (message) {
        message("Dropping nodes ... ")
    }
    if (missing(drop)) {
        r <- FALSE
    } else {
        e <- substitute(drop)
        r <- eval(expr = e, envir = score_data, enclos = parent.frame())

        if (!is.logical(r)) {
            stop("'drop' must be or evaluate to logical")
        }
        r <- r & !is.na(r)
    }
    score_data$keep[r] <- FALSE
    if (message) {
        message(sum(!score_data$keep), " nodes are dropped...")
    }

    ## Drop nodes without available score
    score_data$keep[is.na(score_data[, score_column])] <- FALSE
    if (message) {
        message(sum(is.na(score_data[, score_column])),
                " nodes with missing score are dropped...")
    }

    ## Search nodes
    ## -------------------------------------------------------------------------
    # This works only on nodes that are not dropped.
    # It compares the score on an internal node with those on the descendant
    # nodes.
    # 1. get_max = TRUE & parent_first = TRUE
    #   An internal node is selected if its score is not lower than the max
    #   score of its descendants.
    # 2. get_max = TRUE & parent_first = FALSE
    #   An internal node is selected if its score is above the max score
    #   of its descendants.
    # 3. get_max = FALSE & parent_first = TRUE
    #   An internal node is selected if its score is lower or equal to the min
    #   score of its descendants.
    # 4. get_max = FALSE & parent_first = FALSE
    #   An internal node is selected if its score is lower than the min score of
    #   its descendants.

    if (message) {
        message("Searching candidate nodes... ")
    }

    nodeKeep <- score_data[[node_column]][score_data$keep]
    treeData <- TreeSummarizedExperiment::printNode(tree = tree, type = "all")
    nodeIn <- treeData$nodeNum[!treeData$isLeaf]
    nodeI <- intersect(nodeKeep, nodeIn)

    if (message) {
        message("Searching the descendant nodes of the candidate nodes... ")
    }

    descI <- TreeSummarizedExperiment::findDescendant(
        tree = tree, node = nodeI, only.leaf = FALSE,
        self.include = FALSE, use.alias = TRUE)
    chlI <- findChild(tree = tree, node = nodeI)

    if (message) {
        message("Comparing nodes... ")
        if (length(descI)) {
            pb <- utils::txtProgressBar(0, length(descI), style = 3)
        }
    }

    for (i in seq_along(descI)) {
        node.i <- nodeI[i]
        child.i <- chlI[[i]]
        desc.i <- descI[[i]]

        row.p <- match(node.i, score_data[, node_column])
        row.c <- match(child.i, score_data[, node_column])
        row.d <- match(desc.i, score_data[, node_column])

        ## Extract the values on the parent and on the children
        score.p <- score_data[row.p, score_column]
        score.c <- score_data[row.c, score_column]
        score.d <- score_data[row.d, score_column]

        ## If more than half of the direct child nodes has NA score, don't take
        ## the parent node
        na_c <- sum(is.na(score.c))/length(score.c)
        na_d <- sum(is.na(score.d))/length(score.d)
        if (na_c >= 0.5 & na_d != 1) {
            score_data$keep[row.p] <- FALSE
        } else {
            ## If the direct child nodes all have score values, do the following
            ## comparison
            if (get_max) {
                score.p[is.na(score.p)] <- -Inf
                score.d[is.na(score.d)] <- -Inf
                if (parent_first) {
                    isKeep <- all(score.p >= score.d)
                } else {
                    isKeep <- all(score.p > score.d)
                }
            } else {
                score.p[is.na(score.p)] <- Inf
                score.d[is.na(score.d)] <- Inf
                if (parent_first) {
                    isKeep <- all(score.p <= score.d)
                } else {
                    isKeep <- all(score.p < score.d)
                }
            }

            if (isKeep) {
                score_data$keep[row.d] <- FALSE
            } else {
                score_data$keep[row.p] <- FALSE
            }
        }
        if (message) {
            utils::setTxtProgressBar(pb, i)
        }
    }
    return(score_data)
}


