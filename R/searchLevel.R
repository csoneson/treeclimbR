#' search for a target level on the tree via a specified score
#'
#' \code{searchLevel} searches for the target level of the tree via a specified
#' score. The score value needs to be provided for each node of the tree.
#'
#' @param scoreData a data frame to provide scores for all nodes of the tree.
#'   The data frame should have at least five columns, four of them naming
#'   \code{nodeLab}, \code{nodeLab_alias}, \code{nodeNum} and \code{isLeaf}.
#'   These four columns provide the information about nodes. The score column
#'   could be named freely.
#' @param tree a phylo object to provide the hierarchical structure of nodes.
#' @param filter a logical expression indicating elements or rows to keep:
#'   missing values are taken as false.
#' @param scoreCol the name of the score column. The search is based on the
#'   score value.
#' @param nodeCol the name of the column that stores the node number.
#' @param searchMax a logical value, TRUE or FALSE. If TRUE, search for nodes
#'   that has higher score value than its descendants; otherwise, search for
#'   nodes that has lower score value than its descendants.
#' @param parentPrior a logical value, TRUE or FALSE. If TRUE, the parent
#'   node is selected if tied values occur on the parent node and some of the
#'   chidren nodes.
#' @param message a logical value, TRUE or FALSE. The default is FALSE. If TRUE,
#'   the progress message is printed out.
#'
#' @importFrom methods is
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @return a value
#' @export
#' @author Ruizhu huang
#' @examples
#'
#' library(TreeSummarizedExperiment)
#' library(ggtree)
#' data("tinyTree")
#' ggtree(tinyTree, size = 2) +
#' geom_text2(aes(label = node), color = "darkblue",
#'            hjust = -0.5, vjust = 0.7, size = 4) +
#'     geom_text2(aes(label = label), color = "darkorange",
#'                hjust = -0.1, vjust = -0.7, size = 4)
#' set.seed(1)
#' count <- matrix(rnbinom(300,size=1,mu=10),nrow=10)
#' colnames(count) <- paste(rep(LETTERS[1:3], each = 10),
#'                          rep(1:10,3), sep = "_")
#' rownames(count) <- tinyTree$tip.label
#' count[1, ] <- 0
#' rowInf <- DataFrame(var1 = sample(letters[1:3], 10, replace = TRUE),
#'                     var2 = sample(c(TRUE, FALSE), 10, replace = TRUE))
#' colInf <- DataFrame(gg = factor(sample(1:3, 30, replace = TRUE)),
#'                     group = rep(LETTERS[1:3], each = 10))
#' lse <- TreeSummarizedExperiment(assays = list(count),
#'                                 rowData = rowInf,
#'                                 colData = colInf,
#'                                 rowTree = tinyTree)
#'
#' tse <- aggData(x = lse, onRow = TRUE)
#'
#' out <- runEdgeR(tse = tse, onRow = TRUE)
#' 
#' final <- searchLevel(tree = rowTree(tse),
#'                      scoreData = data.frame(out),
#'                      filter =  FDR > 0.05,
#'                      scoreCol = "PValue",
#'                      nodeCol = "nodeNum",
#'                      searchMax = FALSE,
#'                      parentPrior = TRUE, 
#'                      message = FALSE)
#'                      
#' # No node is kept because there is no difference between two groups
#' final$nodeNum[final$keep]
#' 

searchLevel <- function(tree, scoreData, filter, scoreCol,
                        nodeCol,
                        searchMax, parentPrior = TRUE,
                        message = FALSE) {

    # ------------ check inputs -----------------------------
    if (!is(tree, "phylo")) {
        stop("tree requires a phylo object.")
    }

    if (!is(scoreData, "data.frame")) {
        if (is(scoreData, "DataFrame")) {
            scoreData <- data.frame(scoreData)
        } else {
            stop("scoreData requires a data frame.")
        }
    }

    # add a column to store the result of search: keep
    if (any(colnames(scoreData) == "keep")) {
        stop("The result will be output in the 'keep' column;
             Please use other name for the current 'keep' column.")
    }

    if (message) {
        message("Preparing data ... ")
    }

    # ------------ fiter nodes -----------------------------
    scoreData$keep <- TRUE

    # filter nodes using the specified standard
    if (message) {
        message("filtering nodes ... ")
    }
    if (missing(filter)) {
        r <- FALSE
    } else {
        e <- substitute(filter)
        r <- eval(e, scoreData, parent.frame())

        if (!is.logical(r))
            stop("'filter' must be or evaluate to logical")
        r <- r & !is.na(r)
    }
    scoreData$keep[r] <- FALSE
    scoreData$keep[is.na(scoreData[, scoreCol])] <- FALSE
    if (message) {
        message(sum(!scoreData$keep), " nodes are filtered out...")
    }
    # # ------------ search nodes -----------------------------
    # This works only on nodes that are not filtered out.
    # It compares the score on an internal node with those on the descendant
    # nodes.
    # 1. searchMax = TRUE & parentPrior = TRUE
    #   An internal node is selected if its score is greater or equal to the max
    #   score of its descendants.
    # 2. searchMax = TRUE & parentPrior = FALSE
    #   An internal node is selected if its score is greater than the max score
    #   of its descendants.
    # 3. searchMax = FALSE & parentPrior = TRUE
    #   An internal node is selected if its score is lower or equal to the min
    #   score of its descendants.
    # 4. searchMax = FALSE & parentPrior = FALSE
    #   An internal node is selected if its score is lower than the min score of
    #   its descendants.

    if (message) {
        message("searching candidate nodes... ")
    }

    treeData <- printNode(tree = tree, type = "all")
    nodeI <- treeData$nodeNum[!treeData$isLeaf & scoreData$keep]
    if (message) {
        message("searching the descendant nodes of the candidate nodes... ")
    }
    descI <- findOS(tree = tree, node = nodeI, only.leaf = FALSE,
                    self.include = FALSE, use.alias = TRUE)

    if (message) {
        message("comparing nodes ... ")
        if (length(descI)) {
        pb <- txtProgressBar(0, length(descI), style = 3)
        }
    }

    for (i in seq_along(descI)) {
        node.i <- nodeI[i]
        desc.i <- descI[[i]]

        row.p <- match(node.i, scoreData[, nodeCol])
        row.d <- match(desc.i, scoreData[, nodeCol])

        # extract the values on the parent and on the children
        score.p <- scoreData[row.p, scoreCol]
        score.d <- scoreData[row.d, scoreCol]

        if (searchMax) {
            score.p[is.na(score.p)] <- -Inf
            score.d[is.na(score.d)] <- -Inf
            if (parentPrior) {
                isKeep <- all(score.p >= score.d)
            } else {
                isKeep <- all(score.p > score.d)
            }
        } else {
            score.p[is.na(score.p)] <- Inf
            score.d[is.na(score.d)] <- Inf
            if (parentPrior) {
                isKeep <- all(score.p <= score.d)
            } else {
                isKeep <- all(score.p < score.d)
            }
        }

        if (isKeep) {
            scoreData$keep[row.d] <- FALSE
        } else {
            scoreData$keep[row.p] <- FALSE
        }
        if (message) {
            setTxtProgressBar(pb, i)
            # message(i, " out of ", length(nodeI),
            #         " finished", "\r", appendLF = FALSE)
            # flush.console()
        }

        }

    return(scoreData)

}


