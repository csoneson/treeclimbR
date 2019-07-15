#' search for a target level on the tree via a specified score
#'
#' \code{getLevel} searches for the target level of the tree via a specified
#' score. The score value needs to be provided for each node of the tree.
#'
#' @param score_data a data frame to provide scores for all nodes of the tree.
#'   The data frame should have at least five columns, four of them naming
#'   \code{nodeLab}, \code{nodeLab_alias}, \code{nodeNum} and \code{isLeaf}.
#'   These four columns provide the information about nodes. The score column
#'   could be named freely.
#' @param tree a phylo object to provide the hierarchical structure of nodes.
#' @param drop a logical expression indicating elements or rows to keep:
#'   missing values are taken as false.
#' @param score_column the name of the score column. The search is based on the
#'   score value.
#' @param node_column the name of the column that stores the node number.
#' @param get_max a logical value, TRUE or FALSE. If TRUE, search for nodes
#'   that has higher score value than its descendants; otherwise, search for
#'   nodes that has lower score value than its descendants.
#' @param parent_first a logical value, TRUE or FALSE. If TRUE, the parent
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
#' final <- getLevel(tree = rowTree(tse),
#'                      score_data = data.frame(out),
#'                      drop =  FDR > 0.05,
#'                      score_column = "PValue",
#'                      node_column = "nodeNum",
#'                      get_max = FALSE,
#'                      parent_first = TRUE, 
#'                      message = FALSE)
#'                      
#' # No node is kept because there is no difference between two groups
#' final$nodeNum[final$keep]
#' 

getLevel <- function(tree, score_data, drop, score_column,
                     node_column, get_max, parent_first = TRUE,
                     message = FALSE) {
    
    # ------------ check inputs -----------------------------
    if (!is(tree, "phylo")) {
        stop("tree requires a phylo object.")
    }
    
    if (!is(score_data, "data.frame")) {
        if (is(score_data, "DataFrame")) {
            score_data <- data.frame(score_data)
        } else {
            stop("score_data requires a data frame.")
        }
    }
    
    # add a column to store the result of search: keep
    if (any(colnames(score_data) == "keep")) {
        stop("The result will be output in the 'keep' column;
             Please use other name for the current 'keep' column.")
    }
    
    if (message) {
        message("Preparing data ... ")
    }
    
    # ------------ fiter nodes -----------------------------
    score_data$keep <- TRUE
    
    # drop nodes using the specified standard
    if (message) {
        message("dropping nodes ... ")
    }
    if (missing(drop)) {
        r <- FALSE
    } else {
        e <- substitute(drop)
        r <- eval(e, score_data, parent.frame())
        
        if (!is.logical(r))
            stop("'drop' must be or evaluate to logical")
        r <- r & !is.na(r)
    }
    score_data$keep[r] <- FALSE
    score_data$keep[is.na(score_data[, score_column])] <- FALSE
    if (message) {
        message(sum(!score_data$keep), " nodes are ruled out...")
    }
    # # ------------ search nodes -----------------------------
    # This works only on nodes that are not ruled out.
    # It compares the score on an internal node with those on the descendant
    # nodes.
    # 1. get_max = TRUE & parent_first = TRUE
    #   An internal node is selected if its score is greater or equal to the max
    #   score of its descendants.
    # 2. get_max = TRUE & parent_first = FALSE
    #   An internal node is selected if its score is greater than the max score
    #   of its descendants.
    # 3. get_max = FALSE & parent_first = TRUE
    #   An internal node is selected if its score is lower or equal to the min
    #   score of its descendants.
    # 4. get_max = FALSE & parent_first = FALSE
    #   An internal node is selected if its score is lower than the min score of
    #   its descendants.
    
    if (message) {
        message("searching candidate nodes... ")
    }
    
    treeData <- printNode(tree = tree, type = "all")
    nodeI <- treeData$nodeNum[!treeData$isLeaf & score_data$keep]
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
        
        row.p <- match(node.i, score_data[, node_column])
        row.d <- match(desc.i, score_data[, node_column])
        
        # extract the values on the parent and on the children
        score.p <- score_data[row.p, score_column]
        score.d <- score_data[row.d, score_column]
        
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
        if (message) {
            setTxtProgressBar(pb, i)
            # message(i, " out of ", length(nodeI),
            #         " finished", "\r", appendLF = FALSE)
            # flush.console()
        }
        
    }
    
    return(score_data)
    
    }


