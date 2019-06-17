#' convert a TreeSummarizedExperiment object to a data frame
#' 
#' \code{convertDf} converts a TreeSummarizedExperiment object to a data frame.
#' The data frame includes observed values of entities on each sample and the
#' sample information that is originally stored in the column data (if samples
#' are in the column dimension of \code{assays} table) or the row data (if
#' samples are in the row dimension of \code{assays} table).
#' 
#' @param tse A TreeSummarizedExperiment object.
#' @param onRow A logical value, TRUE or FALSE. If entities are in rows and
#'   samples are in columns, then TRUE is required; if the other way around,
#'   then FALSE is required.
#' @param assayNum A numeric value to specify which table in the \code{assays}
#'   slot should be used. The first table is used by default.
#' 
#' @importFrom SummarizedExperiment assays
#' @import TreeSummarizedExperiment
#' @importFrom dplyr mutate "%>%" mutate_if
#' @importFrom tidyr gather
#' 
#' @export
#' @return A data frame
#' @author Ruizhu Huang
#' @examples   
#' 
#' library(ape)
#' library(TreeSummarizedExperiment)
#' library(ggtree)
#' 
#' set.seed(1)
#' treeR <- rtree(5)
#'
#' 
#' ggtree(treeR, size = 2) +
#'     geom_text2(aes(label = node), color = "darkblue",
#'                hjust = -0.5, vjust = 0.7, size = 6) +
#'     geom_text2(aes(label = label), color = "darkorange",
#'                hjust = -0.1, vjust = -0.7, size = 6)
#'                
#' toyTable <- matrix(rnbinom(20, size = 1, mu = 10), nrow = 5)
#' colnames(toyTable) <- paste(rep(LETTERS[1:2], each = 2),
#'                             rep(1:2, 2), sep = "_")
#' rownames(toyTable) <- treeR$tip.label
#'
#' toyTable
#'
#' # the column data
#' colInf <- DataFrame(gg = c(1, 2, 3, 3),
#'                     gr = rep(LETTERS[1:2], each = 2),
#'                     row.names = colnames(toyTable))
#' colInf
#'
#' # the toy tree

#'
#' 
#' tse <- TreeSummarizedExperiment(assays = list(toyTable),
#'                                 colData = colInf,
#'                                 rowTree = treeR)
#'
#' df <- convertDf(tse = tse, onRow = TRUE)
#'  

convertDf <- function(tse, onRow = TRUE, 
                   assayNum = NULL) {
    
    # The input data should be a TreeSummarizedExperiment object
    # It should includes only data on the leaf nodes
    if (!is(tse, "TreeSummarizedExperiment")) {
        stop("tse should be a TreeSummarizedExperiment.")
    }
    
    # If not specified, the first table in the assays is used
    if (is.null(assayNum)) {
        assayNum <- 1
    }
    
    
    # analysis step
    if (onRow) {
        # the count table
        count <- assays(tse)[[assayNum]]
        
        # the link data
        ld <- rowLinks(tse)
        
        # the annotated data for samples
        spd <- data.frame(colData(tse))
        
        indDim <- "row"
    } else {
        # the count table
        count <- t(assays(tse)[[assayNum]])
        
        # extract link data
        ld <- colLinks(tse)
        
        # the annotated data for samples
        spd <- data.frame(rowData(tse))
        
        indDim <- "col"
       
    }
    addC <- "int_add_col"
    reqName <- c(addC, "node", "level", "isLeaf", "sample", "value")
    if (any(colnames(spd) %in% reqName)) {
        stop("Columns: (", 
             paste(reqName, collapse = ","), ") will be generated in output;", 
             "please use other column names in the ", indDim, "Data of tse")
    }
    
    spd[[addC]] <- colnames(count)
    # prepare data in the right form to run logistic model or logistic mixed
    # effects model
    
    # prepare data in the right form to run logistic model or logistic mixed
    # effects model
    df <- data.frame(count, check.names = FALSE) %>%
        mutate(node = ld$nodeNum) %>%
        mutate(level = ld$levelNum) %>%
        mutate(isLeaf = ld$isLeaf) %>%
        gather(sample, value, colnames(count), factor_key = FALSE) %>%
        mutate_if(is.factor, as.character)
    if (is.null(df$level)) {
        df <- df %>%
            mutate(level = node) 
    }
    
    
    sampInf <- spd[match(df$sample, spd[[addC]]), 
                   setdiff(colnames(spd), addC), drop = FALSE]
    rownames(sampInf) <- rownames(df) 
    dff <- cbind(df, sampInf)
        
    
    return(dff)
}
