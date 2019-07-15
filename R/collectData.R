#' collect data for each node of the tree
#'
#' \code{collectData} generates data for each node of the tree by collecting
#' data on its descendant leaves. For leaf node, the collected data is from
#' itself. Users could decide on which dimension (row or column) the collection
#' should perform.
#'
#' @param x A TreeSummarizedExperiment object.
#' @param onRow A logical value, TRUE or FALSE.
#' @param message A logcical value, TRUE or FALSE.
#' @import TreeSummarizedExperiment
#' @importFrom methods is
#' @export
#'
#' @return A TreeSummarizedExperiment objects
#' @author Ruizhu HUANG
#'
#' @examples
#' #' 
#' library(ape)
#' library(ggtree)
#' set.seed(1)
#' 
#' # tree
#' exTree <- rtree(5)
#' ggtree(exTree) + geom_text2(aes(label = node), hjust = -0.5)
#' 
#' # toy table
#' toyTable <- matrix(1:25, nrow = 5)
#' rownames(toyTable) <- exTree$tip.label
#' 
#' # row data
#' rowD <- DataFrame(letters = sample(letters[1:2], 5, replace = TRUE))
#' # tse
#' lse <- TreeSummarizedExperiment(assays = list(toyTable),
#'                                 rowTree = exTree,
#'                                 rowData = rowD)
#' 
#' # collect data for each node
#' tse <- collectData(x = lse)
#' 
#' 
#' 

collectData <- function(x, rowLevel = NULL, colLevel = NULL,
                        message = FALSE) {
    
    # The input data should be on the leaf level
    if (!is.null(rowLevel)) {
        
        linkD <- rowLinks(x)
        treeD <- rowTree(x)
        if (is.character(rowLevel)) {
            level <- transNode(tree = treeD, node = rowLevel)  
        } else {
            level <- rowLevel
        }
        
        # aggregation
        out <- aggValue(x = x, rowLevel = level, 
                        FUN = function(x){x}, 
                        message = message)
        
    } else {
        linkD <- colLinks(x)
        treeD <- colTree(x)
        if (is.character(colLevel)) {
            level <- transNode(tree = treeD, node = colLevel)  
        } else {
            level <- colLevel
        }
        # aggregation
        out <- aggValue(x = x, colLevel = level, 
                        FUN = function(x){x}, 
                        message = message)
        
    }
    
    # The descendants of specified level
    desd <- findOS(tree = treeD, node = level, only.leaf = TRUE, 
                   self.include = TRUE)
    
    # create additional columns in link data to show leaf nodes where the
    # data is from
    infLeaf <- !level %in% treeD$edge[, 1]
    linkLevel <- lapply(seq_along(level), FUN = function(i) {
        
        #ii <- match(desd[[i]], linkD$nodeNum)
        ii <- lapply(desd[[i]], FUN = function(x) {
            which(linkD$nodeNum == x)
        })
        ii <- unlist(ii)
        link.i <- linkD[ii, ]
        # linkX <- LinkDataFrame(nodeLab = link.i$nodeLab,
        #                        nodeLab_alias = link.i$nodeLab_alias,
        #                        nodeNum = link.i$nodeNum, 
        #                        isLeaf = rep(infLeaf[i], nrow(link.i)), 
        #                        levelNum = rep(level[i], nrow(link.i)), 
        #                        levelLab = transNode(tree = treeD,
        #                                             node = level[i], 
        #                                             use.alias = FALSE),
        #                        levelLab_alias = transNode(tree = treeD, 
        #                                                   node = level[i], 
        #                                                   use.alias = TRUE))
        linkX <- LinkDataFrame(nodeLab = transNode(tree = treeD,
                                                   node = level[i], 
                                                   use.alias = FALSE),
                               nodeNum = rep(level[i], nrow(link.i)),
                               isLeaf = rep(infLeaf[i], nrow(link.i)), 
                               nodeLab_alias = transNode(tree = treeD, 
                                                          node = level[i], 
                                                          use.alias = TRUE),
                               from_nodeLab = link.i$nodeLab,
                               from_labAlias = link.i$nodeLab_alias,
                               from_nodeNum = link.i$nodeNum)
        return(linkX)
    })
    
    linkND <- do.call(rbind, linkLevel)
    
    
    if (!is.null(rowLevel)) {
        outF <- BiocGenerics:::replaceSlots(out, 
                                            metadata = metadata(x),
                                            rowLinks = linkND)
    } else {
        outF <- BiocGenerics:::replaceSlots(out, 
                                            metadata = metadata(x),
                                            colLinks = linkND)
    }
    
   return(outF) 
}





