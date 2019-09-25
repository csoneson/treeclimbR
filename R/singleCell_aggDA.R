#' calculate the number of cells
#' \code{aggDA} calculate number of cells per node-sample combination.
#' 
#' @param d_sce A \code{SingleCellExperiment} object
#' @param tree A phylo tree.
#' @param sample_id The column name of sample id.
#' @param group_id The column name of group id.
#' @param cluster_id The column name of clusters.
#' @param FUN The function
#' 
#' @importFrom SummarizedExperiment assays colData
#' @import TreeSummarizedExperiment
#' @importFrom dplyr "%>%" select mutate group_by summarize add_row
#' @export
#' 
#' @return A TreeSummarizedExperiment object. In matrix of \code{assays},
#'   samples in columns and clusters (or nodes) in rows.
#' 
#' @examples 
#' 



aggDA <- function(d_sce, tree, 
                  sample_id = "sample_id",
                  group_id = "group_id",
                  cluster_id = "cluster_id",
                  FUN = sum) {
    # cell information
    cell_info <- colData(d_sce)[, c(sample_id, group_id, cluster_id)]
    colnames(cell_info) <- c("sample_id", "group_id", "cluster_id")
    
    # the number of cells in each combinaton of sample_id and cluster_id 
    count_info <- cell_info %>%
        data.frame() %>%
        group_by(sample_id, cluster_id) %>%
        summarize(count = n())
    
    # split by samples:
    # clusters in rows, samples in columns
    list1 <- split(count_info, f = count_info$sample_id)
    
    # add rows for missing clusters
    nodes <- showNode(tree = tree, only.leaf = TRUE)
    clusters <- factor(names(nodes), levels = names(nodes))
    
    list2 <- lapply(list1, FUN = function(x) {
        lv <- clusters
        lvx <- unique(x$cluster_id)
        
        lvd <- setdiff(lv, lvx)
        if (length(lvd)) {
            xx <- x %>%
                ungroup() %>%
                add_row(sample_id = unique(x$sample_id),
                        cluster_id = lvd,
                        count = 0)
        } else {
            xx <- x
        }
        
        xx <- xx %>%
            mutate(cluster_id = factor(cluster_id, 
                                       levels = levels(clusters))) %>%
            arrange(cluster_id)
        return(xx)
    })
    
    # cluster lab
    lab <- list2[[1]]$cluster_id
    
    # sample info
    sample_info <- cell_info %>%
        data.frame() %>%
        select(sample_id, group_id) %>%
        distinct()
    
    # count table
    count <- lapply(list2, FUN = function(x) {
        xx <- cbind(x$count)
        colnames(xx) <- unique(x$sample_id)
        return(xx)
    })
    count <- do.call(cbind, count)
    count <- count[, sample_info$sample_id, drop = FALSE]
    rownames(count) <- lab
    
    # TreeSE
    lse <- TreeSummarizedExperiment(assays = list(count),
                                    rowTree = tree)
    tse <- aggValue(x = lse, rowLevel = nodes, FUN = FUN)
    return(tse)
    
}