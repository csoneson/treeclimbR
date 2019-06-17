#' create a heatmap
#' 
#' create a heatmap
#' 
#' @param obj A TreeSummarizedExperiment
#' @param scaleColumn TRUE or FALSE. If TRUE, the column of assay table is
#'   scaled.
#' @param groupCol The name of the column that contains the group information of
#'   samples in \code{colData} of \code{obj}.
#' @param ... Arguments that are accepted in
#'   \code{\link[ComplexHeatmap]{Heatmap}}.
#' @export
#' @importFrom S4Vectors metadata
#' @import TreeSummarizedExperiment
#' @importFrom ggtree fortify
#' @importFrom ComplexHeatmap Heatmap rowAnnotation
#' @export
#' @return a heatmap
#' 

viewHeatmap <- function(obj, groupCol = "group", 
                        scaleColumn = FALSE, ...) {
    # category
    fc <- metadata(obj)$FC
    cg <- list("increase" = names(fc[fc > 1]), 
               "decrease" = names(fc[fc < 1]),
               "other" = names(fc[fc == 1]))
    
    # data on the leaf level
    rLink <- rowLinks(obj)
    obj <- obj[rLink$isLeaf, ]
    
    # counts
    count <- assays(obj)[[1]]
    rownames(count) <- rowLinks(obj)$nodeLab
    
    # tree
    rTree <- rowTree(obj)
    
    # scenario
    scene <- metadata(obj)$scenario
    
    # group
    gr <- colData(obj)[groupCol]
    
    hp <- .viewHP(tree = rTree, count = count, 
                  category = cg, scaleColumn = scaleColumn,
                  group = gr, name = scene, 
                  ...)
    draw(hp)
}
.viewHP <- function(tree, count, category, 
                    scaleColumn, group, ...) {
    
    d <- fortify(tree)
    td <- subset(d, isTip)
    leafLab <- with(td, label[order(y, decreasing=T)])
     
    
    # treeDr <- as.dendrogram(force.ultrametric(tree))
    
    
    sdm <- apply(count, 2, FUN = function(x){
        if (scaleColumn) {
            scale(x)
        } else {x}
    })
    rownames(sdm) <- rownames(count)
    cn <- colnames(sdm)
    
    # ht <- Heatmap(sdm, cluster_columns = FALSE, cluster_rows = treeDr, 
    #               show_column_names = FALSE)
    ht <- Heatmap(sdm, row_order = leafLab, 
                  cluster_columns = FALSE, 
                  show_column_names = FALSE,
                  show_row_names = FALSE,
                  column_split = group, ...)
    # oo <- row_order(ht)
    # rowNam <- rownames(count)[oo]
    # a data frame to store category information
    df <- lapply(seq_along(category), FUN = function(x) {
        xx <- data.frame(name = category[[x]],
                         category = names(category)[[x]], 
                         stringsAsFactors = FALSE)
        return(xx)
    })
    df <- do.call(rbind, df)
    df <- df[match(rownames(sdm), df$name), ]
    
    ha <- rowAnnotation(change = df$category, 
                        col = list(change = c("decrease" = "blue", 
                                              "other" = "grey",
                                              "increase" = "orange")))
                                              
    out <- ha + ht
    return(out)
}


# set.seed(123)
# mat = matrix(rnorm(100), 10)
# rownames(mat) = paste0("R", 1:10)
# colnames(mat) = paste0("C", 1:10)
# column_ha = HeatmapAnnotation(foo1 = runif(10), bar1 = anno_barplot(runif(10)))
# row_ha = rowAnnotation(foo2 = runif(10), bar2 = anno_barplot(1:10))
# Heatmap(mat, name = "mat", top_annotation = column_ha, right_annotation = row_ha)
# h <- Heatmap(mat)
# h + row_ha
# draw(row_ha)
# 
# 
# library(ggplot2)
# library(ggtree)
# set.seed(102)
# 
# p <- ggtree(tinyTree)
# cc <- count
# gheatmap(p, cc) + scal
# 
# sim_tree <- rtree(48, rooted = FALSE)
# 
# ggsim <- ggtree(sim_tree)
# 
# df <- data.frame("category" = rep(c("a", "b", "c", "d"), 3))
# rownames(df) = sample(sim_tree$tip.label, replace=F, 12)
# gheatmap(ggsim, df[, "category", drop = FALSE], width = 0.1, color = "black") 
# 
# set.seed(1)
# exTree <- ape::rtree(10)
# 
# # visualize the tree
# ggtree(exTree, branch.length = "none") + geom_tiplab() 
# 
# count <- matrix(1:30, ncol = 3)
# rownames(count) <- exTree$tip.label
# 
# newTree <- as.dendrogram(phytools::force.ultrametric(exTree))
# Heatmap(count, cluster_rows = newTree)
# 
# # order the rows of matrix based on the leaf labels of the original tree
# d <- fortify(exTree)
# td <- subset(d, isTip)
# leafLab <- with(td, label[order(y, decreasing=T)])
# Heatmap(count, row_order = leafLab)
# 
# 
# 
# 
