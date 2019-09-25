#' build a tree from clusters 
#' \code{buildTree2} creates a tree (\code{phylo} class). 
#' 
#' Distances among clusters are calculated based on the median values of
#' features (genes) of clusters in reduced dimensions (PCA)
#' 
#' @param d_se A
#'   \code{\link[TreeSummarizedExperiment]{TreeSummarizedExperiment}} object
#'   output from \code{\link{prepareData}}.
#' @param column_cluster A vector of cluster labels. It should have the same
#'   length as the number of cells in \code{data_HVG}.
#' @param npcs A numeric value. The number of principal components used to
#'   calculated distances among clusters. Default (NULL) takes the minumum value
#'   of the number of clusters and the number of highly variable genes.
#' @param dist_method the distance measure to be used. This must be one of
#'   "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski".
#'   Any unambiguous substring can be given. Please refer to \code{method} in
#'   \code{\link[stats]{dist}}
#' @param hclust_method the agglomeration method to be used. This should be (an
#'   unambiguous abbreviation of) one of "ward.D", "ward.D2", "single",
#'   "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC)
#'   or "centroid" (= UPGMC). Please refer to \code{method} in
#'   \code{\link[stats]{hclust}}
#' 
#' @importFrom dplyr "%>%" mutate group_by summarize_all
#' @importFrom tibble column_to_rownames
#' @importFrom stats hclust dist prcomp
#' @importFrom ape as.phylo
#' @import TreeSummarizedExperiment
#' @export
#' @return A phylo object
#' @examples 
#' 
#' library(ggtree)
#' # data
#' ncell <- 500
#' ngene <- 20
#' dat <- matrix(rpois(ngene*ncell, lambda = 20), nrow = ngene)
#' 
#' # clusters
#' clus <- sample(1:50, 500, replace = TRUE)
#' 
#' # tree
#' tr <- buildTree2(data_HVG = dat, cluster_id = clus)
#' 
#' ggtree(tr, branch.length = "none") +
#'    geom_tiplab(size = 2) +
#'    xlim(0, 20)

# buildTree2 <- function(d_se,
#                        column_cluster,
#                        npcs = 10,
#                        dist_method = "euclidean",
#                        hclust_method = "average") {
#     
#     # scaled data on highly variable genes
#     data_HVG <- assays(d_se)$scale.data
#     # data_HVG <- metadata(d_se)$hvg_scaleData
#     
#     # cluster_id
#     cluster_id <- colData(d_se)[[column_cluster]]
#     
#     # median values of HVG on each cluster
#     med <- t(data_HVG) %>%
#         data.frame() %>%
#         mutate(cluster_id = cluster_id) %>%
#         group_by(cluster_id) %>%
#         summarize_all(median, na.rm = TRUE) %>%
#        # mutate(cluster_id = paste("cluster", cluster_id, sep = "_")) %>%
#         column_to_rownames("cluster_id")
#     
#     ll <- lapply(med, FUN = function(x) {
#         length(unique(x))
#     })
#     ll <- unlist(ll)
#     isD <- ll == 1
#     if(sum(isD)) {
#         warning("Remove ", sum(isD), 
#                 " gene(s) with the same median value in all clusters")
#         med <- med[, !isD]
#     }
#     
#     
#     # run pca over the median value of HVG on each cluster
#     med_pca <- prcomp(med, scale = TRUE)
#     x_pca <- med_pca$x
#     
#     
#     # calculate distance 
#     if (is.null(npcs)) {
#         npcs <- length(unique(cluster_id))
#     }
#     dst <- dist(x = x_pca[, seq_len(npcs)], method = dist_method)
#     tree_p <- hclust(d = dst, method = hclust_method)
#     tree_p <- as.phylo(tree_p)
#     
#     
#     return(tree_p)
# }

buildTree2 <- function(d_se,
                       column_cluster,
                       npcs = 10,
                       dist_method = "euclidean",
                       hclust_method = "average") {
    
    # scaled data on highly variable genes
    data_HVG <- metadata(d_se)$hvg_scaleData
    
    # data_HVG <- metadata(d_se)$hvg_scaleData
    
    # cluster_id
    cluster_id <- colData(d_se)[[column_cluster]]
    
    # median values of HVG on each cluster
    med <- t(data_HVG) %>%
        data.frame() %>%
        mutate(cluster_id = cluster_id) %>%
        group_by(cluster_id) %>%
        summarize_all(median, na.rm = TRUE) %>%
        # mutate(cluster_id = paste("cluster", cluster_id, sep = "_")) %>%
        column_to_rownames("cluster_id")
    
    ll <- lapply(med, FUN = function(x) {
        length(unique(x))
    })
    ll <- unlist(ll)
    isD <- ll == 1
    if(sum(isD)) {
        warning("Remove ", sum(isD), 
                " gene(s) with the same median value in all clusters")
        med <- med[, !isD]
    }
    
    
    # run pca over the median value of HVG on each cluster
    med_pca <- prcomp(med, scale = TRUE)
    x_pca <- med_pca$x
    
    
    # calculate distance 
    if (is.null(npcs)) {
        npcs <- length(unique(cluster_id))
    }
    dst <- dist(x = x_pca[, seq_len(npcs)], method = dist_method)
    tree_p <- hclust(d = dst, method = hclust_method)
    tree_p <- as.phylo(tree_p)
    
    
    return(tree_p)
}
