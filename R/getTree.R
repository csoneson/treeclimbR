#' build a tree from clusters 
#' \code{getTree} creates a tree (\code{phylo} class). 
#' 
#' Distances among clusters are calculated based on the median values of
#' features (genes) of clusters in reduced dimensions (PCA)
#' 
#' @param SE A
#'   \code{\link[SummarizedExperiment]{SummarizedExperiment}} object
#'   output from \code{\link{medianByClusterMarker}}.
#' @param column_cluster The name of the column that stores the cluster id.
#' @param use_marker A logical vector that has the same length as the number of
#'   columns in \code{x} to specify markers to build the tree.
#' @param npcs A numeric value. The number of principal components used to
#'   calculated distances among clusters. Default (NULL) takes the minumum value
#'   between the number of clusters and the number of markers.
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
#' @importFrom SummarizedExperiment rowData assays
#' @export
#' @return A phylo object
#' @examples 
#' 
#' library(ggtree)
#' 
#' # data
#' ncell <- 500
#' ngene <- 20
#' dat <- matrix(rpois(ngene*ncell, lambda = 20), ncol = ngene)
#' 
#' # clusters
#' clus <- sample(letters, 500, replace = TRUE)
#' rowD <- data.frame(cluster_id = clus,
#'                    stringsAsFactors = FALSE)
#' d_se <- SummarizedExperiment(assays = list(dat),
#'                              rowData = rowD)
#' dx <- medianByClusterMarker(SE = d_se, marker_in_column = TRUE,
#'                             column_cluster = "cluster_id")
#'              
#' # tree
#' tr <- getTree(SE = dx, column_cluster = "cluster_id",
#'               npcs = 20)
#' # leaf nodes are labelled with cluster id
#' tr$tip.label
#' ggtree(tr, branch.length = "none") +
#'    geom_tiplab(size = 2) +
#'    xlim(0, 20)

getTree <- function(SE,
                    column_cluster = "cluster_id",
                    use_marker = NULL,
                    npcs = NULL,
                    dist_method = "euclidean",
                    hclust_method = "average") {
    
    # decide the marker to be used
    if (!is.null(use_marker)) {
       SE <- SE[, use_marker] 
    }
    
    # cluster_id
    cluster_id <- rowData(SE)[[column_cluster]]
    
    # data of type markers
    med <- assays(SE)[[1]]
    rownames(med) <- cluster_id
    
    
    # remove markers that have the same value in all clusters
    ll <- apply(med, 2, FUN = function(x) {
        length(unique(x))
    })
    isD <- ll == 1
    if(sum(isD)) {
        warning("Remove ", sum(isD), 
                " gene(s) with the same median value in all clusters")
        med <- med[, !isD]
    }
    
    
    # run pca 
    med_pca <- prcomp(med, scale = TRUE)
    x_pca <- med_pca$x
    
    
    # calculate distance 
    if (is.null(npcs)) {
        npcs <- min(length(unique(cluster_id)), ncol(med))
    }
    
    # build tree
    dst <- dist(x = x_pca[, seq_len(npcs)], method = dist_method)
    tree_p <- hclust(d = dst, method = hclust_method)
    tree_p <- as.phylo(tree_p)
    
    # # output
    # 
    # out <- TreeSummarizedExperiment(assays = list(med),
    #                                 rowTree = tree_p, 
    #                                 rowNodeLab = cluster_id,
    #                                 rowData = rowData(x))
    
    return(tree_p)
}


