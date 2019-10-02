#' prepare data
#' 
#' \code{prepareData} is a wrapper function that includes functions
#' \code{\link[Seurat]{RunPCA}}, \code{\link[Seurat]{FindNeighbours}}, and
#' \code{\link[Seurat]{FindClusters}}. It clusters cells at specified
#' resolutions, and outputs  a
#' \code{\link[SingleCellExperiment]{SingleCellExperiment}} object. In the
#' output, the \code{assay} stores the original counts that was stored in the
#' slot \code{d_seurat@assays$RNA@counts}. The cluster results is stored in the
#' column with prefix \code{reso_} followed by the specified \code{reso}. The
#' integrated data is stored in \code{metadata}
#'  
#' @param d_seurat A \code{Seurat} object to provide data that is normalized,
#'   scaled and integrated from package \code{Seurat}. Marker genes of cell
#'   types are selected as the integrated features and their scale data will be
#'   used in next steop to generate the tree.
#' @param assay The name of assay that will be used to run
#'   \code{\link[Seurat]{runPCA}} and perform clustering.
#' @param npcs Please refer to \code{npcs} from \code{\link[Seurat]{RunPCA}}.
#' @param dims Please refer to \code{dims} from
#'   \code{\link[Seurat]{FindNeighbours}}.
#' @param reso A numeric vector. It provides the resolutions that the clustering
#'   is generated at. Please refer to \code{resolution} from
#'   \code{\link[Seurat]{FindClusters}}.
#' @param seed A numeric value. It provides a seed number. Please refer to
#'   \code{random.seed} from \code{\link[Seurat]{FindClusters}}.
#'
#' @importFrom Seurat RunPCA FindNeighbors FindClusters
#' @import TreeSummarizedExperiment
#' 
#' @export
#' @author Ruizhu Huang
#' @examples 
#' 
#' reso <- 1:2
#' 

# prepareData <- function(d_seurat, reso, seed, 
#                         npcs = 30, dims = 1:9) {
#     
#     # run PCA
#     d_seurat <- RunPCA(d_seurat, npcs = npcs, verbose = FALSE)
#     
#     # find neighbours
#     d_seurat <- FindNeighbors(d_seurat, reduction = "pca", 
#                               dims = dims, verbose = FALSE)
#     # clustering
#     for (r in reso) {
#         d_seurat <- FindClusters(d_seurat, resolution = r,
#                                  random.seed = seed, verbose = FALSE)
#     }
#     
#     
#     
#     # the cell information
#     cell_info <- d_seurat@meta.data
#     colnames(cell_info) <- gsub(
#         pattern = "integrated_snn_res.",
#         replacement = "res_", 
#         colnames(cell_info))
#     
#     # the highly variable gene (HVG)
#     hvg <- d_seurat@assays$integrated@var.features
#     
#     # the normalized and scaled data of HVG
#     hvg_s <- d_seurat@assays$integrated@scale.data
#     
#     # original counts for HVG
#     hvg_o <- d_seurat@assays$RNA@counts[hvg, ]
#     
#     
#     tse <- TreeSummarizedExperiment(
#         assays = list(counts = hvg_o, scale.data = hvg_s),
#         colData = cell_info)
#     
#     
#     # counts of all genes
#     # g_o <- d_seurat@assays$RNA@counts
#     # g_s <- d_seurat@assays$RNA@scale.data
#     # output sce
#     
#     
#     # tse <- TreeSummarizedExperiment(
#     #     assays = list(counts = g_o),
#     #     colData = cell_info, 
#     #     metadata = list(hvg = hvg, hvg_scaleData = hvg_s, hvg_counts = hvg_o))
#     return(tse)
#     
# }

prepareData <- function(d_seurat, reso, seed, 
                        assay = "integrated", 
                        npcs = 30, dims = 1:9) {
    
    # run PCA
    d_seurat <- RunPCA(d_seurat, assay = assay, npcs = npcs, 
                       verbose = FALSE)
    
    # find neighbours
    d_seurat <- FindNeighbors(d_seurat, reduction = "pca", 
                              dims = dims, verbose = FALSE)
    # clustering
    for (r in reso) {
        d_seurat <- FindClusters(d_seurat, resolution = r,
                                 random.seed = seed, verbose = FALSE)
    }
    
    
    
    # the cell information
    cell_info <- d_seurat@meta.data
    colnames(cell_info) <- gsub(
        pattern = "integrated_snn_res.",
        replacement = "reso_", 
        colnames(cell_info))
    
    # the integrated features
    mk <- d_seurat@assays$integrated@var.features
    
    # the normalized and scaled data of integrated features
    data_mk <- d_seurat@assays$integrated@scale.data
    
    # original counts for HVG
    # data_o <- d_seurat@assays$RNA@counts[mk, ]
    
    
    # tse <- TreeSummarizedExperiment(
    #     assays = list(counts = hvg_o, scale.data = hvg_s),
    #     colData = cell_info)
    
    
    # counts of all genes
    g_o <- d_seurat@assays$RNA@counts
    #g_s <- d_seurat@assays$RNA@scale.data
    
    # output sce
    tse <- TreeSummarizedExperiment(
        assays = list(counts = g_o),
        colData = cell_info,
        metadata = list(type_marker = mk, 
                        data_type_marker = data_mk))
    return(tse)
    
}
