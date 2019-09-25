#' aggregate data
#' \code{aggDS} aggregate data based on the tree.
#' 
#' @param d_sce A \code{SingleCellExperiment} object
#' @param tree A phylo tree.
#' @param assay The name of a \code{assays}.
#' @param sample_id The column name of sample id.
#' @param group_id The column name of group id.
#' @param cluster_id The column name of clusters.
#' @param FUN The function
#' 
#' @importFrom SummarizedExperiment assays colData
#' @import TreeSummarizedExperiment
#' @importFrom dplyr "%>%" select mutate
#' @importFrom tibble column_to_rownames
#' @importFrom S4Vectors t split
#' @export
#' 
#' @return A singleCellExperiment object. In \code{assays}, each matrix is data of a node on the tree. Rows of matrix represent genes and columns represent samples.
#' 
#' @examples 
#' 

# aggDS <- function(d_sce, tree, 
#                   assay = "counts", 
#                   sample_id = "sample_id",
#                   group_id = "group_id",
#                   cluster_id = "cluster_id",
#                   FUN = sum, message = FALSE) {
#     
#     if (message) {
#         message("Preparing ...") }
#     # cell information
#     cell_info <- colData(d_sce)[, c(sample_id, group_id, cluster_id)]
#     colnames(cell_info) <- c("sample_id", "group_id", "cluster_id")
#     
#     
#     if (message) {
#         message("Splitting data by sample_id ...") }
#     # a list: each element is a sample
#     #         cells in rows, genes + cell_info in columns
#     asy <- assays(d_sce)[[assay]]
#     # asy <- S4Vectors::t(asy)
#     # dat <- base::cbind(cell_info, asy)
#     # dat_list <- S4Vectors::split(x = dat, f = dat$sample_id)
#     asy <- t(asy)
#     dat <- cbind(cell_info, asy)
#     dat_list <- split(x = dat, f = dat$sample_id)
#     
#     # a list: each element is a sample
#     #         nodes (or clusters) in rows, nodeNum + genes in columns
#     if (message) {
#         message("Generating data on nodes of the tree for each sample...") }
#     
#     node <- showNode(tree = tree, only.leaf = FALSE, use.alias = TRUE)
#     dat_list2 <- lapply(seq_along(dat_list), FUN = function(i, ff) {
#         
#         x <- dat_list[[i]]
#         ays_i <- x %>%
#             data.frame() %>%
#             select(-sample_id, -group_id, -cluster_id) %>%
#             as.matrix()
#         lab_i <- as.character(x$cluster_id)
#         lse_i <- TreeSummarizedExperiment(assays = list(counts = ays_i),
#                                           rowTree = tree, rowNodeLab = lab_i)
#         tse_i <- aggValue(x = lse_i, rowLevel = node, FUN = ff)
#         out_i <- assays(tse_i)[[1]]
#         out_ci <- out_i %>%
#             data.frame() %>%
#             mutate(node = rowLinks(tse_i)$nodeNum,
#                    sample_id = names(dat_list)[i])
#         return(out_ci)
#         
#         if (message) {
#             message(i, " out of ", length(dat_list),
#                     " samples finished", "\r", appendLF = FALSE)
#             flush.console()
#         }
# 
#     }, ff = FUN)
#     # system.time({
#     #     dat_test <- bind_rows(dat_list2)
#     # })
#     dat_2 <- do.call(rbind, dat_list2)
#     
#     # a list: each element is a node (or cluster)
#     #         samples in rows, genes in columns
#     dat_list3 <- split(x = dat_2, f = node)
#     names(dat_list3) <- names(node)
# 
#     # a list: each element is a node (or cluster)
#     #         genes in rows, samples in columns
#     
#     dat_list4 <- lapply(dat_list3, FUN = function(x) {
#         rownames(x) <- x$sample_id
#         xx <- x %>%
#             select(-node, -sample_id) %>%
#             as.matrix() %>%
#             t()
#         return(xx)
#     })
# 
#     
#     
#     # sample information
#     if (message) {
#         message("Working on sample data ...") }
#     sample_info <- colData(d_sce)[, c(sample_id, group_id)] %>%
#         data.frame() %>%
#         distinct() %>%
#         mutate(group_id = factor(group_id)) %>%
#         column_to_rownames(sample_id)
#     sample_info <- sample_info[colnames(dat_list4[[1]]), ,drop = FALSE]
# 
#     # meta data:
#     if (message) {
#         message("Working on metadata ...") }
#     cell_n <- table(colData(d_sce)[[sample_id]])
#     experiment_info <- colData(d_sce)[, c(sample_id, group_id)] %>%
#         data.frame() %>%
#         distinct() %>%
#         mutate(n_cells = cell_n[sample_id])
# 
#     agg_pars <- list(assay = assay,
#                      by = c(cluster_id, sample_id),
#                      fun = sum,
#                      scale = FALSE)
# 
#     cell_tab <- table(colData(d_sce)[[cluster_id]],
#                       colData(d_sce)[[sample_id]])
#     rn <- transNode(tree = tree,
#                     node = rownames(cell_tab))
# 
# 
#     desd <- findOS(tree = tree, node = node,
#                    only.leaf = TRUE, self.include = TRUE)
#     cell_tab <- as.matrix(cell_tab)
#     clus_tab <- lapply(desd, FUN = function(x) {
#         ii <- match(x, rn)
#         xi <- cell_tab[ii, , drop = FALSE]
#         apply(xi, 2, sum)
# 
#     })
#     clus_tab <- do.call(rbind, clus_tab)
#     rownames(clus_tab) <- names(node)
# 
#     meta <- list(experiment_info = experiment_info,
#                  agg_pars = agg_pars,
#                  n_cells = clus_tab )
#     # output
#     if (message) {
#         message("Output data ...") }
#     sce <- SingleCellExperiment(assays = dat_list4,
#                                 colData = sample_info,
#                                 metadata = meta)
#     return(sce)
#     
#     
# }


# aggDS <- function(d_sce, tree, 
#                   assay = "counts", 
#                   sample_id = "sample_id",
#                   group_id = "group_id",
#                   cluster_id = "cluster_id",
#                   FUN = sum, message = FALSE) {
#     
#     if (message) {
#         message("Preparing ...") }
#     
#     # cell information
#     cell_info <- colData(d_sce)[, c(sample_id, group_id, cluster_id)]
#     colnames(cell_info) <- c("sample_id", "group_id", "cluster_id")
#     cell_info$node <- transNode(tree = tree, 
#                                 node = as.character(cell_info$cluster_id))
#     
#     
#     if (message) {
#         message("Splitting data by sample_id ...") }
#     
#     # a list: each element is a sample
#     #         cells in rows, genes + cell_info in columns
#     asy <- assays(d_sce)[[assay]]
#     dat <- cbind(cell_info, t(asy))
#     dat_list <- split(x = dat, f = dat$sample_id)
#     
#     # a list: each element is a sample
#     #         nodes (or clusters) in rows, nodeNum + genes in columns
#     if (message) {
#         message("Generating data on nodes of the tree for each sample...") }
#     
#     node <- showNode(tree = tree, only.leaf = FALSE, use.alias = TRUE)
#     desd_list <- findOS(tree = tree, node = node, 
#                         only.leaf = TRUE, self.include = TRUE)
#     
#     
#     
#     dat_list2 <- lapply(seq_along(dat_list), FUN = function(i, ff) {
#         x <- dat_list[[i]]
#         ays_i <- x[, !colnames(x) %in% c(sample_id, group_id, cluster_id)]
#         lab_i <- as.character(x$cluster_id)
#         lse_i <- TreeSummarizedExperiment(assays = list(counts = ays_i),
#                                           rowTree = tree, rowNodeLab = lab_i)
#         tse_i <- aggValue(x = lse_i, rowLevel = node, FUN = ff)
#         out_i <- assays(tse_i)[[1]]
#         out_ci <- cbind(out_i,
#                         node = rowLinks(tse_i)$nodeNum)
#         
#         rownames(out_ci) <- rep(names(dat_list)[i], nrow(out_ci))
#         
#         if (message) {
#             message(i, " out of ", length(dat_list),
#                     " samples finished", "\r", appendLF = FALSE)
#             flush.console()
#         }
#         return(out_ci)
#         
#         
#         
#     }, ff = FUN)
#     # system.time({
#     #     dat_test <- bind_rows(dat_list2)
#     # })
#     dat_2 <- do.call(rbind, dat_list2)
#     
#     # the row indicator for each node
#     # a list: each element is data at a node
#     #         genes in rows, samples in columns
#     ri <- lapply(node, 
#                  FUN = function(x) {
#                      dat_2[, "node"] ==x })
#     dat_list3 <- lapply(ri, FUN = function(x) {
#         xx <- dat_2[x, !colnames(dat_2) %in% "node"]
#         t(xx)
#     })
#     
#     
#     
#     
#     # sample information
#     if (message) {
#         message("Working on sample data ...") }
#     sample_info <- colData(d_sce)[, c(sample_id, group_id)] %>%
#         data.frame() %>%
#         distinct() %>%
#         mutate(group_id = factor(group_id)) %>%
#         column_to_rownames(sample_id)
#     sample_info <- sample_info[colnames(dat_list3[[1]]), ,drop = FALSE]
#     
#     # meta data:
#     if (message) {
#         message("Working on metadata ...") }
#     cell_n <- table(colData(d_sce)[[sample_id]])
#     experiment_info <- colData(d_sce)[, c(sample_id, group_id)] %>%
#         data.frame() %>%
#         distinct() %>%
#         mutate(n_cells = cell_n[sample_id])
#     
#     agg_pars <- list(assay = assay,
#                      by = c(cluster_id, sample_id),
#                      fun = sum,
#                      scale = FALSE)
#     
#     cell_tab <- table(colData(d_sce)[[cluster_id]],
#                       colData(d_sce)[[sample_id]])
#     rn <- transNode(tree = tree,
#                     node = rownames(cell_tab))
#     
#     
#     desd <- findOS(tree = tree, node = node,
#                    only.leaf = TRUE, self.include = TRUE)
#     cell_tab <- as.matrix(cell_tab)
#     clus_tab <- lapply(desd, FUN = function(x) {
#         ii <- match(x, rn)
#         xi <- cell_tab[ii, , drop = FALSE]
#         apply(xi, 2, sum)
#         
#     })
#     clus_tab <- do.call(rbind, clus_tab)
#     rownames(clus_tab) <- names(node)
#     
#     meta <- list(experiment_info = experiment_info,
#                  agg_pars = agg_pars,
#                  n_cells = clus_tab )
#     # output
#     if (message) {
#         message("Output data ...") }
#     sce <- SingleCellExperiment(assays = dat_list3,
#                                 colData = sample_info,
#                                 metadata = meta)
#     return(sce)
#     
#     
# }



aggDS <- function(d_sce, tree, 
                  assay = "counts", 
                  sample_id = "sample_id",
                  group_id = "group_id",
                  cluster_id = "cluster_id",
                  FUN = sum, message = FALSE) {
    
    if (message) {
        message("Preparing ...") }
    
    # cell information
    if (message) {
        message("Extracting cell information ...") }
    cell_info <- colData(d_sce)[, c(sample_id, group_id, cluster_id)]
    colnames(cell_info) <- c("sample_id", "group_id", "cluster_id")
    cell_info$node <- transNode(tree = tree, 
                                node = as.character(cell_info$cluster_id))
    
    # tree information
    if (message) {
        message("Extracting tree information ...") }
    node <- showNode(tree = tree, only.leaf = FALSE, use.alias = TRUE)
    desd_list <- findOS(tree = tree, node = node, 
                        only.leaf = TRUE, self.include = TRUE)
    names(desd_list) <- names(node)
    
    ## indicators
    if (message) {
        message("Grouping cells by samples and nodes ...") }
    # an indicator: cells belong to each node of the tree
    ind_list1 <- lapply(desd_list, FUN = function(x) {
        which(cell_info$node %in% x)
    })
    
    # an indicator: cells belong to each sample
    sp <- unique(cell_info$sample_id)
    ind_list2 <- lapply(sp, FUN = function(x) {
        which(cell_info$sample_id %in% x)
    })
    names(ind_list2) <- sp
    
    # an indicator: cells to each combination of nodes and samples
    ind_list3 <- lapply(ind_list1, FUN = function(x) {
        xx <- lapply(ind_list2, intersect, y = x)
    })
    
    
    if (message) {
        message("Preparing data on each node ...") }
    ## reshape data
    asy <- assays(d_sce)[[assay]]
    dat_list <- vector("list", length(ind_list3))
    names(dat_list) <- names(ind_list3)
    for (i in seq_along(dat_list)) {
        
        ind_i <- ind_list3[[i]]
        sp_i <- names(ind_i)
        asy_i <- lapply(ind_i, FUN = function(x) {
            # genes in rows, cells in columns
            xx <- asy[, x, drop = FALSE]
            
            # aggregation
            ax <- apply(xx, 1, FUN = FUN)
        })
        asy_i <- do.call(cbind, asy_i)
        colnames(asy_i) <- sp_i
        dat_list[[i]] <- asy_i
        
        if (message) {
            message(i, " out of ", length(dat_list),
                    " nodes finished", "\r", appendLF = FALSE)
            flush.console()
        }
    }
    
    # sample information
    if (message) {
        message("Working on sample information ...") }
    sample_info <- colData(d_sce)[, c(sample_id, group_id)] %>%
        data.frame() %>%
        distinct() %>%
        mutate(group_id = factor(group_id)) %>%
        column_to_rownames(sample_id)
    sample_info <- sample_info[colnames(dat_list[[1]]), ,drop = FALSE]
    
    # meta data:
    if (message) {
        message("Working on metadata ...") }
    cell_n <- table(cell_info[[sample_id]])
    experiment_info <- cell_info %>%
        data.frame() %>%
        select(sample_id, group_id) %>%
        distinct() %>%
        mutate(n_cells = cell_n[sample_id])
    
    agg_pars <- list(assay = assay,
                     by = c(cluster_id, sample_id),
                     fun = deparse(substitute(FUN)),
                     scale = FALSE)
    # on leaf level
    cell_tab <- table(cell_info[[cluster_id]],
                      cell_info[[sample_id]])
    rn <- transNode(tree = tree,
                    node = rownames(cell_tab))
    
    # on all levels
    desd <- findOS(tree = tree, node = node,
                   only.leaf = TRUE, self.include = TRUE)
    cell_tab <- as.matrix(cell_tab)
    clus_tab <- lapply(desd, FUN = function(x) {
        ii <- match(x, rn)
        xi <- cell_tab[ii, , drop = FALSE]
        apply(xi, 2, sum)
        
    })
    clus_tab <- do.call(rbind, clus_tab)
    rownames(clus_tab) <- names(node)
    clus_tab <- clus_tab[, colnames(dat_list[[1]])]
    
    meta <- list(experiment_info = experiment_info,
                 agg_pars = agg_pars,
                 n_cells = clus_tab )
    # output
    if (message) {
        message("Output data ...") }
    sce <- SingleCellExperiment(assays = dat_list,
                                colData = sample_info,
                                metadata = meta)
    return(sce)
    
    
}
