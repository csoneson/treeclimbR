#' Test for differential state using edgeR
#'
#' \code{runDA} tests differential state of entities using functions from the
#' \code{\link{edgeR}}. This adapts \code{\link{edgerWrp}} to accept input as a
#' \code{\link[SummarizedExperiment]{SummarizedExperiment}} (\strong{SE}) object
#' instead of \code{matrix}. To clarify, the count table is stored in the
#' \code{assays} table of \strong{SE} object. Each matrix in \code{assays}
#' stores data corresponding to a node of the tree. Samples are in columns and
#' features, e.g., genes or antibodies, are in rows. The sample information is
#' in \code{colData}. The tree that stores the hierarchical information about
#' the \code{assays} is provided via the argument \code{tree}. \code{assays} are
#' named with the alias labels of nodes that could be mapped to the provided
#' tree. The experimental design is specified by a design matrix and provide to
#' model via the argument \code{design}. Alternatively, users could specify the
#' name of columns (from \code{colData}) that are passed to
#' \code{\link[stats]{model.matrix}} to generate the design matrix.
#' 
#' @param SE A SummarizedExperiment object.
#' @param tree A phylo object. Each matrix in \strong{assays} (of \code{SE})
#'   stores data on one node of the tree.
#' @param option  \code{"glm"} or \code{"glmQL"}. If \code{"glm"},
#'   \code{\link[edgeR]{glmFit}} and \code{\link[edgeR]{glmLRT}} are used;
#'   otherwise, \code{\link[edgeR]{glmQLFit}} and
#'   \code{\link[edgeR]{glmQLFTest}} are used. Details about the difference
#'   between two options are in the help page of \code{\link[edgeR]{glmQLFit}}.  
#' @param design A numeric matrix. It must be of full column rank. Defaults to
#'   use all columns of sample annotation data to create the design matrix. The
#'   sample annotation data is stored in the \code{colData} of \code{tse} when
#'   \code{onRow = TRUE}; otherwise it is in the \code{rowData}. Note: Users
#'   should check whether the default created design matrix is exactly what they
#'   want or create their own design matrix using
#'   \code{\link[stats]{model.matrix}}.
#' @param contrast A numeric vector specifying one contrast of
#'   the linear model coefficients to be tested equal to zero. Its length
#'   must equal to the number of columns of design. If NULL, the last
#'   coefficient will be tested equal to zero.
#' @param filter_min_count A numeric value. It's passed to \strong{min.count} of
#'   \code{\link[edgeR]{filterByExpr}}.
#' @param filter_min_total_count A numeric value. It's passed to
#'   \strong{min.total.count} of \code{\link[edgeR]{filterByExpr}}.
#' @param filter_large_n A numeric value. It's passed to \strong{large.n} of
#'   \code{\link[edgeR]{filterByExpr}}.
#' @param filter_min_prop A numeric value. It's passed to \strong{min.prop} of
#'   \code{\link[edgeR]{filterByExpr}}.
#' @param normalize A logical value, TRUE or FALSE. The default is TRUE.
#' @param normalize_method Normalization method to be used. See
#'   \code{\link[edgeR]{calcNormFactors}} for more details.
#'   It could be "bonferroni", "holm", "hochberg", "hommel", "BH", or "BY". This
#'   is passed to \code{adjust.method} of \code{\link[edgeR]{topTags}}
#' @param group_column The column name of group 
#' @param design_terms The names of columns from \code{colData} (if samples in
#'   columns) that are used to generate design matrix. This is ignored if
#'   \strong{design} is provided.
#' @param message A logical value, TRUE or FALSE. If TRUE, the running process
#'   is printed out.
#' @param ... More arguments to pass to \code{\link[edgeR]{glmFit}}
#'   (\code{option = "glm"} or \code{\link[edgeR]{glmQLFit}} (\code{option =
#'   "glmQL"}).
#' @importFrom edgeR DGEList calcNormFactors estimateDisp glmFit glmLRT
#'   glmQLFit glmQLFTest cpm filterByExpr
#' @importFrom methods is
#' @importFrom SummarizedExperiment assays colData rowData
#' @importFrom stats model.matrix
#' @export
#' @author Ruizhu Huang
#' @return a list includes \strong{edgeR_results}, \strong{tree}, and
#'   \strong{nodes_drop}.
#' \describe{
#'   \item{edgeR_results} {A list, each of its element is the output of
#'   \code{\link[edgeR]{glmQLFTest}} or \code{\link[edgeR]{glmLRT}} depends on
#'   the specified \code{option}. The list is named by the alias label of nodes}
#'   \item{tree} {The hiearchical structure of entities that was stored in the
#'      input \code{tse}}
#'   \item{nodes_drop} {A vector storing the alias node labels of entities. 
#'      These nodes have no sufficient data to run analysis. }
#' }
#' 
#' @examples
#' 
#' library(TreeSummarizedExperiment)
#' 
#' data(tinyTree)
#' 
#' set.seed(1)
#' # data: genes in rows, cells in columns
#' ncell <- 500
#' ngene <- 20
#' dat <- matrix(rpois(ngene*ncell, lambda = 20), nrow = ngene)
#' 
#' # cells belong to different clusters
#' clus <- sample(tinyTree$tip.label, 500, replace = TRUE)
#' sid <- sample(1:4, 500, replace = TRUE)
#' gid <- ifelse(sid %in% 1:2, "A", "B")
#' colD <- data.frame(cluster_id = clus,
#'                    sample_id = sid,
#'                    group = gid,
#'                    stringsAsFactors = FALSE)
#' # TreeSummarizedExperiment object
#' d_tse <- TreeSummarizedExperiment(
#'             assays = list(dat),
#'             colData = colD,
#'             colTree = tinyTree, 
#'             colNodeLab = clus)
#' d_se <- aggDS(TSE = d_tse, 
#'               assay = 1, sample_id = "sample_id",
#'               group_id = "group",
#'               cluster_id = "cluster_id", FUN = sum)
#' 
#' res <- runDS(SE = d_se, tree = tinyTree, 
#'              option = "glm", group_column = "group",
#'              design_terms = "group")
#' topNodes(res, type = "DS")

runDS <- function(SE, tree, 
                  option = c("glm", "glmQL"),
                  design = NULL, contrast = NULL, 
                  filter_min_count = 10, 
                  filter_min_total_count = 15,
                  filter_large_n = 10,
                  filter_min_prop = 0.7, 
                  normalize = TRUE, 
                  normalize_method = "TMM",
                  group_column = "group_id", 
                  design_terms = "group_id",
                  message = TRUE, ...) {
    
    # the alias of nodes
    alias <- names(assays(SE))
    res <- vector("list", length(alias))
    names(res) <- alias
    nodes_list <- res
    
    for (i in seq_along(alias)) {
        #message(i)
        if (message) {
            message(i, " out of ", length(alias),
                    " nodes finished", "\r", appendLF = FALSE)
            flush.console()}
        
        # extract: count, sample information, node information
        count <- assays(SE)[[i]]
        sp_info <- colData(SE)
        
        
        # ---------------------------- sufficient data? -----------------------
        # check whether there are samples from different groups to run analysis
        sp_keep <- colSums(count, na.rm = TRUE) > 0
        if (sum(!sp_keep)) {
            sp_info <- sp_info[sp_keep, ,drop = FALSE]
         } 
        sp_gr <- unique(sp_info[, group_column])
        
        ## if not, return NULL
        if (length(sp_gr) == 1) {
            res[[i]] <- NULL 
            nodes_list[[i]] <- alias[i]
        } else {
        ## if yes, run analysis
            res[[i]] <- .DS(SE = SE,
                            assay = alias[i],
                            option = option,
                            design = design, contrast = contrast, 
                            filter_min_count = filter_min_count, 
                            filter_min_total_count = 
                                filter_min_total_count,
                            filter_large_n = filter_large_n,
                            filter_min_prop = filter_min_prop, 
                            normalize = normalize, 
                            normalize_method = normalize_method,
                            group_column = group_column, 
                            design_terms = design_terms)
            
            
            
            
        }
        out <- list(edgeR_results = res,
                    tree = tree,
                    nodes_drop = unlist(nodes_list))  
    }
    
   return(out)
}

.DS <- function(SE, feature_on_row = TRUE, assay,
                option = c("glm", "glmQL"),
                design = NULL, contrast = NULL, 
                filter_min_count = 10, 
                filter_min_total_count = 15,
                filter_large_n = 10,
                filter_min_prop = 0.7, 
                normalize = TRUE, normalize_method = "TMM",
                group_column = "group", 
                design_terms = "group", ...) {
    
    # require a SummarizedExperiment object
    if (!is(SE, "SummarizedExperiment")) {
        stop("SE should be a SummarizedExperiment.")
    }
    
    # decide the model to use
    option <- match.arg(option)
    
    # extract: count, sample information, node information
    count <- assays(SE)[[assay]]
    sp_info <- colData(SE)
    sp_min_0 <- min(table(sp_info[[group_column]]))
    
    sp_keep <- colSums(count, na.rm = TRUE) > 0
    if (sum(!sp_keep)) {
        count <- count[, sp_keep]
        sp_info <- sp_info[sp_keep, ,drop = FALSE]
    } 
    sp_min_1 <- min(table(sp_info[[group_column]]))
    # ---------------------------- design matrix -----------------------
    if (is.null(design)) {
        formula <- as.formula(paste("~", paste(design_terms, collapse = "+")))
        design <- model.matrix(formula, data = data.frame(sp_info))
    }
    
    ### ------------ filter lowly expressed features --------------------------
    # This is from (https://f1000research.com/articles/5-1438/v2).
    # The filtering is on count-per-million (CPM)
    ### -----------------------------------------------------------------------
    # Samples with zero library size are removed manually above, so we need to
    # recalculate the min.prop here
    filter_min_prop <- filter_min_prop * sp_min_0/sp_min_1
    keep <- filterByExpr(count, design = design, 
                         min.count = filter_min_count, 
                         min.total.count = filter_min_total_count,
                         large.n = filter_large_n,
                         min.prop = filter_min_prop)
    
    count_keep <- count[keep, ]
    isLow <- !keep
    feature_drop <- rownames(count)[isLow]
    # ---------------------------- run edgeR -----------------------
    lrt <- edgerWrp(count = count_keep, lib_size = NULL , 
                    option = option,
                    design = design, contrast = contrast,
                    normalize = normalize, 
                    normalize_method = normalize_method, ...)
    
    # -------------------------- output  ---------------------------------
    return(lrt)
}


        