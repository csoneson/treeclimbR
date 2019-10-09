#' Test for differential abundance using edgeR
#'
#' \code{runDA} tests differential abundance of entities using functions from
#' the \code{\link{edgeR}}. This adapts \code{\link{edgerWrp}} to accept input
#' as a \code{\link[TreeSummarizedExperiment]{TreeSummarizedExperiment}}
#' (\strong{TSE}) object instead of \code{matrix}. To clarify, the count table
#' is stored in the \code{assays} table of \strong{TSE} object. Features could
#' be stored either in rows or columns. Here, we consider features are in rows.
#' Then, samples are in columns and the sample information is in \code{colData}.
#' The tree that stores the hierarchical information about features is in
#' \code{rowTree}. Each row of \code{assays} table could be mapped to a node of
#' the tree. Data on rows that are mapped to internal nodes is generated from
#' data on leaf nodes. Normalization for samples is automatically performed by
#' \code{edgeR} package and the library size is calculated using features that
#' are mapped to leaf nodes. More details are in Vignette.
#'
#' The experimental design is specified by a design matrix and provide to model
#' via the argument \code{design}. More details about the calculation of
#' normalization factor could be found from
#' \code{\link[edgeR]{calcNormFactors}}.
#' 
#' @param tse A TreeSummarizedExperiment object.
#' @param feature_on_row A logical value, TRUE or FALSE. If TRUE (default),
#'   features or entities (e.g. genes, OTUs) are in rows of the \code{assays}
#'   tables, and samples are in columns; otherwise, the other way around.
#' @param assay A numeric value to specify which table from \code{assays} is
#'   used for analysis.
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
#' @param design_term The names of columns from \code{colData} (if samples in
#'   columns) that are used to generate design matrix. This is ignored if
#'   \strong{design} is provided.
#' @param ... More arguments to pass to \code{\link[edgeR]{glmFit}}
#'   (\code{option = "glm"} or \code{\link[edgeR]{glmQLFit}} (\code{option =
#'   "glmQL"}).
#' @import TreeSummarizedExperiment
#' @importFrom edgeR DGEList calcNormFactors estimateDisp glmFit glmLRT
#'   glmQLFit glmQLFTest cpm filterByExpr
#' @importFrom methods is
#' @importFrom SummarizedExperiment assays colData rowData
#' @export
#' @author Ruizhu Huang
#' @return a list includes \strong{edgeR_results}, \strong{tree}, and
#'   \strong{nodes_drop}.
#' \describe{
#'   \item{edgeR_results} {The output of \code{\link[edgeR]{glmQLFTest}} or
#'      \code{\link[edgeR]{glmLRT}} depends on the specified \code{option}.}
#'   \item{tree} {The hiearchical structure of entities that was stored in the
#'      input \code{tse}}
#'   \item{nodes_drop} {A vector storing the alias node labels of entities. 
#'      These entities are filtered before analysis due to low counts. }
#' }
#' 
#' @examples
#' library(TreeSummarizedExperiment)
#' library(treeAGG2)
#' set.seed(1)
#' count <- matrix(rnbinom(300,size=1,mu=10),nrow=10)
#' colnames(count) <- paste(rep(LETTERS[1:3], each = 10), rep(1:10,3), sep = "_")
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
#' nodes <- showNode(tree = tinyTree, only.leaf = FALSE)
#' tse <- aggValue(x = lse, rowLevel = nodes)
#' 
#' dd <- model.matrix( ~ group, data = colInf)
#' out <- runDA(tse = tse, feature_on_row = TRUE,
#'              assay = 1, option = "glmQL",
#'              design = dd, contrast = NULL, 
#'              normalize = TRUE, 
#'              group_column = "group")

runDA <- function(tse, feature_on_row = TRUE, assay = NULL,
                  option = c("glm", "glmQL"),
                  design = NULL, contrast = NULL, 
                  filter_min_count = 10, 
                  filter_min_total_count = 15,
                  filter_large_n = 10,
                  filter_min_prop = 0.7, 
                  normalize = TRUE, normalize_method = "TMM",
                  adjust_method = "BH", 
                  group_column = "group", 
                  design_terms = "group", ...) {
    
    # require a TreeSummarizedExperiment object
    if (!is(tse, "TreeSummarizedExperiment")) {
        stop("tse should be a TreeSummarizedExperiment.")
    }
    
    # If not specified, the first assays is used
    if (is.null(assay)) {
        assay <- 1
    }
    
    # decide the model to use
    option <- match.arg(option)
    
    # extract: count, sample information, node information
    if (feature_on_row) {
        # count
        count <- assays(tse)[[assay]]
        
        # node information
        ld <- rowLinks(tse)
        
        # sample information
        sp_info <- colData(tse)
        
        # tree
        tree <- rowTree(tse)
    } else {
        # the count table
        count <- t(assays(tse)[[assayNum]])
        
        # extract link data
        ld <- colLinks(tse)
        
        # sample information
        sp_info <- rowData(tse)
        
        # tree
        tree <- colTree(tse)
    }
    
    # add rownames
    rownames(ld) <- rownames(count) <- ld$nodeLab_alias
    
    
    # ---------------------------- design matrix -----------------------
    if (is.null(design)) {
        formula <- as.formula(paste("~", paste(design_terms, collapse = "+")))
        design <- model.matrix(formula, data = data.frame(sp_info))
    }
    
    ### ------------ filter lowly expressed features --------------------------
    # This is from (https://f1000research.com/articles/5-1438/v2).
    # The filtering is on count-per-million (CPM)
    ### -----------------------------------------------------------------------
    # the library size (use data on leaf level)
    # the cutoff
    tip_count <- count[ld$isLeaf, ]
    lib_size <- apply(tip_count, 2, sum)
    keep <- filterByExpr(count, design = design, 
                         lib.size = lib_size,
                         min.count = filter_min_count, 
                         min.total.count = filter_min_total_count,
                         large.n = filter_large_n,
                         min.prop = filter_min_prop)
     
    count_keep <- count[keep, ]
    isLow <- !keep
    feature_drop <- rownames(count)[isLow]
    # ---------------------------- run edgeR -----------------------
    lrt <- edgerWrp(count = count_keep, lib_size = lib_size , 
                    option = option,
                    design = design, contrast = contrast,
                    normalize = normalize, 
                    normalize_method = normalize_method, ...)
    
    # -------------------------- output  ---------------------------------
    if (sum(isLow)) {
        out <- list(edgeR_results = lrt,
                    tree = tree,
                    nodes_drop = feature_drop)
    } else {
        out <- list(edgeR_results = lrt,
                    tree = tree,
                    nodes_drop = NULL)
    }
    
   return(out)
    
}



