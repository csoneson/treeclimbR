#' Test for differential abundance using edgeR
#'
#' Test for differential abundance of entities using functions from
#' the \code{\link{edgeR}} package. This adapts \code{\link{edgerWrp}} to
#' accept input as a
#' \code{\link[TreeSummarizedExperiment]{TreeSummarizedExperiment}}
#' (TSE) object instead of a \code{matrix}. Features could be
#' represented in either rows or columns. By default, features are in the rows.
#' Then, samples are in columns and the sample information is in \code{colData}.
#' The tree that stores the hierarchical information about features is in
#' \code{rowTree}. Each row of the \code{assays} can be mapped to a node of
#' the tree. Data on rows that are mapped to internal nodes is generated from
#' data on leaf nodes. Normalization for samples is automatically performed by
#' \code{edgeR} and the library size is calculated using features that
#' are mapped to leaf nodes.
#'
#' The experimental design is specified by a design matrix and provided
#' via the argument \code{design}. More details about the calculation of
#' normalization factor could be found from
#' \code{\link[edgeR]{calcNormFactors}}.
#'
#' @author Ruizhu Huang
#' @export
#'
#' @param TSE A \code{TreeSummarizedExperiment} object.
#' @param feature_on_row A logical scalar. If \code{TRUE} (default),
#'     features or entities (e.g. genes, OTUs) are in rows of the \code{assays}
#'     tables, and samples are in columns; otherwise, it's the other way around.
#' @param assay A numeric index or assay name to specify which assay from
#'     \code{assays} is used for analysis.
#' @param option Either \code{"glm"} or \code{"glmQL"}. If \code{"glm"},
#'     \code{\link[edgeR]{glmFit}} and \code{\link[edgeR]{glmLRT}} are used;
#'     otherwise, \code{\link[edgeR]{glmQLFit}} and
#'     \code{\link[edgeR]{glmQLFTest}} are used. Details about the difference
#'     between two options are in the help page of
#'     \code{\link[edgeR]{glmQLFit}}.
#' @param design A numeric design matrix. If \code{NULL}, all columns of
#'     the sample annotation will be used to create the design matrix.
#' @param contrast A numeric vector specifying one contrast of
#'     the linear model coefficients to be tested equal to zero. Its length
#'     must equal to the number of columns of design. If \code{NULL}, the last
#'     coefficient will be tested equal to zero.
#' @param filter_min_count A numeric value, passed to \strong{min.count} of
#'     \code{\link[edgeR]{filterByExpr}}.
#' @param filter_min_total_count A numeric value, passed to
#'     \strong{min.total.count} of \code{\link[edgeR]{filterByExpr}}.
#' @param filter_large_n A numeric value, passed to \strong{large.n} of
#'     \code{\link[edgeR]{filterByExpr}}.
#' @param filter_min_prop A numeric value, passed to \strong{min.prop} of
#'     \code{\link[edgeR]{filterByExpr}}.
#' @param normalize A logical scalar indicating whether to estimate
#'     normalization factors (using \code{\link[edgeR]{calcNormFactors}}).
#' @param normalize_method Normalization method to be used. See
#'     \code{\link[edgeR]{calcNormFactors}} for more details.
#' @param group_column The name of the column in the sample annotation
#'     providing group labels for samples (currently not used).
#' @param design_terms The names of columns from the sample annotation
#'     that will be used to generate the design matrix. This is ignored if
#'     \strong{design} is provided.
#' @param ... More arguments to pass to \code{\link[edgeR]{glmFit}}
#'     (\code{option = "glm"} or \code{\link[edgeR]{glmQLFit}}
#'     (\code{option = "glmQL"}).
#'
#' @returns A list with entries \strong{edgeR_results}, \strong{tree}, and
#' \strong{nodes_drop}.
#' \describe{
#'     \item{edgeR_results}{The output of \code{\link[edgeR]{glmQLFTest}} or
#'          \code{\link[edgeR]{glmLRT}} depending on the specified
#'          \code{option}.}
#'     \item{tree}{The hierarchical structure of entities that was stored in the
#'          input \code{TSE}.}
#'     \item{nodes_drop}{A vector storing the alias node labels of entities
#'          that are filtered before analysis due to low counts.}
#' }
#'
#' @importFrom TreeSummarizedExperiment TreeSummarizedExperiment rowLinks
#'     colLinks colTree rowTree
#' @importFrom edgeR filterByExpr
#' @importFrom SummarizedExperiment assays colData rowData assayNames
#'
#' @examples
#' suppressPackageStartupMessages({
#'     library(TreeSummarizedExperiment)
#' })
#'
#' ## Load example data set
#' lse <- readRDS(system.file("extdata", "da_sim_100_30_18de.rds",
#'                            package = "treeclimbR"))
#'
#' ## Aggregate counts on internal nodes
#' nodes <- showNode(tree = tinyTree, only.leaf = FALSE)
#' tse <- aggTSE(x = lse, rowLevel = nodes)
#'
#' dd <- model.matrix(~ group, data = colData(tse))
#' out <- runDA(TSE = tse, feature_on_row = TRUE,
#'              assay = 1, option = "glmQL",
#'              design = dd, contrast = NULL,
#'              normalize = TRUE, filter_min_count = 2)
#' names(out)
#' out$nodes_drop
#' edgeR::topTags(out$edgeR_results, sort.by = "PValue")
#'
runDA <- function(TSE, feature_on_row = TRUE, assay = NULL,
                  option = c("glm", "glmQL"),
                  design = NULL, contrast = NULL,
                  filter_min_count = 10,
                  filter_min_total_count = 15,
                  filter_large_n = 10,
                  filter_min_prop = 0.7,
                  normalize = TRUE, normalize_method = "TMM",
                  group_column = "group",
                  design_terms = "group", ...) {

    ## Check arguments
    ## -------------------------------------------------------------------------
    .assertVector(x = TSE, type = "TreeSummarizedExperiment")
    .assertScalar(x = feature_on_row, type = "logical")
    .assertVector(x = design, type = "matrix", allowNULL = TRUE)
    .assertVector(x = contrast, type = "numeric", allowNULL = TRUE)
    .assertScalar(x = filter_min_count, type = "numeric")
    .assertScalar(x = filter_min_total_count, type = "numeric")
    .assertScalar(x = filter_large_n, type = "numeric")
    .assertScalar(x = filter_min_prop, type = "numeric")
    .assertScalar(x = normalize, type = "logical")
    .assertScalar(x = normalize_method, type = "character")
    .assertVector(x = design_terms, type = "character", allowNULL = TRUE)

    ## If not specified, the first assays is used
    if (is.null(assay)) {
        assay <- 1
    } else {
        if (is.character(assay)) {
            stopifnot(assay %in% SummarizedExperiment::assayNames(TSE))
        } else if (is.numeric(assay)) {
            stopifnot(assay <= length(SummarizedExperiment::assays(TSE)))
        } else {
            stop("'assay' must be of class 'numeric' or 'character'")
        }
    }

    ## Get the model to use
    option <- match.arg(option)

    ## Extract required information
    ## -------------------------------------------------------------------------
    if (feature_on_row) {
        ## counts
        count <- SummarizedExperiment::assays(TSE)[[assay]]
        ## node information
        ld <- TreeSummarizedExperiment::rowLinks(TSE)
        ## sample information
        sp_info <- SummarizedExperiment::colData(TSE)
        ## tree
        tree <- TreeSummarizedExperiment::rowTree(TSE)
    } else {
        ## counts
        count <- t(SummarizedExperiment::assays(TSE)[[assay]])
        ## node information
        ld <- TreeSummarizedExperiment::colLinks(TSE)
        ## sample information
        sp_info <- SummarizedExperiment::rowData(TSE)
        ## tree
        tree <- TreeSummarizedExperiment::colTree(TSE)
    }

    ## Add rownames
    rownames(ld) <- rownames(count) <- ld$nodeLab_alias

    ## Design matrix
    ## -------------------------------------------------------------------------
    if (is.null(design)) {
        formula <- as.formula(paste("~", paste(design_terms, collapse = "+")))
        design <- model.matrix(formula, data = data.frame(sp_info))
    }

    ## Filter lowly abundant features
    ## -------------------------------------------------------------------------
    ## Use data on leaf level to get library size
    tip_count <- count[ld$isLeaf, ]
    lib_size <- apply(tip_count, 2, sum)
    keep <- edgeR::filterByExpr(count, design = design,
                                lib.size = lib_size,
                                min.count = filter_min_count,
                                min.total.count = filter_min_total_count,
                                large.n = filter_large_n,
                                min.prop = filter_min_prop)

    count_keep <- count[keep, , drop = FALSE]
    isLow <- !keep
    feature_drop <- rownames(count)[isLow]

    ## Run edgeR
    ## -------------------------------------------------------------------------
    lrt <- edgerWrp(count = count_keep, lib_size = lib_size ,
                    option = option,
                    design = design, contrast = contrast,
                    normalize = normalize,
                    normalize_method = normalize_method, ...)

    ## Assemble output
    ## -------------------------------------------------------------------------
    if (sum(isLow)) {
        out <- list(edgeR_results = lrt,
                    nodes_drop = feature_drop,
                    tree = tree)
    } else {
        out <- list(edgeR_results = lrt,
                    nodes_drop = NULL,
                    tree = tree)
    }

   return(out)
}



