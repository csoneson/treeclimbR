#' Wrapper applying an edgeR differential analysis workflow
#'
#' \code{edgerWrp} is a wrapper using functions from the \code{\link{edgeR}}
#' (Robinson et al. 2010, \emph{Bioinformatics}; McCarthy et al. 2012,
#' \emph{Nucleic Acids Research}) to fit models and perform a moderated test for
#' each entity.
#'
#' The function performs the following steps:
#' \itemize{
#' \item Create a \code{\link[edgeR]{DGEList}} object. If \code{lib_size} is
#' given, set the library sizes to these values, otherwise use the column sums
#' of the count matrix.
#' \item If \code{normalize} is \code{TRUE}, estimate normalization factors
#' using \code{\link[edgeR]{calcNormFactors}}.
#' \item Estimate dispersions with \code{\link[edgeR]{estimateDisp}}.
#' \item Depending on the value of \code{option}, apply either the LRT or
#' QLF edgeR workflows (i.e., either \code{\link[edgeR]{glmFit}} +
#' \code{\link[edgeR]{glmLRT}} or \code{\link[edgeR]{glmQLFit}} +
#' \code{\link[edgeR]{glmQLFTest}}), testing for the specified contrast.
#' }
#'
#' @author Ruizhu Huang
#' @export
#'
#' @param count A matrix with features (e.g., genes or microbes) in rows and
#'     samples in columns.
#' @param lib_size A numeric vector with library sizes for each sample. If
#'     \code{NULL} (default), the column sums of \code{count} are used.
#' @param option Either \code{"glm"} or \code{"glmQL"}. If \code{"glm"},
#'     \code{\link[edgeR]{glmFit}} and \code{\link[edgeR]{glmLRT}} are used;
#'     otherwise, \code{\link[edgeR]{glmQLFit}} and
#'     \code{\link[edgeR]{glmQLFTest}} are used. Details about the difference
#'     between the two options can be found in the help pages of
#'     \code{\link[edgeR]{glmQLFit}}.
#' @param design A numeric design matrix, e.g. created by
#'     \code{\link[stats]{model.matrix}}. Please refer to \code{design} in
#'     \code{\link[edgeR]{glmQLFit}} and \code{\link[edgeR]{glmFit}} for
#'     more details.
#' @param contrast A numeric vector specifying one contrast of
#'     the linear model coefficients to be tested. Its length
#'     must equal the number of columns of \code{design}. If \code{NULL}, the
#'     last coefficient will be tested. Please refer to \code{contrast}
#'     in \code{\link[edgeR]{glmQLFTest}} and \code{\link[edgeR]{glmLRT}} for
#'     more details.
#' @param normalize A logical scalar, specifying whether normalization factors
#'     should be calculated (using \code{\link[edgeR]{calcNormFactors}}).
#' @param normalize_method Normalization method to be used. Please refer to
#'     \code{method} in \code{\link[edgeR]{calcNormFactors}} for more details.
#' @param ... More arguments to pass to \code{\link[edgeR]{glmFit}}
#'     (if \code{option = "glm"} or \code{\link[edgeR]{glmQLFit}}
#'     (if \code{option = "glmQL"}).
#'
#' @importFrom edgeR DGEList calcNormFactors estimateDisp glmFit glmLRT
#'     glmQLFit glmQLFTest
#'
#' @return The output of \code{\link[edgeR]{glmQLFTest}} or
#'     \code{\link[edgeR]{glmLRT}} depending on the specified \code{option}.
#'
#' @examples
#' ## Read example data
#' x <- readRDS(system.file("extdata/da_sim_100_30_18de.rds",
#'                          package = "treeclimbR"))
#'
#' ## Run differential abundance analysis
#' out <- edgerWrp(count = assay(x), option = "glm",
#'                 design = model.matrix(~ group, data = colData(x)),
#'                 contrast = c(0, 1))
#'
edgerWrp <- function(count, lib_size = NULL, option = c("glm", "glmQL"),
                     design, contrast = NULL, normalize = TRUE,
                     normalize_method = "TMM", ...) {

    ## Check input arguments
    ## -------------------------------------------------------------------------
    .assertVector(x = count, type = "matrix")
    .assertVector(x = lib_size, type = "numeric", allowNULL = TRUE,
                  len = ncol(count))
    .assertVector(x = design, type = "matrix")
    .assertVector(x = contrast, type = "numeric", allowNULL = TRUE,
                  len = ncol(design))
    .assertScalar(x = normalize, type = "logical")
    .assertScalar(x = normalize_method, type = "character")

    ## If the library size isn't given, the column sum is used
    ## -------------------------------------------------------------------------
    if (is.null(lib_size)) {
        lib_size <- apply(count, 2, sum)
    }

    ## Create the DGEList
    ## -------------------------------------------------------------------------
    y <- edgeR::DGEList(count, remove.zeros = FALSE)
    y$samples$lib.size <- lib_size

    ## Normalization
    ## -------------------------------------------------------------------------
    if (normalize) {
        y <- edgeR::calcNormFactors(y, method = normalize_method)
    }

    ## Estimate dispersion
    ## -------------------------------------------------------------------------
    y <- edgeR::estimateDisp(y, design = design)

    ## Fit model and test contrast
    ## -------------------------------------------------------------------------
    option <- match.arg(option)
    if (option == "glm") {
        fit <- edgeR::glmFit(y, design = design, ...)
        lrt <- edgeR::glmLRT(fit, contrast = contrast)
    } else {
        fit <- edgeR::glmQLFit(y, design = design, ...)
        lrt <- edgeR::glmQLFTest(fit, contrast = contrast)
    }

    return(lrt)
}
