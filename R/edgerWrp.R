#' A wrapper function of edgeR
#'
#' \code{edgerWrp} is a wrapper using functions from the \code{\link{edgeR}}
#' (Robinson et al. 2010, \emph{Bioinformatics}; McCarthy et al. 2012,
#' \emph{Nucleic Acids Research}) to fit models and calculate moderated test for
#' each entity. We have used \code{\link[edgeR]{estimateDisp}} to estimate the
#' dispersion. The statistical methods implemented in the \code{edgeR} package
#' were originally designed for the analysis of gene expression data such as
#' RNA-sequencing counts. Here, we apply these methods to counts that might be
#' the abundance of  microbes or cells.
#'
#' The experimental design must be specified using a design matrix. The
#' customized design matrix could be given by \code{design}.
#'
#' Normalization for samples is automatically performed by \code{edgeR} package.
#' More details about the calculation of normalization factor could be found
#' from \code{\link[edgeR]{calcNormFactors}}. 
#'
#' @param count A matrix with features (e.g., genes or microbes) on rows and
#'   samples on columns.
#' @param lib_size A numeric vector to specify the library size. If NULL
#'   (default), the column sum of \code{count} is used.
#' @param option  \code{"glm"} or \code{"glmQL"}. If \code{"glm"},
#'   \code{\link[edgeR]{glmFit}} and \code{\link[edgeR]{glmLRT}} are used;
#'   otherwise, \code{\link[edgeR]{glmQLFit}} and
#'   \code{\link[edgeR]{glmQLFTest}} are used. Details about the difference
#'   between two options are in the help page of \code{\link[edgeR]{glmQLFit}}.
#' @param design A numeric matrix. It must be of full column rank created by
#'   \code{\link[stats]{model.matrix}}. Please refer to \code{design} in
#'   \code{\link[edgeR]{glmQLFit}} and \code{\link[edgeR]{glmFit}}
#' @param contrast A numeric vector specifying one contrast of
#'   the linear model coefficients to be tested equal to zero. Its length
#'   must equal to the number of columns of design. If NULL, the last
#'   coefficient will be tested equal to zero. Please refer to \code{contrast} in
#'   \code{\link[edgeR]{glmQLFTest}} and \code{\link[edgeR]{glmLRT}}.
#' @param normalize A logical value, TRUE or FALSE, to specify whether
#'   \code{count} should be normalized. The default is TRUE.
#' @param normalize_method Normalization method to be used. Please refer to
#'   \code{method} in \code{\link[edgeR]{calcNormFactors}} for more details.
#' @param ... More arguments to pass to \code{\link[edgeR]{glmFit}}
#'   (\code{option = "glm"} or \code{\link[edgeR]{glmQLFit}} (\code{option =
#'   "glmQL"}).
#' @importFrom edgeR DGEList calcNormFactors estimateDisp glmFit glmLRT
#'   glmQLFit glmQLFTest
#' @export 
#' @return The output of \code{\link[edgeR]{glmQLFTest}} or
#'   \code{\link[edgeR]{glmLRT}} depends on the specified \code{option}.
#' @examples 
#' 
#' nlibs <- 6
#' ngenes <- 100
#' dispersion.true <- 0.1
#' 
#' # Make first 5 gene respond to covariate x
#' x <- rep(1:2, each = 3)
#' design <- model.matrix(~ x)
#' beta.true <- cbind(Beta1 = 2, Beta2 = c(rep(2, 5), rep(0, ngenes-5)))
#' mu.true <- 2^(beta.true %*% t(design))
#' 
#' # Generate count data
#' y <- rnbinom(ngenes*nlibs, mu = mu.true, size = 1/dispersion.true)
#' y <- matrix(y, ngenes, nlibs)
#' colnames(y) <- paste(rep(LETTERS[1:2], each = 3),
#'                      rep(1:3, 2), sep = "_")
#' rownames(y) <- paste("gene", 1:ngenes, sep=".")
#' 
#' # analysis
#' out <- edgerWrp(count = y, option = "glm", 
#'                 design = design, normalize = FALSE)
#' 

edgerWrp <- function(count, lib_size = NULL , 
                     option = c("glm", "glmQL"),
                     design = NULL, contrast = NULL,
                     normalize = TRUE, 
                     normalize_method = "TMM", ...) {
    
    # if the library size isn't given, the column sum is used
    if (is.null(lib_size)) {
        lib_size <- apply(count, 2, sum)
    }
    
    # create the DGEList
    y <- DGEList(count, remove.zeros = FALSE)
    y$samples$lib.size <- lib_size
    
    # normalisation
    if (normalize) {
        y <- calcNormFactors(y, method = normalize_method)
    }
    
    # estimate dispersion
    y <- estimateDisp(y, design = design)
    
    # model and test contrast
    option <- match.arg(option)
    if (option == "glm") {
        fit <- glmFit(y, design = design)
        lrt <- glmLRT(fit, contrast = contrast)
    } else {
        fit <- glmQLFit(y, design = design)
        lrt <- glmQLFTest(fit, contrast = contrast)
    }
    
    return(lrt)
}