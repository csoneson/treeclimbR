## If data is a matrix, the estimated parameters pi and theta are organized as a
## list and returned as output.
## If data is a list, we confirm that this list includes two elements pi and
## theta, and directly return data as output
#' @keywords internal
#' @noRd
#'
#' @importFrom dirmult dirmult
#'
.estimateA <- function(obj) {
    if (is.list(obj)) {
        return(.estimateB(obj))
    } else {
        if (is.null(rownames(obj))) {
            stop("obj must have rownames")
        }
        estP <- rep(0, nrow(obj))
        names(estP) <- rownames(obj)

        DirMultOutput <- dirmult::dirmult(data = t(obj))

        ## Tip proportion
        estP[names(DirMultOutput$pi)] <- DirMultOutput$pi

        ## alpha parameter for Dirichlet distribution
        theta <- DirMultOutput$theta
        parList <- list(pi = estP, theta = theta)
    }

    return(parList)
}

## If obj is a list, we confirm that this list includes two elements pi and
## theta, and directly return obj as output
#' @keywords internal
#' @noRd
#'
.estimateB <- function(obj) {
    .assertVector(x = obj, type = "list")

    ind <- setequal(names(obj), c("pi", "theta"))
    if (!ind) {
        stop("obj is a list; it should contain pi and theta")
    }
    return(obj)
}

## If obj is a TreeSummarizedExperiment, the estimated parameters pi and theta
## are organized as a list and stored in metadata with name assays.par
#' @keywords internal
#' @noRd
#'
#' @importFrom TreeSummarizedExperiment rowLinks
#' @importFrom S4Vectors metadata
#' @importFrom SummarizedExperiment assays
#' @importFrom dirmult dirmult
#'
.estimateC <- function(obj, assay) {
    .assertVector(x = obj, type = "TreeSummarizedExperiment")

    ## Node label
    nodeLab <- TreeSummarizedExperiment::rowLinks(obj)$nodeLab
    if (is.null(nodeLab)) {
        nodeLab <- rownames(obj)
    }

    ## Estimate parameters
    pars <- S4Vectors::metadata(obj)$assays.par
    ind <- setequal(names(pars), c("pi", "theta"))

    if (!ind) {
        ## Estimate parameters
        tdat <- t(SummarizedExperiment::assays(obj)[[assay]])
        colnames(tdat) <- nodeLab
        estP <- rep(0, ncol(tdat))
        names(estP) <- nodeLab
        DirMultOutput <- dirmult::dirmult(data = tdat)

        ## Tip proportion
        estP[names(DirMultOutput$pi)] <- DirMultOutput$pi

        ## alpha parameter for Dirichlet distribution
        theta <- DirMultOutput$theta
        S4Vectors::metadata(obj)$assays.par <- list(pi = estP, theta = theta)
    }

    return(obj)
}

#' Parameter estimation for Dirichlet-multinomial distribution
#'
#' \code{parEstimate} is a wrapper of the function
#' \code{\link[dirmult]{dirmult}} with default settings for \code{init},
#' \code{initscalar}, \code{epsilon}, \code{trace} and \code{mode}. It allows
#' the input \code{obj} to be either a \code{matrix} or a
#' \code{TreeSummarizedExperiment} and outputs the estimated values of
#' \code{pi} and \code{theta}.
#'
#' @author Ruizhu Huang, Charlotte Soneson
#' @export
#'
#' @param obj A matrix or \code{TreeSummarizedExperiment}, with samples in the
#'     columns and entities in the rows.
#' @param assay If \code{obj} is a \code{TreeSummarizedExperiment}, the name
#'     or index of the assay to use to estimate Dirichlet multinomial
#'     parameters. If \code{NULL}, the first assay will be used.
#'
#' @returns A list including the estimates of \code{pi} (a vector with one
#' element per row in \code{obj}) and \code{theta} (a scalar).
#'
#' @importFrom methods is
#' @importFrom SummarizedExperiment assayNames
#'
#' @examples
#' suppressPackageStartupMessages({
#'     library(TreeSummarizedExperiment)
#' })
#'
#' set.seed(1L)
#' y <- matrix(rnbinom(200, size = 1, mu = 10), nrow = 10)
#' colnames(y) <- paste("S", seq_len(20), sep = "")
#' rownames(y) <- tinyTree$tip.label
#' toy_tse <- TreeSummarizedExperiment(rowTree = tinyTree,
#'                                     assays = list(y))
#' res <- parEstimate(obj = toy_tse)
#' metadata(res)$assays.par
#'
parEstimate <- function(obj, assay = NULL) {

    ## Check arguments
    ## -------------------------------------------------------------------------
    stopifnot(methods::is(obj, "matrix") |
                  methods::is(obj, "list") |
                  methods::is(obj, "TreeSummarizedExperiment"))
    stopifnot(length(assay) == 1 || is.null(assay))
    if (methods::is(obj, "TreeSummarizedExperiment") && is.character(assay)) {
        stopifnot(assay %in% SummarizedExperiment::assayNames(obj))
    }
    if (methods::is(obj, "TreeSummarizedExperiment") && is.null(assay)) {
        assay <- 1
    }

    ## Estimate parameters
    ## -------------------------------------------------------------------------
    if (is.matrix(obj)) {
        out <- .estimateA(obj = obj)
    }

    if (is.list(obj)) {
        out <- .estimateB(obj = obj)
    }

    if (is(obj, "TreeSummarizedExperiment")) {
        out <- .estimateC(obj = obj, assay = assay)
    }

    return(out)
}
