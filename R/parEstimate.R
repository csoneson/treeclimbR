# If data is a matrix, the estimated parameters pi and theta is organized as a
# list and returned as output.
# If data is a list, we confirme this list includes two elements pi and theta,
# and directly return data as output
.estimateA <- function(obj) {

    if (is.list(obj)) {
        ind <- setequal(names(obj), c("pi", "theta"))
        if (!ind) {
            stop("obj is a list; it should provide pi and theta")
        }
        parList <- obj
        } else {
            if (is.null(rownames(obj))) {
                stop("The rownames of obj are required.")
            }
            estP <- rep(0, nrow(obj))
            names(estP) <- rownames(obj)

            DirMultOutput <- dirmult(data = t(obj))
            # tip proportion
            estP[names(DirMultOutput$pi)] <- DirMultOutput$pi

            # parameter alpha for dirichlet distribution
            theta <- DirMultOutput$theta
            parList <-  list(pi = estP, theta = theta)
        }

    return(parList)
}

# If obj is a list, we confirme this list includes two elements pi and theta,
# and directly return obj as output
.estimateB <- function(obj) {
    ind <- setequal(names(obj), c("pi", "theta"))
    if (!ind) {
        stop("obj is a list; it should provide pi and theta")
    }
    parList <- obj


    return(parList)
    }

# If obj is a TreeSummarizedExperiment, the estimated parameters pi and theta
# is organized as a list and store in metadata with name assays.par
.estimateC <- function(obj) {

    # node label
    nodeLab <- rowLinks(obj)$nodeLab
    if (is.null(nodeLab)) {
        nodeLab <- rownames(obj)
    }

    # estimate parameters
    pars <- metadata(obj)$assays.par
    ind <- setequal(names(pars), c("pi", "theta"))

    if (!ind) {
        tdat <- t(assays(obj)[[1]])
        colnames(tdat) <- nodeLab
        estP <- rep(0, ncol(tdat))
        names(estP) <- nodeLab
        DirMultOutput <- dirmult(data = tdat)
        # tip proportion
        estP[names(DirMultOutput$pi)] <- DirMultOutput$pi


        # parameter alpha for dirichlet distribution
        theta <- DirMultOutput$theta
        metadata(obj)$assays.par <-  list(pi = estP, theta = theta)
    }

    return(obj)

}

#' Parameter estimation in Dirichlet-multinomial distribution
#'
#' \code{parEstimate} is a wrapper function generated from the function
#' \code{\link[dirmult]{dirmult}} with default setting on \code{init},
#' \code{initscalar}, \code{epsilon}, \code{trace} and \code{mode}. It allows
#' the input \code{obj} to accept \code{matrix} or
#' \code{leafSummarizedExperiment} and output the value \code{pi} and
#' \code{theta}.
#'
#' @param obj A matrix or leafSummarizedExperiment. Samples in the columns and
#'   entities in the rows.
#'
#' @importFrom dirmult dirmult
#' @importFrom methods is
#' @import TreeSummarizedExperiment
#' @importFrom SummarizedExperiment assays
#' @importFrom S4Vectors metadata
#' @export
#' @return A list including \dQuote{pi} and \dQuote{theta}
#'
#' @author Ruizhu Huang
#' @examples
#'
#'
#' set.seed(1)
#' y <- matrix(rnbinom(100,size=1,mu=10),nrow=10)
#' colnames(y) <- paste("S", 1:10, sep = "")
#' rownames(y) <- tinyTree$tip.label
#'
#'
#' toy_tse <- TreeSummarizedExperiment(rowTree = tinyTree,
#'                                     assays = list(y))
#' res <- parEstimate(obj = toy_tse)
#'
#' metadata(res)$assays.par
#'
parEstimate <- function(obj) {

    stopifnot(any(class(obj) %in% c("matrix", "list",
                                 "TreeSummarizedExperiment")))

    if (is.matrix(obj)) {
        out <- .estimateA(obj = obj)
    }

    if (is.list(obj)) {
        out <- .estimateB(obj = obj)
    }

    if (is(obj, "TreeSummarizedExperiment")) {
        out <- .estimateC(obj = obj)
    }
    return(out)
}
