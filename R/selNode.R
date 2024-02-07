#' Select branches meeting certain criteria
#'
#' Select branches in a tree meeting the specified criteria in terms of number
#' of leaves and the count proportion. Note that only internal branch nodes
#' are considered - no individual leaves will be returned.
#'
#' @author Ruizhu Huang, Charlotte Soneson
#' @export
#'
#' @param pr A named numeric vector to provide proportions of entities. If
#'     this is provided, \code{obj} and \code{data} will be ignored.
#' @param obj A \code{TreeSummarizedExperiment} object. Only used if \code{pr}
#'     is \code{NULL}.
#' @param assay The index or name of the assay of \code{obj} to use for
#'     estimating node count proportions. Only used if \code{obj} is not
#'     \code{NULL}.
#' @param tree A \code{phylo} object. If \code{obj} is used as input, the tree
#'     will be extracted from the \code{rowTree} of \code{obj}.
#' @param data Either a count table with entities in rows and samples in
#'     columns, or a list with \code{pi} and \code{theta} estimates (the
#'     output of \code{\link{parEstimate}}). Only used if \code{pr} and
#'     \code{obj} are \code{NULL}.
#' @param minTip the minimum number of leaves in the selected branch.
#' @param maxTip The maximum number of leaves in the selected branch.
#' @param minPr The minimum count proportion of the selected branch in a sample.
#'     A value between 0 and 1.
#' @param maxPr The maximum count proportion of the selected branch in a sample.
#'     A value between 0 and 1.
#' @param skip A character vector of node labels. These nodes can not be
#'     descendants or the ancestors of the selected branch.
#' @param all A logical scalar. If \code{FALSE} (default), the branch node of a
#'     single branch, which meets the requirements and has the minimum count
#'     proportion of branches meeting the requirements, is returned; otherwise
#'     branch nodes of all branches meeting the requirements are returned.
#'
#' @returns A \code{data.frame} with node information for the selected
#'     internal node(s).
#'
#' @importFrom methods is
#' @importFrom TreeSummarizedExperiment rowTree convertNode findDescendant
#' @importFrom S4Vectors metadata
#'
#' @examples
#' ## Generate example data
#' library(TreeSummarizedExperiment)
#' set.seed(1)
#' data(tinyTree)
#' toyTable <- matrix(rnbinom(40, size = 1, mu = 10), nrow = 10)
#' colnames(toyTable) <- paste(rep(LETTERS[seq_len(2)], each = 2),
#'                             rep(seq_len(2), 2), sep = "_")
#' rownames(toyTable) <- tinyTree$tip.label
#'
#' ## Estimate entity proportions from count matrix under a Dirichlet
#' ## Multinomial framework, and use this as the input for selNode
#' dat <- parEstimate(obj = toyTable)
#' selNode(tree = tinyTree, data = dat, all = TRUE)
#' selNode(tree = tinyTree, data = dat,
#'         minTip = 4, maxTip = 9, minPr = 0, maxPr = 0.8, all = TRUE)
#'
#' ## Alternatively, directly provide the proportions vector
#' selNode(tree = tinyTree, pr = dat$pi, all = TRUE)
#'
#' ## Return only branch with lowest proportion among valid ones
#' selNode(tree = tinyTree, pr = dat$pi, all = FALSE)
#'
#' ## Start instead from a TreeSummarizedExperiment object
#' lse <- TreeSummarizedExperiment(rowTree = tinyTree,
#'                                 assays = list(counts = toyTable))
#' selNode(obj = lse, assay = "counts", all = TRUE)
#'
#' ## Don't allow node 1 to be included
#' selNode(obj = lse, assay = "counts", skip = 1, all = TRUE)
#'
selNode <- function(pr = NULL, obj = NULL, assay = 1, data = NULL, tree = NULL,
                    minTip = 0, maxTip = Inf, minPr = 0, maxPr = 1,
                    skip = NULL, all = FALSE) {

    ## Check arguments and extract tree and proportions vector
    ## -------------------------------------------------------------------------
    .assertScalar(x = minTip, type = "numeric", rngIncl = c(0, Inf))
    .assertScalar(x = maxTip, type = "numeric", rngIncl = c(0, Inf))
    .assertScalar(x = minPr, type = "numeric", rngIncl = c(0, 1))
    .assertScalar(x = maxPr, type = "numeric", rngIncl = c(0, 1))
    .assertScalar(x = all, type = "logical")

    if (!is.null(pr)) {
        ## If pr is provided, use that
        ## pr must be a numeric proportions vector, and the tree must exist
        if (!is.null(obj)) {
            message("Ignoring obj when input is a proportions vector")
        }
        if (!is.null(data)) {
            message("Ignoring data when input is a proportions vector")
        }
        .assertVector(x = pr, type = "numeric", rngIncl = c(0, 1))
        .assertVector(x = tree, type = "phylo")
        pars <- pr
    } else {
        ## If pr is not provided, first check for a TSE
        if (methods::is(obj, "TreeSummarizedExperiment")) {
            if (!is.null(tree)) {
                message("Ignoring tree when input is a TSE")
            }
            if (!is.null(data)) {
                message("Ignoring data when input is a TSE")
            }
            ## Estimate parameters
            obj <- parEstimate(obj = obj, assay = assay)
            tree <- TreeSummarizedExperiment::rowTree(obj)
            data <- S4Vectors::metadata(obj)$assays.par
            pars <- data$pi
        } else if (!is.null(data)) {
            ## pr not provided, obj not provided -> use data
            .assertVector(x = tree, type = "phylo")
            tree <- tree
            data <- data
            pars <- parEstimate(obj = data)$pi
        } else {
            stop("No valid input was found - please provide either pr + tree, ",
                 "obj, or data + tree.")
        }
    }

    ## Proportions vector must be named and have one entry for each tree tip
    stopifnot(!is.null(names(pars)))
    stopifnot(length(pars) == length(tree$tip.label))
    stopifnot(all(names(pars) %in% tree$tip.label))

    ## Get number of descendant leaves for internal nodes
    ## -------------------------------------------------------------------------
    leaf <- setdiff(tree$edge[, 2], tree$edge[, 1])
    nodI <- setdiff(tree$edge[, 1], leaf)

    nodeLab <- TreeSummarizedExperiment::convertNode(
        tree = tree, node = nodI, use.alias = FALSE, message = FALSE)
    nodeLab_alias <- TreeSummarizedExperiment::convertNode(
        tree = tree, node = nodI, use.alias = TRUE, message = FALSE)

    ## Get descendant leaves
    tipI <- TreeSummarizedExperiment::findDescendant(
        tree = tree, node = nodI, only.leaf = TRUE, self.include = TRUE)
    names(tipI) <- nodeLab_alias

    ## Get number of descendant leaves
    numI <- unlist(lapply(tipI, length))

    ## Get proportions for internal nodes
    ## -------------------------------------------------------------------------
    vnum <- TreeSummarizedExperiment::convertNode(
        tree = tree, node = names(pars), message = FALSE)

    ## Proportion for each node
    nodP <- vapply(tipI, FUN = function(x) {
        sum(pars[match(x, vnum)])
    }, FUN.VALUE = NA_real_)

    ## Sample
    ## -------------------------------------------------------------------------
    if (any(duplicated(nodeLab))) {
        tt <- cbind.data.frame(
            nodeNum = TreeSummarizedExperiment::convertNode(tree = tree,
                                                            node = names(nodP),
                                                            message = FALSE),
            nodeLab = nodeLab,
            nodeLab_alias = nodeLab_alias,
            proportion = nodP,
            numTip = numI)
    } else {
        tt <- cbind.data.frame(
            nodeNum = TreeSummarizedExperiment::convertNode(tree = tree,
                                                            node = names(nodP),
                                                            message = FALSE),
            nodeLab = nodeLab,
            proportion = nodP,
            numTip = numI)
    }

    if (maxPr < min(tt$proportion)) {
        stop("All nodes have proportions exceeding maxPr. The lowest node ",
             "proportion is ", signif(min(tt$proportion), 2))
    }

    ## Only consider nodes with enough tips and desired proportion level
    st <- tt[tt$numTip >= minTip & tt$numTip <= maxTip &
                 tt$proportion >= minPr & tt$proportion <= maxPr, ]
    if (nrow(st) == 0) {
        stop("No nodes fulfill the requirements; try other settings ",
             "for tip numbers or proportions.")
    }

    ## Remove nodes containing nodes that should be skipped
    if (!is.null(skip)) {
        if (is.character(skip)) {
            skip <- TreeSummarizedExperiment::convertNode(
                tree = tree, node = skip, message = FALSE)
        }
        tipS <- TreeSummarizedExperiment::findDescendant(
            tree = tree, node = skip, only.leaf = TRUE, self.include = TRUE)
        tipS <- unlist(tipS)

        rmp <- vapply(st$nodeNum, FUN = function(x) {
            tx <- TreeSummarizedExperiment::findDescendant(
                node = x, tree = tree, only.leaf = TRUE, self.include = TRUE)
            ix <- intersect(tipS, unlist(tx))
            length(ix) == 0
        }, FUN.VALUE = TRUE)
        new.st <- st[rmp, ]
    } else {
        new.st <- st
    }

    if (nrow(new.st) == 0) {
        stop("No nodes fulfill the requirements; try other settings ",
             "for tip numbers or proportions.")
    }

    ## Return the branch with the lowest proportion if all = FALSE
    ## -------------------------------------------------------------------------
    if (all) {
        final <- new.st
    } else {
        final <- new.st[which.min(new.st$proportion), ]
    }

    return(final)
}
