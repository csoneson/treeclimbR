#' Simulate different scenarios of abundance change in entities
#'
#' Simulate a data set with different abundance patterns for entities under
#' different conditions. These entities have their corresponding nodes on a
#' tree.
#'
#' @author Ruizhu Huang, Charlotte Soneson
#' @export
#'
#' @param tree A \code{phylo} object. Only used when \code{obj} is \code{NULL}.
#' @param data A count matrix with entities corresponding to tree leaves
#'     in the rows and samples in the columns. Only used when \code{obj} is
#'     \code{NULL}.
#' @param obj A \code{TreeSummarizedExperiment} object with observed data to
#'     use as the input for the simulation. If \code{NULL}, \code{data} and \
#'     \code{tree} must be provided instead.
#' @param assay If \code{obj} is not \code{NULL}, a numeric index or
#'     character scalar indicating which assay of the object to use as the
#'     basis for simulation. If \code{assay} is \code{NULL}, the first
#'     assay in the object is used.
#' @param scenario The simulation scenario, either \dQuote{BS}, \dQuote{US},
#'     or \dQuote{SS} (see \bold{Details}). The default is \dQuote{BS}.
#' @param from.A,from.B The branch node labels of branches A and B for which the
#'     signal will be swapped. By default, both are \code{NULL}, in which case
#'     they will be chosen based on the restrictions provided (\code{minTip.A},
#'     \code{maxTip.A}, \code{minTip.B}, \code{maxTip.B}, \code{minPr.A},
#'     \code{maxPr.A}, \code{ratio}). Note: If \code{from.A} is \code{NULL},
#'     \code{from.B} is also set to \code{NULL}.
#' @param minTip.A The minimum number of leaves allowed in branch A.
#' @param maxTip.A The maximum number of leaves allowed in branch A.
#' @param minTip.B The minimum number of leaves allowed in branch B.
#' @param maxTip.B The maximum number of leaves allowed in branch B.
#' @param minPr.A A numeric value in [0, 1]. The minimum abundance
#'     proportion of leaves in branch A.
#' @param maxPr.A A numeric value  in [0, 1]. The maximum abundance
#'     proportion of leaves in branch A.
#' @param ratio A numeric value. The proportion ratio of branch B to branch A.
#'     This value is used to select branches(see \bold{Details}). If there are
#'     no branches having exactly this ratio, the pair with the value closest to
#'     \code{ratio} will be selected.
#' @param adjB A numeric value in [0, 1] (only for \code{scenario}
#'     \dQuote{SS}), or \code{NULL}. If \code{NULL}, branch A and the selected
#'     part of branch B swap their proportions. If a numeric value, e.g. 0.1,
#'     then the counts for the selected part of branch B decreases to 10% of
#'     the original value, and this decrease is added to branch A. For example,
#'     assume there are two experimental conditions (C1 & C2), branch A has
#'     a count of 10 and branch B has a count of 40 in C1. If adjB is set to
#'     0.1, then in C2 branch B becomes 4 and branch A 46
#'     so that the total count of the two branches stays the same.
#' @param pct The percentage of leaves in branch B that have differential
#'     abundance under different conditions (only for scenario \dQuote{SS}).
#' @param nSam A numeric vector of length 2, indicating the sample size for each
#'     of the two simulated conditions.
#' @param mu,size The parameters of the Negative Binomial distribution. (see mu
#'     and size in \code{\link[stats:NegBinomial]{rnbinom}}). These parameters
#'     are used to generate the library size for each simulated sample. If
#'     \code{size} is not specified, \code{mu} should be a vector of numbers
#'     from which the library size is sampled with replacement.
#' @param n A numeric value to specify how many count tables would be generated
#'     with the same settings. The default is 1, i.e., one count table would be
#'     obtained at the end. If greater than 1, the output is a list of matrices.
#' @param FUN A function to calculate the aggregated count at each internal
#'     node based on its descendant leaves (e.g., \code{sum}, \code{mean}).
#'     The argument of the function should be a numeric vector with the counts
#'     of an internal node's descendant leaves.
#' @param message A logical scalar, indicating whether progress messages
#'     should be printed to the console.
#'
#' @returns a TreeSummarizedExperiment object.
#' \itemize{
#'     \item \strong{assays} A list of count matrices, with entities in rows and
#'     samples in columns. Each row can be mapped to a node of the tree.
#'     \item \strong{rowData} Annotation data for the rows.
#'     \item \strong{colData} Annotation data for the columns.
#'     \item \strong{rowTree} The tree structure of entities.
#'     \item \strong{rowLinks} The link between rows and nodes on the tree.
#'     \item \strong{metadata} More details about the simulation.
#'     \itemize{
#'         \item \strong{FC} the fold change of entities corresponding to the
#'         tree leaves.
#'         \item \strong{Branch} the information about two selected branches.
#'         \itemize{
#'             \item \strong{A} The branch node label (or number) of branch A.
#'             \item \strong{B} The branch node label (or number) of branch B.
#'             \item \strong{ratio} The count proportion ratio of branch B to
#'             branch A.
#'             \item \strong{A_tips} The number of leaves on branch A.
#'             \item \strong{B_tips} The number of leaves on branch B.
#'             \item \strong{A_prop} The count proportion of branch A.
#'             \item \strong{B_prop} The count proportion of branch B.
#'         }
#'     }
#' }
#'
#' @details Simulate a count table for entities which are corresponding to the
#'     nodes of a tree. The entities are in rows and the samples from different
#'     groups or conditions are in columns. The library size of each sample is
#'     sampled from a Negative Binomial distribution with mean and size
#'     specified by the arguments \code{mu} and \code{size}. The counts of
#'     entities, that are mapped to the leaf nodes, in a sample are assumed
#'     to follow a Dirichlet-Multinomial distribution. The parameters for
#'     the Dirichlet-Multinomial distribution are estimated from a real data set
#'     specified by \code{data} via the function \code{dirmult} (see
#'     \code{\link[dirmult]{dirmult}}). To generate different abundance patterns
#'     under different conditions, we provide three different scenarios,
#'     \dQuote{BS}, \dQuote{US}, and \dQuote{SS} (specified via
#'     \code{scenario}).
#'     \itemize{ \item BS: two branches are selected to swap their proportions,
#'     and leaves on the same branch have the same fold change.
#'     \item US: two branches are selected to swap their proportions. Leaves
#'     in the same branch have different fold changes but same
#'     direction (either increase or decrease).
#'     \item SS: two branches are selected. One branch has its proportion
#'     swapped with the proportion of some leaves from the other branch.}
#'
#' @importFrom dirmult dirmult
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom methods is
#' @importFrom SummarizedExperiment assays assayNames
#' @importFrom TreeSummarizedExperiment TreeSummarizedExperiment
#'
#' @examples
#' ## Generate data to use as the starting point (this would usually be a
#' ## real data set)
#' set.seed(1L)
#' y <- matrix(rnbinom(120, size = 1, mu = 10), nrow = 10)
#' colnames(y) <- paste("S", seq_len(12), sep = "")
#' rownames(y) <- tinyTree$tip.label
#'
#' toy_lse <- TreeSummarizedExperiment(rowTree = tinyTree,
#'                                     assays = list(counts = y))
#' simData(obj = toy_lse, ratio = 2, scenario = "BS", pct = 0.5)
#'
simData <- function(tree = NULL, data = NULL, obj = NULL, assay = NULL,
                    scenario = "BS", from.A = NULL, from.B = NULL,
                    minTip.A = 0, maxTip.A = Inf, minTip.B = 0, maxTip.B = Inf,
                    minPr.A = 0, maxPr.A = 1, ratio = 4, adjB = NULL,
                    pct = 0.6, nSam = c(50, 50), mu = 10000, size = NULL,
                    n = 1, FUN = sum, message = FALSE) {

    ## Check arguments
    ## -------------------------------------------------------------------------
    .assertVector(x = tree, type = "phylo", allowNULL = TRUE)
    .assertVector(x = obj, type = "TreeSummarizedExperiment", allowNULL = TRUE)
    stopifnot((!is.null(tree) && !is.null(data)) || !is.null(obj))
    stopifnot(is.null(assay) || length(assay) == 1)
    if (!is.null(obj) && is.character(assay)) {
        stopifnot(assay %in% SummarizedExperiment::assayNames(obj))
    }
    .assertScalar(x = scenario, type = "character",
                  validValues = c("BS", "SS", "US"))
    .assertScalar(x = minTip.A, type = "numeric")
    .assertScalar(x = maxTip.A, type = "numeric")
    .assertScalar(x = minTip.B, type = "numeric")
    .assertScalar(x = maxTip.B, type = "numeric")
    .assertScalar(x = minPr.A, type = "numeric", rngIncl = c(0, 1))
    .assertScalar(x = maxPr.A, type = "numeric", rngIncl = c(0, 1))
    .assertScalar(x = ratio, type = "numeric")
    .assertScalar(x = adjB, type = "numeric", rngIncl = c(0, 1),
                  allowNULL = TRUE)
    .assertScalar(x = pct, type = "numeric", rngIncl = c(0, 1),
                  allowNULL = TRUE)
    .assertVector(x = nSam, type = "numeric")
    .assertVector(x = mu, type = "numeric")
    .assertScalar(x = size, type = "numeric", allowNULL = TRUE)
    .assertScalar(x = n, type = "numeric")
    .assertVector(x = FUN, type = "function")
    .assertScalar(x = message, type = "logical")

    if (is.null(assay)) {
        assay <- 1
    }

    ## Simulate data
    ## -------------------------------------------------------------------------
    if (is.null(obj)) {
        ## Tree and data provided
        out <- .doData(tree = tree, data = data,
                       scenario = scenario,
                       from.A = from.A, from.B = from.B,
                       minTip.A = minTip.A, maxTip.A = maxTip.A,
                       minTip.B = minTip.B, maxTip.B = maxTip.B,
                       minPr.A = minPr.A, maxPr.A = maxPr.A,
                       ratio = ratio, adjB = adjB, pct = pct,
                       nSam = nSam, mu = mu, size = size,
                       n = n, FUN = FUN)
    } else {
        ## TreeSummarizedExperiment provided
        pars <- obj$pars
        if (is.null(pars)) {
            obj <- parEstimate(obj = obj, assay = assay)
            pars <- S4Vectors::metadata(obj)$assays.par
            tree <- TreeSummarizedExperiment::rowTree(obj)
        }

        out <- .doData(tree = tree, data = pars,
                       scenario = scenario,
                       from.A = from.A, from.B = from.B,
                       minTip.A = minTip.A, maxTip.A = maxTip.A,
                       minTip.B = minTip.B, maxTip.B = maxTip.B,
                       minPr.A = minPr.A, maxPr.A = maxPr.A,
                       ratio = ratio, adjB = adjB, pct = pct,
                       nSam = nSam, mu = mu, size = size,
                       n = n, FUN = FUN)
    }

    return(out)
}

#' Simulate a count table
#'
#' @author Ruizhu Huang
#' @keywords internal
#' @noRd
#'
#' @returns A list of objects:
#' \itemize{
#'     \item{FC}{The fold change of entities corresponding to the tree leaves.}
#'     \item{Count}{A list of count matrices or a count matrix.
#'     Entities on the row and samples in the column. Each count matrix includes
#'     entities corresponding to all nodes on the tree structure.}
#'     \item{Branch}{Information about two selected branches.}
#'     \describe{
#'         \item{A}{The branch node label of branch A}
#'         \item{B}{The branch node label of branch B}
#'         \item{ratio}{The count proportion ratio of branch B to branch A}
#'         \item{A_tips}{The number of leaves on branch A}
#'         \item{B_tips}{The number of leaves on branch B}
#'         \item{A_prop}{The count proportion of branch A (not above 1)}
#'         \item{B_prop}{The count proportion of branch B (not above 1)}
#'     }
#' }
#'
#' @importFrom dirmult dirmult
#'
.doData <- function(tree = NULL, data = NULL, scenario = "BS",
                    from.A = NULL, from.B = NULL, minTip.A = 0, maxTip.A = Inf,
                    minTip.B = 0, maxTip.B = Inf, minPr.A = 0, maxPr.A = 1,
                    ratio = 2, adjB = NULL, pct = 0.6, nSam = c(50, 50),
                    mu = 50, size = 10000, n = 1, FUN = sum) {

    ## Check arguments
    ## -------------------------------------------------------------------------
    .assertVector(x = tree, type = "phylo")
    ## data should be either a list (from parEstimate) or a count matrix
    if (!is.list(data)) {
        if (!is.matrix(data)) {
            stop("data should be a matrix.")
        } else {
            if (!setequal(rownames(data), tree$tip.label)) {
                stop("The rownames of data do not match the leaf labels.")
            }
        }
    }

    ## Check validity of from.A and from.B
    sn <- TreeSummarizedExperiment::showNode(tree, use.alias = TRUE)
    if (!is.null(from.A)) {
        if (is.character(from.A)) {
            is_alias <- from.A %in% names(sn)
            if (!is_alias) {
                ## Check if from.A is a node label and convert to alias
                is_label <- from.A %in% tree$node.label
                if (!is_label) {
                    stop("The provided from.A is not a node in the tree")
                } else {
                    ## Convert to alias
                    tmp <- TreeSummarizedExperiment::convertNode(
                        tree, node = from.A)
                    from.A <- names(sn)[match(tmp, sn)]
                }
            }
        } else if (is.numeric(from.A)) {
            ## Node number - check that it's in the tree
            if (!from.A %in% sn) {
                stop("The provided from.A is not a node in the tree")
            }
        } else {
            stop("from.A must be a character or numeric value")
        }
    }
    if (!is.null(from.B)) {
        if (is.character(from.B)) {
            is_alias <- from.B %in% names(sn)
            if (!is_alias) {
                ## Check if from.B is a node label and convert to alias
                is_label <- from.B %in% tree$node.label
                if (!is_label) {
                    stop("The provided from.B is not a node in the tree")
                } else {
                    ## Convert to alias
                    tmp <- TreeSummarizedExperiment::convertNode(
                        tree, node = from.B)
                    from.B <- names(sn)[match(tmp, sn)]
                }
            }
        } else if (is.numeric(from.B)) {
            ## Node number - check that it's in the tree
            if (!from.B %in% sn) {
                stop("The provided from.B is not a node in the tree")
            }
        } else {
            stop("from.B must be a character or numeric value")
        }
    }

    ## Estimate parameters for DM distribution
    ## -------------------------------------------------------------------------
    data <- parEstimate(obj = data)

    ## Pick branches
    ## -------------------------------------------------------------------------
    if (!is.null(from.A) && !is.null(from.B)) {
        pk <- .infLoc(tree = tree, data = data,
                      from.A = from.A, from.B = from.B)
    } else {
        pk <- .pickLoc(tree = tree, data = data,
                       from.A  = from.A, minTip.A = minTip.A,
                       maxTip.A = maxTip.A, minTip.B = minTip.B,
                       maxTip.B = maxTip.B, minPr.A = minPr.A,
                       maxPr.A = maxPr.A, ratio = ratio)
    }

    ## Estimate fold changes
    ## -------------------------------------------------------------------------
    beta <- .doFC(tree = tree, data = data, scenario = scenario,
                  branchA = pk$A, branchB = pk$B, ratio = pk$ratio,
                  adjB = adjB, pct = pct)

    ## Generate count matrix/matrices and TreeSummarizedExperiment objects
    ## -------------------------------------------------------------------------
    count <- .doCount(data = data, FC = beta, nSam = nSam, mu = mu,
                      size = size, n = n)

    if (is.list(count)) {
        grpDat <- data.frame(group = substr(colnames(count[[1]]), 1, 2))
        lse <- TreeSummarizedExperiment::TreeSummarizedExperiment(
            assays = count,
            metadata = list(FC = beta, branch = pk, scenario = scenario),
            colData = grpDat,
            rowTree = tree,
            rowNodeLab = rownames(count[[1]]))
    } else if (is.matrix(count)) {
        grpDat <- data.frame(group = substr(colnames(count), 1, 2))
        lse <- TreeSummarizedExperiment::TreeSummarizedExperiment(
            assays = list(counts = count),
            metadata = list(FC = beta, branch = pk, scenario = scenario),
            colData = grpDat,
            rowTree = tree,
            rowNodeLab = rownames(count))
    }

    return(lse)
}


#' Select branches to swap
#'
#' @author Ruizhu Huang
#' @keywords internal
#' @noRd
#'
#' @returns A \code{data.frame} with one row
#'
#' @importFrom TreeSummarizedExperiment convertNode findDescendant
#'
.pickLoc <- function(tree = NULL, data = NULL, from.A = NULL,
                     minTip.A = 0, maxTip.A = Inf, minTip.B = 0, maxTip.B = Inf,
                     minPr.A = 0, maxPr.A = 1, ratio = 1) {

    ## Estimate tip proportions, rename using the node label alias
    ## -------------------------------------------------------------------------
    pars <- parEstimate(obj = data)$pi
    nam1 <- names(pars)
    val1 <- TreeSummarizedExperiment::convertNode(
        tree = tree, node = nam1, message = FALSE)
    nam2 <- TreeSummarizedExperiment::convertNode(
        tree = tree, node = val1, use.alias = TRUE, message = FALSE)
    names(pars) <- nam2

    ## Proportions for internal nodes
    ## -------------------------------------------------------------------------
    leaf <- setdiff(tree$edge[, 2], tree$edge[, 1])
    leaf <- sort(leaf)
    nodI <- setdiff(tree$edge[, 1], leaf)
    nodI <- sort(nodI)
    desI <- TreeSummarizedExperiment::findDescendant(
        tree = tree, node = nodI, only.leaf = TRUE, self.include = TRUE,
        use.alias = TRUE)
    desI <- lapply(desI, FUN = function(x) {
        TreeSummarizedExperiment::convertNode(
            tree = tree, node = x, use.alias = TRUE, message = FALSE)
    })
    names(desI) <- TreeSummarizedExperiment::convertNode(
        tree = tree, node = nodI, use.alias = TRUE, message = FALSE)
    nodP <- mapply(function(x, y) {
        sum(x[y])
    }, x = list(pars), y = desI)

    ## Abundance proportions and number of descendant leaves
    ## -------------------------------------------------------------------------
    lenI <- unlist(lapply(desI, length))
    tt <- cbind(nodP, lenI)
    rownames(tt) <- TreeSummarizedExperiment::convertNode(
        tree = tree, node = nodI, use.alias = TRUE, message = FALSE)

    if (maxPr.A < min(tt[, 1])) {
        stop("maxPr.A is lower than the minimum value of
             node proportion", signif(min(tt[, 1]), 2), "\n")
    }
    if (minPr.A * ratio > max(tt[, 1])) {
        stop("minPr.A*ratio is above the maximum value of
             node proportion; try lower ratio", signif(min(tt[, 1]), 2), "\n")
    }

    ## Only consider nodes with enough tips and desired proportions
    ## -------------------------------------------------------------------------
    if (is.null(from.A)) {
        tt.sel <- tt
    } else {
        ## If node numbers, change them to node labels
        if (!is.character(from.A)) {
            from.A <- TreeSummarizedExperiment::convertNode(
                tree = tree, node = from.A, use.alias = TRUE, message = FALSE)
        }
        tt.sel <- tt[match(from.A, rownames(tt)), , drop = FALSE]
    }

    st <- tt.sel[tt.sel[, 2] >= minTip.A & tt.sel[, 2] <= maxTip.A &
                 tt.sel[, 1] >= minPr.A & tt.sel[, 1] <= maxPr.A, ,
                 drop = FALSE]
    if (nrow(st) == 0) {
        stop("No nodes fulfill the requirements; try other values for ",
             "minTip.A, maxTip.A, minPr.A, or maxPr.A")
    }
    st2 <- tt[tt[, 2] >= minTip.B & tt[, 2] <= maxTip.B, , drop = FALSE]

    ## Fold changes between any two nodes
    ## -------------------------------------------------------------------------
    mm <- (1/st[, 1]) %o% st2[, 1]
    rownames(mm) <- rownames(st)

    maxM <- max(mm, na.rm = TRUE)
    minM <- min(mm, na.rm = TRUE)

    if (ratio > 1 & ratio > maxM) {
        stop("Could not find two branches which fulfill the requirement; ",
             "try lower ratio, lower minTip.A, or higher maxTip.B\n")
    }
    if (ratio < 1 & ratio < minM) {
        stop("Could not find two branches which fulfill the requirement; ",
             "try higher ratio or lower minTip.B\n")
    }

    nm <- mm
    nm[] <- vapply(
        seq_len(ncol(mm)),
        FUN = function(x) {
            ## each column
            cn <- colnames(mm)
            cx <- cn[x]

            ## all rows
            rn <- rownames(mm)
            tx <- desI[rn]

            cs <- lapply(
                tx,
                FUN = function(x) {
                    length(intersect(x, desI[[cx]])) > 0
                })
            cv <- unlist(cs)
            fm <- mm[, x]
            fm[cv] <- NA
            fm
        },
        FUN.VALUE = numeric(nrow(mm))
    )

    colnames(nm) <- colnames(mm)
    rownames(nm) <- rownames(mm)

    dif <- abs(nm - ratio)
    wi <- which(dif == min(dif, na.rm = TRUE), arr.ind = TRUE)
    si <- sample(seq_len(nrow(wi)), 1)
    an <- rownames(nm)[wi[si, 1]]
    bn <- colnames(nm)[wi[si, 2]]

    ## Assemble information about selected branches
    ## -------------------------------------------------------------------------
    du <- cbind.data.frame(
        "A" = TreeSummarizedExperiment::convertNode(
            tree = tree, node = an, use.alias = FALSE, message = FALSE),
        "B" = TreeSummarizedExperiment::convertNode(
            tree = tree, node = bn, use.alias = FALSE, message = FALSE),
        "ratio" = unique(nm[wi]),
        "A_tips" = tt[an, 2],
        "B_tips" = tt[bn, 2],
        "A_prop" = round(tt[an, 1], digits = 4),
        "B_prop" = round(tt[bn, 1], digits = 4)
    )

    rownames(du) <- NULL
    return(du)
}

#' Provide information about proportions and number of leaves in two branches
#'
#' @keywords internal
#' @author Ruizhu Huang
#' @noRd
#'
#' @returns A \code{data.frame} with one row
#'
#' @importFrom TreeSummarizedExperiment convertNode findDescendant
#'
.infLoc <- function(tree = NULL, data = NULL,
                    from.A = NULL, from.B = NULL) {

    ## Estimate tip proportions, rename using the node label alias
    ## -------------------------------------------------------------------------
    pars <- parEstimate(obj = data)$pi
    nam1 <- names(pars)
    val1 <- TreeSummarizedExperiment::convertNode(
        tree = tree, node = nam1, message = FALSE)
    nam2 <- TreeSummarizedExperiment::convertNode(
        tree = tree, node = val1, use.alias = TRUE, message = FALSE)
    names(pars) <- nam2

    ## Proportions for internal nodes
    ## -------------------------------------------------------------------------
    leaf <- setdiff(tree$edge[, 2], tree$edge[, 1])
    leaf <- sort(leaf)
    nodI <- setdiff(tree$edge[, 1], leaf)
    nodI <- sort(nodI)
    nodA <- c(leaf, nodI)
    desI <- TreeSummarizedExperiment::findDescendant(
        tree = tree, node = nodI, only.leaf = TRUE, self.include = TRUE,
        use.alias = TRUE)
    desI <- lapply(desI, FUN = function(x) {
        TreeSummarizedExperiment::convertNode(
            tree = tree, node = x, use.alias = TRUE, message = FALSE)
    })
    names(desI) <- TreeSummarizedExperiment::convertNode(
        tree = tree, node = nodI, use.alias = TRUE, message = FALSE)
    nodP <- mapply(function(x, y) {
        sum(x[y])
    }, x = list(pars), y = desI)

    ## Abundance proportions and number of descendant leaves
    ## -------------------------------------------------------------------------
    lenI <- unlist(lapply(desI, length))
    tt <- cbind(nodP, lenI)
    rownames(tt) <- TreeSummarizedExperiment::convertNode(
        tree = tree, node = nodI, use.alias = TRUE, message = FALSE)

    ## Get branch names
    ## -------------------------------------------------------------------------
    labA <- ifelse(is.character(from.A), from.A,
                   TreeSummarizedExperiment::convertNode(
                       tree = tree, node = from.A, use.alias = TRUE,
                       message = FALSE))
    labB <- ifelse(is.character(from.B), from.B,
                   TreeSummarizedExperiment::convertNode(
                       tree = tree, node = from.B, use.alias = TRUE,
                       message = FALSE))

    ## Assemble summary information about nodes
    ## -------------------------------------------------------------------------
    rAB <- tt[labB, 1] / tt[labA, 1]
    du <- cbind.data.frame(
        "A" = from.A,
        "B" = from.B,
        "ratio" = rAB,
        "A_tips" = tt[labA, 2],
        "B_tips" = tt[labB, 2],
        "A_prop" = round(tt[labA, 1], digits = 4),
        "B_prop" = round(tt[labB, 1], digits = 4)
    )

    rownames(du) <- NULL
    return(du)
}

#' Generate fold changes
#'
#' @author Ruizhu Huang
#' @keywords internal
#' @noRd
#'
#' @returns a numeric vector with fold changes
#'
#' @importFrom stats runif
#'
.doFC <- function(tree = NULL, data = NULL, scenario = "BS",
                  branchA = NULL, branchB = NULL,
                  ratio = 1, adjB = NULL, pct = 1) {

    ## Find nodes
    ## -------------------------------------------------------------------------
    leaf <- setdiff(tree$edge[, 2], tree$edge[, 1])
    leaf <- sort(leaf)
    nodI <- setdiff(tree$edge[, 1], leaf)
    nodI <- sort(nodI)
    nodA <- c(leaf, nodI)

    ## Initialize beta
    ## -------------------------------------------------------------------------
    beta <- rep(1, length(leaf))
    names(beta) <- TreeSummarizedExperiment::convertNode(
        tree = tree, node = leaf, use.alias = TRUE)

    ## Node labels on branch A
    ## -------------------------------------------------------------------------
    ## leaves
    tip.A <- TreeSummarizedExperiment::findDescendant(
        tree = tree, node = branchA, only.leaf = TRUE, self.include = TRUE,
        use.alias = TRUE)
    tip.A <- TreeSummarizedExperiment::convertNode(
        tree = tree, node = unlist(tip.A), use.alias = TRUE, message = FALSE)

    ## nodes
    nodA.A <- TreeSummarizedExperiment::findDescendant(
        tree = tree, node = branchA, only.leaf = FALSE, self.include = TRUE,
        use.alias = TRUE)
    nodA.A <- TreeSummarizedExperiment::convertNode(
        tree = tree, node = unlist(nodA.A), use.alias = TRUE, message = FALSE)

    ## internal nodes
    nodI.A <- setdiff(nodA.A, tip.A)

    ## descendants of internal nodes
    des.IA <- TreeSummarizedExperiment::findDescendant(
        tree = tree, node = nodI.A, only.leaf = TRUE, self.include = TRUE,
        use.alias = TRUE)
    des.IA <- lapply(des.IA, FUN = function(x) {
        TreeSummarizedExperiment::convertNode(
            tree = tree, node = x, use.alias = TRUE, message = FALSE)
    })

    ## Tip proportions from real data, rename using the alias of node label
    ## -------------------------------------------------------------------------
    pars <- parEstimate(obj = data)$pi
    nam1 <- names(pars)
    val1 <- TreeSummarizedExperiment::convertNode(
        tree = tree, node = nam1, message = FALSE)
    nam2 <- TreeSummarizedExperiment::convertNode(
        tree = tree, node = val1, use.alias = TRUE, message = FALSE)
    names(pars) <- nam2

    ## Swap proportions of two branches
    ## -------------------------------------------------------------------------
    if (scenario == "BS") {
        ## Tips in the same branch have the same fold change
        ## leaves on branch B
        tip.B <- TreeSummarizedExperiment::findDescendant(
            tree = tree, node = branchB, only.leaf = TRUE, self.include = TRUE,
            use.alias = TRUE)
        tip.B <- TreeSummarizedExperiment::convertNode(
            tree = tree, node = unlist(tip.B), use.alias = TRUE,
            message = FALSE)
        beta[tip.A] <- ratio
        beta[tip.B] <- 1/ratio
    } else if (scenario == "US") {
        ## Tips in the same branch have different fold changes but same
        ## direction (either increase or decrease)
        tip.B <- TreeSummarizedExperiment::findDescendant(
            tree = tree, node = branchB, only.leaf = TRUE, self.include = TRUE,
            use.alias = TRUE)
        tip.B <- TreeSummarizedExperiment::convertNode(
            tree = tree, node = unlist(tip.B), use.alias = TRUE,
            message = FALSE)

        ## Proportion in the two branches
        propA <- sum(pars[tip.A])
        propB <- sum(pars[tip.B])

        ## Make sure the ratio is above zero
        if (propB < propA) {
            rA <- 1 - propB / propA
            rB <- 0
        } else {
            rA <- 0
            rB <- 1 - propA / propB
        }

        a1 <- stats::runif(length(tip.A), rA, 1)
        sa <- sum(a1*pars[tip.A])
        a2 <- (propB - propA)/sa
        a3 <- a1 * a2 + 1
        beta[tip.A] <- a3

        b1 <- stats::runif(length(tip.B), rB, 1)
        sb <- sum(b1 * pars[tip.B])
        b2 <- (propA - propB)/sb
        b3 <- b1 * b2 + 1
        beta[tip.B] <- b3
    } else if (scenario == "SS") {
        ## Distribute signal randomly in one branch and evenly in
        ## another branch
        iter <- 1
        while (iter <= 200) {
            ## Select only some tips
            selA <- sample(tip.A, ceiling(length(tip.A) * pct))
            ## Make the selected tips disperse evenly in the branch
            subA <- lapply(des.IA, FUN = function(x) {
                ix <- intersect(x, selA)
                length(ix)/length(x)
            })
            subA <- unlist(subA)
            sumA <- sum(pars[selA])

            ## The abundance proportions of the selected tips
            ## are roughly equal to its proportion in the branch.
            ## Avoid selecting all low or high abundance tips in the branch.

            spr <- sumA / sum(pars[tip.A])
            ind.pr <- spr <= (pct + 0.05) & spr >= (pct - 0.05)

            if (all(subA <= 0.6) & ind.pr) {
                break
            }
            iter <- iter + 1
        }

        if (length(branchB) == 0) {
            stop("No suitable branches.
                 Try another branchA or another max ratio... \n")
        }
        tip.B <- TreeSummarizedExperiment::findDescendant(
            tree = tree, node = branchB, only.leaf = TRUE, self.include = TRUE,
            use.alias = TRUE)
        tip.B <- TreeSummarizedExperiment::convertNode(
            tree = tree, node = unlist(tip.B), use.alias = TRUE,
            message = FALSE)
        sumB <- sum(pars[tip.B])

        if (is.null(adjB)) {
            beta[selA] <- sumB/sumA
            beta[tip.B] <- sumA/sumB
        } else {
            if (!is.numeric(adjB)) {
                stop("adjB should be numeric")
            }
            beta[tip.B] <- adjB
            beta[selA] <- (sumB * (1 - adjB) + sumA) / sumA
        }
    }

    ## Rename beta with the node label instead of the alias of node label
    ## -------------------------------------------------------------------------
    names(beta) <- TreeSummarizedExperiment::convertNode(
        tree = tree, node = leaf, use.alias = FALSE)
    return(beta)
}


#' Generate a count table given simulation details
#'
#' @author Ruizhu Huang
#' @keywords internal
#' @noRd
#'
#' @returns A matrix or a list of matrices.
#'
#' @importFrom dirmult rdirichlet
#' @importFrom stats rmultinom rnbinom
#'
.doCount <- function(data, FC, nSam, mu, size, n) {

    ## Parameters
    ## -------------------------------------------------------------------------
    pars <- parEstimate(obj = data)
    theta <- pars$theta
    gplus <- (1 - theta) / theta

    ## Tip proportions
    ## -------------------------------------------------------------------------
    pr <- pars$pi
    p.c1 <- pr
    p.c2 <- pr * FC[names(p.c1)]

    ## Parameters for Dirichlet distribution
    ## -------------------------------------------------------------------------
    g.c1 <- p.c1 * gplus
    g.c2 <- p.c2 * gplus

    ## Generate count matrices
    ## -------------------------------------------------------------------------
    resList <- lapply(seq_len(n), FUN = function(j) {
        ## condition 1
        n1 <- nSam[1]
        Mp.c1 <- matrix(0, nrow = n1, ncol = length(g.c1))
        rownames(Mp.c1) <- paste("C1_", seq_len(n1), sep = "")
        colnames(Mp.c1) <- names(p.c1)
        Mobs.c1 <- Mp.c1
        if (length(mu) & !length(size)) {
            nSeq1 <- sample(x = mu, size = n1, replace = TRUE)
        } else {
            nSeq1 <- stats::rnbinom(n = n1, mu = mu, size = size)
        }
        for (i in seq_len(n1)) {
            Mp.c1[i, ] <- dirmult::rdirichlet(1, g.c1)[1, ]
            Mobs.c1[i, ] <- stats::rmultinom(1, nSeq1[i],
                                             prob = Mp.c1[i, ])[, 1]
        }

        ## condition 2
        n2 <- nSam[2]
        Mp.c2 <- matrix(0, nrow = n2, ncol = length(g.c2))
        rownames(Mp.c2) <- paste("C2_", seq_len(n2), sep = "")
        colnames(Mp.c2) <- names(p.c2)
        Mobs.c2 <- Mp.c2
        if (length(mu) & !length(size)) {
            nSeq2 <- sample(x = mu, size = n2, replace = TRUE)
        } else {
            nSeq2 <- stats::rnbinom(n = n2, mu = mu, size = size)
        }

        for (i in seq_len(n2)) {
            Mp.c2[i, ] <- dirmult::rdirichlet(1, g.c2)[1, ]
            Mobs.c2[i, ] <- stats::rmultinom(1, nSeq2[i],
                                             prob = Mp.c2[i, ])[, 1]
        }

        cb <- t(rbind(Mobs.c1, Mobs.c2))

        return(cb)
    })

    if (n == 1) {
        count <- do.call(rbind, resList)
    } else {
        count <- resList
    }

    return(count)
}

