#' Simulate different scenarios of abundance change in entities
#'
#' \code{simData} simulates different abundance patterns for entities under
#' different conditions. These entities have their corresponding nodes on a
#' tree. More details about the simulated patterns could be found in the
#' vignette via \code{browseVignettes("treeclimbR")}.
#'
#' @param tree A phylo object. Only use when \code{obj} is NULL.
#' @param data A matrix, representing a table of values, such as count,
#'   collected from real data. It has the entities corresponding to tree leaves
#'   in the row and samples in the column. Only use when \code{obj} is NULL.
#' @param obj A TreeSummarizedExperiment object. More details in
#'   \code{\link[TreeSummarizedExperiment]{TreeSummarizedExperiment-class}}. The
#'   \code{tree} (see argument \code{tree}) is stored as \code{rowTree} and
#'   \code{data} (see argument \code{data}) is stored as the first table of
#'   \code{assays}.
#' @param scenario \dQuote{BS}, \dQuote{US}, or \dQuote{SS} (see
#'   \bold{Details}). Default is \dQuote{BS}.
#' @param from.A,from.B The branch node labels of branches A and B for which the
#'   signal is swapped. Default, both are NULL. In simulation, we select two
#'   branches (A & B) to have differential abundance under different conditions.
#'   One could specify these two branches or let \code{simData} choose. (Note:
#'   If \code{from.A} is NULL, \code{from.B} is set to NULL).
#' @param minTip.A The minimum number of leaves in branch A
#' @param maxTip.A The maximum number of leaves in branch A
#' @param minTip.B The minimum number of leaves in branch B
#' @param maxTip.B The maximum number of leaves in branch B
#' @param minPr.A A numeric value selected from 0 to 1. The minimum abundance
#'   proportion of leaves in branch A
#' @param maxPr.A A numeric value selected from 0 to 1. The maximum abundance
#'   proportion of leaves in branch A
#' @param ratio A numeric value. The proportion ratio of branch B to branch A.
#'   This value is used to select branches(see \bold{Details}). If there are no
#'   branches having exactly this ratio, the pair with the value closest to
#'   \code{ratio} would be selected.
#' @param adjB a numeric value selected from 0 and 1 (only for \code{scenario}
#'   is \dQuote{SS}). Default is NULL. If NULL, branch A and the selected part
#'   of branch B swap their proportions. If a numeric value, e.g. 0.1, then the
#'   selected part of branch B decreases to its one tenth proportion and the
#'   decrease in branch B is added to branch A. For example, assume there are
#'   two experimental conditions (C1 & C2), branch A has 10 and branch B has 40
#'   in C1. If adjB is set to 0.1, then in C2 branch B becomes 4 and branch A 46
#'   so that the total proportion stays the same.
#' @param pct The percentage of leaves in branch B that have differential
#'   abundance under different conditions (only for scenario \dQuote{SS})
#' @param nSam A numeric vector of length 2, containing the sample size for two
#'   different conditions
#' @param mu,size The parameters of the Negative Binomial distribution. (see mu
#'   and size in \code{\link[stats:NegBinomial]{rnbinom}}). Parameters used to
#'   generate the library size for each simulated sample. If \code{size} is not
#'   specified, \code{mu} should be a vector of numbers from which the library
#'   size is sampled from with replacement.
#' @param n A numeric value to specify how many count tables would be generated
#'   with the same settings. Default is one and one count table would be
#'   obtained at the end. If above one, the output is a list of matrices (count
#'   tables). This is useful, when one needs multiple simulations.
#' @param FUN A function to derive the count at each internal node based on its
#'   descendant leaves, e.g. sum, mean. The argument of the function is a
#'   numeric vector with the counts of an internal node's descendant leaves.
#' @param message A logical value, TRUE or FALSE. Default is FALSE. If TRUE, the
#'   progress message is printed out.
#'
#' @importFrom dirmult dirmult
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom methods is
#' @importFrom SummarizedExperiment assays
#' @export
#'
#' @return a TreeSummarizedExperiment object.
#'  \itemize{
#'  \item \strong{assays} a list of count table or a count table. Entities on
#'  the row and samples in the column. Each row could be mapped to a node of the
#'  tree.
#'  \item \strong{rowData} the annotation data for the rows of tables in
#'  \code{assays}.
#'  \item \strong{colData} the annotation data for the columns of tables in
#'  \code{assays}.
#'  \item \strong{rowTree} the tree structure of entities.
#'  \item \strong{rowLinks} the link between rows of \code{assays} and nodes on
#'  the tree.
#'  \item \strong{metadata} more details about the simulation.
#'    \itemize{
#'    \item \strong{FC} the fold change of entities correspondint to the tree
#'    leaves.
#'    \item \strong{Branch} the information about two selected branches.
#'      \itemize{
#'        \item \strong{A} the branch node label (or number) of branch A
#'        \item \strong{B} the branch node label (or number) of branch B
#'        \item \strong{ratio} the count proportion ratio of branch B to branch
#'        A
#'        \item \strong{A_tips} the number of leaves on branch A
#'        \item \strong{B_tips} the number of leaves on branch B
#'        \item \strong{A_prop} the count proportion of branch A
#'        (a value between 0 and 1)
#'        \item \strong{B_prop} the count proportion of branch B
#'        (a value between 0 and 1)
#'      }
#'      }}
#'
#' @details \code{simData} simulates a count table for entities which are
#'   corresponding to the nodes of a tree. The entities are in rows and the
#'   samples from different groups or conditions are in columns. The library
#'   size of each sample is sampled from a Negative Binomial distribution with
#'   mean and size specified by the arguments \code{mu} and \code{size}. The
#'   counts of entities, that are mapped ot the leaf nodes, in a same sample are
#'   assumed to follow a Dirichlet-Multinomial distribution. The parameters for
#'   the Dirichlet-Multinomial distribution are estimated from a real data set
#'   specified by the argument \code{data} via the function \code{dirmult} (see
#'   \code{\link[dirmult]{dirmult}}). To generate different abundance patterns
#'   under different conditions, we provide three different scenarios,
#'   \dQuote{BS}, \dQuote{US}, and \dQuote{SS} (specified via \code{scenario}).
#'   Our vignette provides figures to explain these three scenarios (try
#'   \code{browseVignettes("treeAGG")}). \itemize{ \item BS: two branches are
#'   selected to swap their proportions, and leaves on the same branch have the
#'   same fold change. \item US: two branches are selected to swap their
#'   proportions. Leaves in the same branch have different fold changes but same
#'   direction (either increase or decrease). \item SS: two branches are
#'   selected. One branch has its proportion swapped with the proportion of some
#'   leaves from the other branch.}
#'
#' @author Ruizhu Huang
#'
#' @examples
#'
#' set.seed(1)
#' y <- matrix(rnbinom(100,size=1,mu=10),nrow=10)
#' colnames(y) <- paste("S", 1:10, sep = "")
#' rownames(y) <- tinyTree$tip.label
#'
#'
#' toy_lse <- TreeSummarizedExperiment(rowTree = tinyTree,
#'                                     assays = list(y))
#' res <- parEstimate(obj = toy_lse)
#'
#' set.seed(1122)
#' dat1 <- simData(obj = res, ratio = 2, scenario = "BS", pct = 0.5)
#'
#'
#'


simData <- function(tree = NULL, data = NULL,
                    obj = NULL, scenario = "BS",
                    from.A = NULL, from.B = NULL,
                    minTip.A = 0, maxTip.A = Inf,
                    minTip.B = 0, maxTip.B = Inf,
                    minPr.A = 0, maxPr.A = 1,
                    ratio = 4, adjB = NULL,
                    pct = 0.6, nSam = c(50, 50),
                    mu = 10000, size = NULL,
                    n = 1, FUN = sum, message = FALSE){

    # -------------------------------------------------------------------------
    # provide (tree & data)
    if (is.null(obj)) {
        if (is.null(tree) | is.null(data)) {
            stop("tree or data is not provided")
        } else {
            obj <- .doData(tree = tree, data = data,
                           scenario = scenario,
                           from.A = from.A, from.B = from.B,
                           minTip.A = minTip.A, maxTip.A = maxTip.A,
                           minTip.B = minTip.B, maxTip.B = maxTip.B,
                           minPr.A = minPr.A, maxPr.A = maxPr.A,
                           ratio = ratio, adjB = adjB, pct = pct,
                           nSam = nSam, mu = mu, size = size,
                           n = n, FUN = FUN) }

        # ----------------------------------------------------------------------
        # provide a TreeSummarizedExperiment object
    } else {
        pars <- obj$pars
        if (is.null(pars)) {
            obj <- parEstimate(obj = obj)
            pars <- metadata(obj)$assays.par
            tree <- rowTree(obj)

            if(length(assays(obj)) > 1){
                message("\n more than one table provided in the assays;
                        only the first one would be used. \n")}
            data <- assays(obj)[[1]]
            }

        obj <- .doData(tree = tree, data = pars,
                       scenario = scenario, from.A = from.A,
                       from.B = from.B,
                       minTip.A = minTip.A, maxTip.A = maxTip.A,
                       minTip.B = minTip.B, maxTip.B = maxTip.B,
                       minPr.A = minPr.A, maxPr.A = maxPr.A,
                       ratio = ratio, adjB = adjB, pct = pct,
                       nSam = nSam, mu = mu, size = size,
                       n = n, FUN = FUN)

    }
    return(obj)
}



#' Simulate a count table
#'
#' \code{.doData} creates a count table for all nodes of a tree under two
#' different groups such that the tree would have different abundance patterns
#' in the different conditions.
#'
#' @param tree A phylo object
#' @param data A matrix, representing a count table from real data. It has the
#'   entities corresponding to tree leaves in the row and samples in the column.
#' @param scenario \dQuote{BS}, \dQuote{US}, or \dQuote{SS} (see
#'   \bold{Details}). Default is \dQuote{BS}.
#' @param from.A,from.B The branch node labels of branches A and B for which the
#'   signal is swapped. Default, both are NULL. In simulation, we select two
#'   branches (A & B) to have differential abundance under different conditions.
#'   One could specify these two branches or let \code{.doData} choose. (Note:
#'   If \code{from.A} is NULL, \code{from.B} is set to NULL).
#' @param minTip.A The minimum number of leaves in branch A
#' @param maxTip.A The maximum number of leaves in branch A
#' @param minTip.B The minimum number of leaves in branch B
#' @param maxTip.B The maximum number of leaves in branch B
#' @param minPr.A A numeric value selected from 0 to 1. The minimum abundance
#'   proportion of leaves in branch A
#' @param maxPr.A A numeric value selected from 0 to 1.The maximum abundance
#'   proportion of leaves in branch A
#' @param ratio The proportion ratio of branch B to branch A. This value is used
#'   to select branches(see \bold{Details}). If there are no branches having
#'   exactly this ratio, the pair with the value closest to \code{ratio} would
#'   be selected.
#' @param adjB a numeric value selected from 0 and 1 (only for \code{scenario}
#'   is \dQuote{SS}). Default is NULL. If NULL, branch A and branch B swap their
#'   proportions. If a numeric value, e.g. 0.1, then branch B decreases to its
#'   one tenth proportion and the decrease in branch B is added to branch A. For
#'   example, assume there are two experimental conditions (C1 & C2), branch A
#'   has 10 and branch B has 40 in C1. If adjB is set to 0.1, then in C2 branch
#'   B becomes 4 and branch A 46 so that the total proportion stays the same.
#' @param pct a numeric value selected from 0 and 1. The percentage of leaves in
#'   branch B that have differential abundance under different conditions (only
#'   for scenario \dQuote{SS})
#' @param nSam A numeric vector of length 2, containing the sample size for two
#'   different conditions
#' @param mu,size The parameters of the Negative Binomial distribution. (see mu
#'   and size in \code{\link[stats:NegBinomial]{rnbinom}}). Parameters used to
#'   generate the library size for each simulated sample.
#' @param n A numeric value to specify how many count tables would be generated
#'   with the same settings. Default is one and one count table would be
#'   obtained at the end. If above one, the output of \code{.doData} is a list
#'   of matrices (count tables). This is useful, when one needs multiple
#'   simulations.
#' @param FUN A function to derive the count at each internal node based on its
#'   descendant leaves, e.g. sum, mean. The argument of the function is a
#'   numeric vector with the counts of an internal node's descendant leaves.
#'
#' @importFrom dirmult dirmult
#' @keywords internal
#'
#' @return a list of objects \item{FC}{the fold change of entities correspondint
#'   to the tree leaves.} \item{Count}{a list of count table or a count table.
#'   Entities on the row and samples in the column. Each count table includes
#'   entities corresponding to all nodes on the tree structure.}
#'   \item{Branch}{the information about two selected branches.} \describe{
#'   \item{A}{the branch node label of branch A} \item{B}{the branch node label
#'   of branch B} \item{ratio}{the count proportion ratio of branch B to branch
#'   A} \item{A_tips}{the number of leaves on branch A} \item{B_tips}{the number
#'   of leaves on branch B} \item{A_prop}{the count proportion of branch A (not
#'   above 1)} \item{B_prop}{the count proportion of branch B (not above 1)} }
#'
#' @details \code{.doData} simulates a count table for entities which are
#'   corresponding to the nodes of a tree. The entities are in rows and the
#'   samples from different groups or conditions are in columns. The library
#'   size of each sample is sampled from a Negative Binomial distribution with
#'   mean and size specified by the arguments \code{mu} and \code{size}. The
#'   counts of entities, which are located on the tree leaves, in the same
#'   sample are assumed to follow a Dirichlet-Multinomial distribution. The
#'   parameters for the Dirichlet-Multinomial distribution are estimated from a
#'   real data set specified by the argument \code{data} via the function
#'   \code{dirmult} (see \code{\link[dirmult]{dirmult}}). To generate different
#'   abundance patterns under different conditions, we provide three different
#'   scenarios, \dQuote{BS}, \dQuote{US}, and \dQuote{SS} (specified via
#'   \code{scenario}). \itemize{ \item BS: two branches are selected to swap
#'   their proportions, and leaves on the same branch have the same fold change.
#'   \item US: two branches are selected to swap their proportions. Leaves in
#'   the same branch have different fold changes but same direction (either
#'   increase or decrease). \item SS: two branches are selected. One branch has
#'   its proportion swapped with the proportion of some leaves from the other
#'   branch.}
#' @author Ruizhu Huang
#'
#' @examples
#' \dontrun{
#' if(require(GUniFrac)){
#' data("throat.otu.tab")
#' data("throat.tree")
#'
#' dat <- .doData(tree = throat.tree,
#' data = as.matrix(t(throat.otu.tab)),
#' ratio = 2)
#' }
#'}

.doData <- function(tree = NULL, data = NULL,
                    scenario = "BS",
                    from.A = NULL, from.B = NULL,
                    minTip.A = 0, maxTip.A = Inf,
                    minTip.B = 0, maxTip.B = Inf,
                    minPr.A = 0, maxPr.A = 1,
                    ratio = 2, adjB = NULL,
                    pct = 0.6, nSam = c(50, 50),
                    mu = 50, size = 10000,
                    n = 1, FUN = sum){

    # ---check input is in correct format --------
    if(!is(tree, "phylo")){
        stop("tree requires a phylo object")
    }

    if (!is.list(data)) {
        if (!is.matrix(data)) {
            stop("data requires a matrix")
        } else {
            if (!setequal(rownames(data), tree$tip.label)) {
                stop("The rownames of data mismatch with the leaf labels")
            }
        }
    }


    # estimate parameters for Dirichlet-multinomial distribution
    data <- parEstimate(obj = data)

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


    beta <- .doFC(tree = tree, data = data,
                  scenario = scenario,
                  branchA = pk$A,
                  branchB = pk$B,
                  ratio = pk$`ratio`,
                  adjB = adjB, pct = pct)



    count <- .doCount(data = data, FC = beta,
                      nSam = nSam, mu = mu,
                      size = size, n = n)


    if(is.list(count)) {
        grpDat <- data.frame(group = substr(colnames(count[[1]]), 1, 2))
        lse <- TreeSummarizedExperiment(assays = count,
                                        metadata = list(
                                            FC = beta,
                                            branch = pk,
                                            scenario = scenario),
                                        colData = grpDat,
                                        rowTree = tree,
                                        rowNodeLab = rownames(count[[1]]))
    }

    if(is.matrix(count)) {
        grpDat <- data.frame(group = substr(colnames(count), 1, 2))
        lse <- TreeSummarizedExperiment(assays = list(count),
                                        metadata = list(
                                            FC = beta,
                                            branch = pk,
                                            scenario = scenario),
                                        colData = grpDat,
                                        rowTree = tree,
                                        rowNodeLab = rownames(count))
    }
    #obj <- aggData(x = lse, onRow = TRUE, FUN = sum, message = FALSE)
    return(lse)
}


#' Select branches
#'
#' \code{.pickLoc} selects two branches which meet the criteria specified by
#' the arguments
#'
#' @param tree A phylo object
#' @param data A count table (a matrix or a data frame). It has tree leaves in
#' rows and samples from different conditions in columns.
#' @param from.A The branch node label of branch A. In simulation, we select two
#' branches (A & B) to have differential abundance under different conditions.
#' If from.A is specified, then branch A is fixed. If from.A is NULL, one can
#' find a suitable branch which meets the criteria specified in \code{minTip.A},
#' \code{maxTip.A}, \code{minPr.A} and \code{maxPr.A}.
#' @param minTip.A The minimum number of leaves in branch A
#' @param maxTip.A The maximum number of leaves in branch A
#' @param minPr.A The minimum abundance proportion of leaves in branch A
#' @param maxPr.A The maximum abundance proportion of leaves in branch A
#' @param minTip.B The minimum number of leaves in branch B
#' @param maxTip.B The maximum number of leaves in branch B
#' @param ratio The proportion ratio of branch B to branch A. This value is
#' used to select branches(see \bold{Details}). If there are no branches having
#' exactly this ratio, the pair with the value closest to \code{ratio} would
#' be selected.
#'
#' @return a data frame of one row
#' @author Ruizhu Huang
#' @keywords internal


.pickLoc <- function(tree = NULL, data = NULL,
                     from.A = NULL,
                     minTip.A = 0, maxTip.A = Inf,
                     minTip.B = 0, maxTip.B = Inf,
                     minPr.A = 0, maxPr.A = 1,
                     ratio = 1) {

    # tip proportions estimated from real data
    # rename using the alias of node label
    pars <- parEstimate(obj = data)$pi
    nam1 <- names(pars)
    val1 <- convertNode(tree = tree, node = nam1, message = FALSE)
    nam2 <- convertNode(tree = tree, node = val1, use.alias = TRUE,
                      message = FALSE)
    names(pars) <- nam2
    # df <- data.frame(pr = pars, nodeLab = nam1,
    #                  nodNum = val1, nodeLab_alias = nam2)

    # proportion of internal nodes
    leaf <- setdiff(tree$edge[, 2], tree$edge[, 1])
    leaf <- sort(leaf)
    nodI <- setdiff(tree$edge[, 1], leaf)
    nodI <- sort(nodI)
    desI <- findDescendant(tree = tree, node = nodI,
                   only.leaf = TRUE,
                   self.include = TRUE,
                   use.alias = TRUE)
    desI <- lapply(desI, FUN = function(x) {
        convertNode(tree = tree, node = x,
                  use.alias = TRUE, message = FALSE)})
    names(desI) <- convertNode(tree = tree, node = nodI,
                             use.alias = TRUE, message = FALSE)
    nodP <- mapply(function(x, y) {sum(x[y])},
                   x = list(pars), y = desI)

    # matrix: abundance proprotion & the number of descendant leaves
    lenI <- unlist(lapply(desI, length))
    tt <- cbind(nodP, lenI)
    rownames(tt) <- convertNode(tree = tree,
                              node = nodI,
                              use.alias = TRUE,
                              message = FALSE)

    # return error when the given limits for
    # proportion are outside
    # the range estimated from the real data.
    if (maxPr.A < min(tt[, 1])) {
        stop("maxPr.A is lower than the minimum value of
             node proportion", signif(min(tt[, 1]), 2), "\n")
    }
    if (minPr.A*ratio > max(tt[, 1])) {
        stop("minPr.A*ratio is above the maximum value of
             node proportion; try lower ratio", signif(min(tt[, 1]), 2), "\n")
    }

    # only consider nodes with enough tips and
    # desired proportion level
    if (is.null(from.A)) {
        tt.sel <- tt
    } else {
        # if node numbers, change them to node labels
        if (is.character(from.A)) {
            from.A <- from.A
        } else {
            from.A <- convertNode(tree = tree, node = from.A,
                                use.alias = TRUE,
                                message = FALSE)
        }
        tt.sel <- tt[match(from.A, rownames(tt)), , drop = FALSE]
    }

    st <- tt.sel[tt.sel[, 2] >= minTip.A &
                 tt.sel[, 2] <= maxTip.A &
                 tt.sel[, 1] >= minPr.A &
                 tt.sel[, 1] <= maxPr.A, , drop = FALSE]
    if (nrow(st) == 0) {
        stop("No nodes fullfill the requirements;
             try other values for minTip.A, maxTip.A,
             minPr.A, or maxPr.A")
    }
    st2 <- tt[tt[, 2] >= minTip.B &
              tt[, 2] <= maxTip.B, , drop = FALSE]

    # fold change between any two nodes (large/small)
    mm <- (1/st[, 1]) %o% st2[, 1]
    rownames(mm) <- rownames(st)

    maxM <- max(mm, na.rm = TRUE)
    minM <- min(mm, na.rm = TRUE)

    if (ratio > 1 & ratio > maxM) {
        stop("could not find any two branches which fullfill
             these requirement;
             try lower ratio, lower minTip.A, or higher maxTip.B",
             "\n")
    }

    if (ratio < 1 & ratio < minM) {
        stop("could not find any two branches which fullfill
             these requirement; try higher ratio or lower minTip.B",
             "\n")
    }

    nm <- mm
    nm[] <- vapply(
        seq_len(ncol(mm)),
        FUN = function(x) {
            # each column
            cn <- colnames(mm)
            cx <- cn[x]

            # all rows
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

    du <- cbind.data.frame(
        "A" = convertNode(tree = tree, node = an,
                        use.alias = FALSE, message = FALSE),
        "B" = convertNode(tree = tree, node = bn,
                        use.alias = FALSE, message = FALSE),
        "ratio" = unique(nm[wi]),
        "A_tips" = tt[an, 2],
        "B_tips" = tt[bn, 2],
        "A_prop" = round(tt[an, 1],
                         digits = 4),
        "B_prop" = round(tt[bn, 1],
                         digits = 4),
        stringsAsFactors =  FALSE
    )

    rownames(du) <- NULL
    return(du)

    }

#' Provide the information of two branches
#'
#' \code{.infLoc} is to give information of two branches about the count
#' proportion and the number of leaves
#'
#' @param tree A phylo object
#' @param data A count table (a matrix or a data frame). It has tree leaves in
#' rows and samples from different conditions in  columns.
#' @param from.A,from.B The branch node labels of Branch A, B.
#' @return A data frame of one row
#' @author Ruizhu Huang
#' @keywords internal

.infLoc <- function(tree = NULL, data = NULL,
                    from.A = NULL, from.B = NULL) {

    # tip proportions estimated from real data
    # rename using the alias of node label
    pars <- parEstimate(obj = data)$pi
    nam1 <- names(pars)
    val1 <- convertNode(tree = tree, node = nam1, message = FALSE)
    nam2 <- convertNode(tree = tree, node = val1, use.alias = TRUE,
                      message = FALSE)
    names(pars) <- nam2

    # nodes
    leaf <- setdiff(tree$edge[, 2], tree$edge[, 1])
    leaf <- sort(leaf)
    nodI <- setdiff(tree$edge[, 1], leaf)
    nodI <- sort(nodI)
    nodA <- c(leaf, nodI)

    # find descendants
    desI <- findDescendant(tree = tree, node = nodI,
                   only.leaf = TRUE,
                   self.include = TRUE,
                   use.alias = TRUE)
    desI <- lapply(desI, FUN = function(x) {
        convertNode(tree = tree, node = x,
                  use.alias = TRUE, message = FALSE)})

    nodP <- mapply(function(x, y) {
        sum(x[y])
    }, x = list(pars), y = desI)

    # matrix: abundance proprotion & the number of descendant leaves
    lenI <- unlist(lapply(desI, length))
    tt <- cbind(nodP, lenI)
    rownames(tt) <- convertNode(tree = tree, node = nodI,
                              use.alias = TRUE,
                              message = FALSE)

    # if both branches are given
    labA <- ifelse(is.character(from.A), from.A,
                   convertNode(tree = tree, node = from.A,
                             use.alias = TRUE,
                             message = FALSE))
    labB <- ifelse(is.character(from.B), from.B,
                   convertNode(tree = tree, node = from.B,
                             use.alias = TRUE,
                             message = FALSE))
    rAB <- tt[labB, 1]/tt[labA, 1]
    du <- cbind.data.frame(
        "A" = from.A,
        "B" = from.B,
        "ratio" = rAB,
        "A_tips" = tt[labA, 2],
        "B_tips" = tt[labB, 2],
        "A_prop" = round(tt[labA, 1],
                         digits = 4),
        "B_prop" = round(tt[labB, 1],
                         digits = 4),
        stringsAsFactors =  FALSE)

    rownames(du) <- NULL
    return(du)
}

#' Generate the fold change
#'
#' \code{.doFC} generates fold changes for different scenarios
#'
#' @param tree A phylo object
#' @param data The real data (count table)
#' @param scenario Scenarios (\dQuote{BS}, \dQuote{US}, \dQuote{SS})
#' @param branchA The branch node label of branch A.
#' @param branchB The branh node label of branch B.
#' @param ratio The proportion ratio between \code{branchB} and \code{branchA}
#' (B/A)
#' @param adjB A numeric value between 0 and 1 (only for \code{scenario}
#' is \dQuote{SS}). Default is NULL. If NULL, branch A and branch B swap their
#' proportions. If a numeric value, e.g. 0.1, then branch B decreases to its
#' one tenth proportion and the decrease in branch B is added to branch A.
#' For example, assume there are two experimental conditions (C1 & C2), branch
#' A has 10 and branch B has 40 in C1. If adjB is set to 0.1, then in C2 branch
#' B becomes 4 and branch A 46 so that the total proportion stays the same.
#' @param pct The percentage (in number) of the leaves in branchA that will
#' swap with branchB. The default is 1.

#'
#' @importFrom stats runif
#' @return numeric vector
#' @author Ruizhu Huang
#' @keywords internal

.doFC <- function(tree = NULL, data = NULL, scenario = "BS",
                  branchA = NULL, branchB = NULL,
                  ratio = 1, adjB = NULL, pct = 1) {
    # nodes
    leaf <- setdiff(tree$edge[, 2], tree$edge[, 1])
    leaf <- sort(leaf)
    nodI <- setdiff(tree$edge[, 1], leaf)
    nodI <- sort(nodI)
    nodA <- c(leaf, nodI)

    # beta
    beta <- rep(1, length(leaf))
    names(beta) <- convertNode(tree = tree, node = leaf,
                             use.alias = TRUE)

    ## the label of nodes on branch A
    # leaves
    tip.A <- findDescendant(tree = tree, node = branchA,
                    only.leaf = TRUE, self.include = TRUE,
                    use.alias = TRUE)
    tip.A <- convertNode(tree = tree, node = unlist(tip.A),
                       use.alias = TRUE, message = FALSE)
    # nodes
    nodA.A <- findDescendant(tree = tree, node = branchA,
                     only.leaf = FALSE, self.include = TRUE,
                     use.alias = TRUE)
    nodA.A <- convertNode(tree = tree, node = unlist(nodA.A),
                        use.alias = TRUE, message = FALSE)
    # internal nodes
    nodI.A <- setdiff(nodA.A, tip.A)

    # descendants of internal nodes
    des.IA <- findDescendant(tree = tree, node = nodI.A,
                     only.leaf = TRUE,
                     self.include = TRUE,
                     use.alias = TRUE)
    des.IA <- lapply(des.IA, FUN = function(x) {
        convertNode(tree = tree, node = x,
                  use.alias = TRUE, message = FALSE)
    })

    # tip proportions estimated from real data
    # rename using the alias of node label
    pars <- parEstimate(obj = data)$pi
    nam1 <- names(pars)
    val1 <- convertNode(tree = tree, node = nam1, message = FALSE)
    nam2 <- convertNode(tree = tree, node = val1, use.alias = TRUE,
                      message = FALSE)
    names(pars) <- nam2

    # swap proportion of two branches: tips in the same branch
    # have the same fold change
    if (scenario == "BS") {
        # leaves on branch B
        tip.B <- findDescendant(tree = tree, node = branchB,
                        only.leaf = TRUE, self.include = TRUE,
                        use.alias = TRUE)
        tip.B <- convertNode(tree = tree, node = unlist(tip.B),
                           use.alias = TRUE, message = FALSE)
        beta[tip.A] <- ratio
        beta[tip.B] <- 1/ratio
    }

    # swap proportion of two branches: tips in the same branch
    # have different fold changes but same direction (either
    # increase or decrease)
    if (scenario == "US") {
        tip.B <- findDescendant(tree = tree, node = branchB,
                        only.leaf = TRUE, self.include = TRUE,
                        use.alias = TRUE)
        tip.B <- convertNode(tree = tree, node = unlist(tip.B),
                           use.alias = TRUE, message = FALSE)

        # proportion on two branches
        propA <- sum(pars[tip.A])
        propB <- sum(pars[tip.B])

        # make sure the ratio is above zero
        if (propB < propA) {
            rA <- 1-propB/propA
            rB <- 0
        } else {
            rA <- 0
            rB <- 1-propA/propB
        }

        # swap proportions on two branches and randomly assign a fold change
        # value to the the leaves on a branch (log fold change in the same
        # branch should have the same sign)

        # a1 <- runif(length(tip.A))
        # sa <- sum(a1*pars[tip.A])
        # a2 <- (propB - propA)/sa
        # a3 <- a1 * a2 + 1
        # beta[tip.A] <- a3
        #
        # b1 <- runif(length(tip.B), probA/probB, 1)
        # sb <- sum(b1 * pars[tip.B])
        # b2 <- (propA - propB)/sb
        # b3 <- b1 * b2 + 1
        #

        a1 <- runif(length(tip.A), rA, 1)
        sa <- sum(a1*pars[tip.A])
        a2 <- (propB - propA)/sa
        a3 <- a1 * a2 + 1
        beta[tip.A] <- a3

        # b1 <- propB/propA
        # b2 <- runif(length(tip.B), 1 - b1)
        # b3 <- b1 + b2

        # b3 <- propA/propB
        b1 <- runif(length(tip.B), rB, 1)
        sb <- sum(b1 * pars[tip.B])
        b2 <- (propA - propB)/sb
        b3 <- b1 * b2 + 1
        beta[tip.B] <- b3
    }

    # distribute signal randomly in one branch and evenly in
    # another branch
    # tip proportions estimated from real data

    if (scenario == "SS") {
        iter <- 1
        while (iter <= 200) {
            # select only some tips
            selA <- sample(tip.A, ceiling(length(tip.A)*pct))
            # to make the selected tips disperse evenly in the branch
            subA <- lapply(des.IA, FUN = function(x) {
                ix <- intersect(x, selA)
                length(ix)/length(x)
            })
            subA <- unlist(subA)
            sumA <- sum(pars[selA])

            # the abundance proportion of the selected tips
            # are roughly equal to its number proportion in the branch
            # avoid (select all low or high abundance tips in the branch)

            spr <- sumA/sum(pars[tip.A])
            ind.pr <- spr <= (pct + 0.05) & spr >= (pct - 0.05)

            if(all(subA <= 0.6) & ind.pr){break}
            iter <- iter+1
        }

        if (length(branchB) == 0) {
            stop("No suitable branches.
                 Try another branchA or another max of ratio... \n")
        }
        tip.B <- findDescendant(tree = tree, node = branchB,
                        only.leaf = TRUE, self.include = TRUE,
                        use.alias = TRUE)
        tip.B <- convertNode(tree = tree, node = unlist(tip.B),
                           use.alias = TRUE, message = FALSE)
        sumB <- sum(pars[tip.B])

        if(is.null(adjB)){
            beta[selA] <- sumB/sumA
            beta[tip.B] <- sumA/sumB
        }else{
            if(!is.numeric(adjB)){
                stop("adjB should be numeric")
            }
            beta[tip.B] <- adjB
            beta[selA] <- (sumB*(1-adjB)+sumA)/sumA
        }
        }

    # rename beta with the node label instead of the alias of node label
    names(beta) <- convertNode(tree = tree, node = leaf,
                             use.alias = FALSE)
    return(beta)
}


#' generate a count table
#'
#' \code{.doCount} generates a count table given some available information.
#' Here, the information includes the fold change between two conditions,
#' the parameters of Negative Binomial distritubion (which the sample library
#' size follows), a data table from real data (to estimate the proportion of
#' each entity in the sample), and the number of samples in two different
#' groups or conditions.
#'
#' @param data A matrix or data frame, representing the count table from real
#' data
#' @param FC A numeric vector, representing the fold changes
#' @param nSam A vector of length two containing the number of samples in two
#' different groups or conditions.
#' @param mu,size The parameters of the Negative Binomial distribution (see mu
#' and size in \code{\link[stats:NegBinomial]{rnbinom}}.)
#' @param n A numeric value, representing the number of count tables generated.
#' If above one, the output is a list of matrices.
#'
#' @importFrom dirmult rdirichlet
#' @importFrom stats rmultinom rnbinom
#' @return a matrix or a list of matrices
#' @author Ruizhu Huang
#' @keywords internal

.doCount <- function(data, FC, nSam, mu,
                     size, n) {
    # parameters
    pars <- parEstimate(obj = data)
    theta <- pars$theta
    gplus <- (1 - theta) / theta

    # tip proportion
    pr <- pars$pi
    p.c1 <- pr
    p.c2 <- pr * FC[names(p.c1)]

    # parameters for dirichlet distribution
    g.c1 <- p.c1 * gplus
    g.c2 <- p.c2 * gplus

    resList <- lapply(seq_len(n), FUN = function(j) {
        # condition 1
        n1 <- nSam[1]
        Mp.c1 <- matrix(0, nrow = n1, ncol = length(g.c1))
        rownames(Mp.c1) <- paste("C1_", seq_len(n1), sep = "")
        colnames(Mp.c1) <- names(p.c1)
        Mobs.c1 <- Mp.c1
        if (length(mu) & !length(size) ) {
            nSeq1 <- sample(x = mu, size = n1, replace = TRUE)
        } else {
            nSeq1 <- rnbinom(n = n1, mu = mu, size = size)
        }
        for (i in seq_len(n1)) {
            Mp.c1[i, ] <- rdirichlet(1, g.c1)[1, ]
            Mobs.c1[i, ] <- rmultinom(1, nSeq1[i], prob = Mp.c1[i, ])[, 1]

        }

        # condition 2
        n2 <- nSam[2]
        Mp.c2 <- matrix(0, nrow = n2, ncol = length(g.c2))
        rownames(Mp.c2) <- paste("C2_", seq_len(n2), sep = "")
        colnames(Mp.c2) <- names(p.c2)
        Mobs.c2 <- Mp.c2
        if (length(mu) & !length(size) ) {
            nSeq2 <- sample(x = mu, size = n2, replace = TRUE)
        } else {
            nSeq2 <- rnbinom(n = n2, mu = mu, size = size)
        }

        for (i in seq_len(n2)) {
            Mp.c2[i, ] <- rdirichlet(1, g.c2)[1, ]
            Mobs.c2[i, ] <- rmultinom(1, nSeq2[i], prob = Mp.c2[i, ])[, 1]
        }

        cb <- t(rbind(Mobs.c1, Mobs.c2))

        return(cb)}
    )

    if (n == 1) {
        count <-  do.call(rbind, resList)
    } else {
        count <- resList
    }

    return(count)
}

