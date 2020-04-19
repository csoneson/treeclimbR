#' Simulate count tables of entities (multinomial distribution)
#'
#' \code{simMult} simulates count tables of entities from multinomial
#' distribution based on the provided proportions. The entities are in the table
#' rows and samples are in the columns. These entities have their corresponding
#' nodes on the tree. More details about the simulated patterns could be found
#' in the vignette via \code{browseVignettes("treeAGG")}.
#' @param pr A named numeric values to provide the proportions of entities in a
#'   sample. It should be named using the leaf labels.
#' @param libSize A numeric values to provide the total counts of samples. The
#'   length should be one or equal to the sum of \code{nSam}.
#' @param tree A phylo object. Only use when \code{obj} is NULL.
#' @param scenario \dQuote{S1}, \dQuote{S2}, or \dQuote{S3} (see
#'   \bold{Details}). Default is \dQuote{S1}.
#' @param from.A,from.B The branch node labels of branches A and B for which the
#'   signal is swapped. Default, both are NULL. In simulation, we select two
#'   branches (A & B) to have differential abundance under different conditions.
#'   One could specify these two branches or let \code{simData} choose. (Note: If
#'   \code{from.A} is NULL, \code{from.B} is set to NULL).
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
#'   is \dQuote{S3}). Default is NULL. If NULL, branch A and the selected part
#'   of branch B swap their proportions. If a numeric value, e.g. 0.1, then the
#'   selected part of branch B decreases to its one tenth proportion and the
#'   decrease in branch B is added to branch A. For example, assume there are
#'   two experimental conditions (C1 & C2), branch A has 10 and branch B has 40
#'   in C1. If adjB is set to 0.1, then in C2 branch B becomes 4 and branch A 46
#'   so that the total proportion stays the same.
#' @param pct The percentage of leaves in branch B that have differential
#'   abundance under different conditions (only for scenario \dQuote{S3})
#' @param nSam A numeric vector of length 2, containing the sample size for two
#'   different conditions
#' @param n A numeric value to specify how many count tables would be generated
#'   with the same settings. Default is one and one count table would be
#'   obtained at the end. If above one, the output is a list of matrices (count
#'   tables). This is useful, when one needs multiple simulations.
#' @param message A logical value, TRUE or FALSE. Default is FALSE. If TRUE, the
#'   progress message is printed out.
#'
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom methods is
#' @importFrom SummarizedExperiment assays
#' @export
#'
#' @return a TreeSummarizedExperiment object.
#'  \itemize{
#'  \item \strong{assays} a list of count table or a count table. Entities on
#'  the row and samples in the column. Each row could be mapped to a leaf node
#'  of the tree.
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
#' @details \code{simMult} simulates a count table for entities which are
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
#'   \dQuote{S1}, \dQuote{S2}, and \dQuote{S3} (specified via \code{scenario}).
#'   Our vignette provides figures to explain these three scenarios (try
#'   \code{browseVignettes("treeAGG")}). \itemize{ \item S1: two branches are
#'   selected to swap their proportions, and leaves on the same branch have the
#'   same fold change. \item S2: two branches are selected to swap their
#'   proportions. Leaves in the same branch have different fold changes but same
#'   direction (either increase or decrease). \item S3: two branches are
#'   selected. One branch has its proportion swapped with the proportion of some
#'   leaves from the other branch.}
#'
#' @author Ruizhu Huang
#'
#' @examples
#' 
#' library(TreeSummarizedExperiment)
#' library(ComplexHeatmap)
#' library(ggtree)
#' 
#' data("tinyTree")
#' 
#' set.seed(1)
#' n <- 10
#' prob <- runif(n = 10)
#' prob <- prob/sum(prob)
#' names(prob) <- tinyTree$tip.label
#' 
#' tse <- simMult(pr = prob, libSize = 500, tree = tinyTree)
#' viewSim(tse) + geom_tiplab()
#' 
#' # heatmap
#' 
#' count <- assays(tse)[[1]]
#' rownames(count) <- rowLinks(tse)$nodeLab
#' 


simMult <- function(pr, libSize, tree, scenario = "S1",
                    from.A = NULL, from.B = NULL,
                    minTip.A = 0, maxTip.A = Inf,
                    minTip.B = 0, maxTip.B = Inf,
                    minPr.A = 0, maxPr.A = 1,
                    ratio = 2, adjB = NULL,
                    pct = 0.6, nSam = c(50, 50),
                    n = 1, message = FALSE) {
    
    if (length(libSize) == 1) {
        libSize <- rep(libSize, 2)
    }
    
    # select branches
    data = list(pi = pr, theta = NULL)
    
    leaf <- unique(setdiff(tree$edge[, 2], tree$edge[, 1]))
    leaf <- sort(leaf)
    if (message) {
        message("prepare branches ...")
    }
    if (scenario == "S0") {
        
        if (message) {
            message("calculate the fold change ...")
        }
        
        beta <- rep(1, length(leaf))
        names(beta) <- transNode(tree = tree, node = leaf, use.alias = FALSE)
        pk <- NULL
    } else {
       # scenario S1, S2, S3
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
        
        if (message) {
            message("calculate the fold change ...")
        }
        # generate the fold change
        beta <- .doFC(tree = tree, data = data,
                      scenario = scenario,
                      branchA = pk$A,
                      branchB = pk$B,
                      ratio = pk$`ratio`,
                      adjB = adjB, pct = pct) 
    }
    
    
    
    
    # generate counts
    if (message) {
        message("generate counts ...")
    }
    if (!all.equal(names(pr), tree$tip.label)) {
        stop("pr should be named with the tip labels of the tree.")
    }
    p.c1 <- pr
    p.c2 <- pr * beta[names(p.c1)]
    resList <- lapply(seq_len(n), FUN = function(j) {
        Mp.c1 <- lapply(seq_len(nSam[1]), FUN = function(x) {
            lib <- sample(x = libSize, size = 1)
            rmultinom(n = 1, size = lib, prob = p.c1)
            })
        Mp.c1 <- do.call(cbind, Mp.c1)
        Mp.c2 <- lapply(seq_len(nSam[2]), FUN = function(x) {
            lib <- sample(x = libSize, size = 1)
            rmultinom(n = 1, size = lib, prob = p.c2)})
        Mp.c2 <- do.call(cbind, Mp.c2)
        # Mp.c2 <- rmultinom(n = nSam[2], size = libSize, prob = p.c2)
        Mp <- cbind(Mp.c1, Mp.c2)
        rownames(Mp) <- names(pr)
        colnames(Mp) <- c(paste0("C1_", seq_len(nSam[1])), 
                          paste0("C2_", seq_len(nSam[2])))
        return(Mp)
    })
    if (n == 1) {
        count <-  do.call(rbind, resList)
    } else {
        count <- resList
    }
    
    # output as a TreeSummarizedExperiment object
    if (message) {
        message("output the TreeSummarizedExperiment object ...")
    }
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
    return(lse)
    
}



