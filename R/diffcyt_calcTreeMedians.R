#' Calculate cluster medians on the tree structure \code{calcTreeMedians}
#' calculates cluster medians (median expression for each cluster-sample-marker
#' combination). This is a tree version of \code{\link[diffcyt]{calcMedians}}.
#' The clusters used in \code{\link[diffcyt]{calcMedians}} are clusters on the
#' leaf level of the tree here. More details about the data could be found in
#' \code{\link[diffcyt]{calcMedians}}.
#' 
#' @param d_se Data object from previous steps, in
#'   \code{\link[TreeSummarizedExperiment]{TreeSummarizedExperiment}} format,
#'   containing cluster labels as a column in the row meta-data (from
#'   \code{\link[diffcyt]{generateClusters}}). Column meta-data is assumed to
#'   contain a factor marker_class.
#' @param tree A phylo object from \code{\link{buildTree}}
#' @param message A logic value, TRUE or FALSE
#' @return A TreeSummarizedExperiment object.  Calculate median marker
#'   expression for each cluster (each node of the tree) and sample (i.e.
#'   medians for each cluster-sample-marker combination).
#'
#'   The data object is assumed to contain a factor \code{marker_class} in the
#'   column meta-data (see \code{\link[diffcyt]{prepareData}}), which indicates the
#'   protein marker class for each column of data (\code{"type"},
#'   \code{"state"}, or \code{"none"}).
#'
#'   The cluster medians are required for testing for differential states within
#'   cell populations, and for plotting purposes.
#'
#'   Variables \code{id_type_markers} and \code{id_state_markers} are saved in
#'   the \code{metadata} slot of the output object. These can be used to
#'   identify the 'cell type' and 'cell state' markers in the list of
#'   \code{assays} in the output
#'   \code{\link[TreeSummarizedExperiment]{TreeSummarizedExperiment}} object,
#'   which is useful in later steps of the 'diffcyt' pipeline.
#'
#'   Results are returned as a new
#'   \code{\link[TreeSummarizedExperiment]{TreeSummarizedExperiment}} object,
#'   where rows = clusters, columns = samples, sheets (\code{assays} slot) =
#'   markers. Note that there is a separate table of values (\code{assay}) for
#'   each marker. The \code{metadata} slot also contains variables
#'   \code{id_type_markers} and \code{id_state_markers}, which can be used to
#'   identify the sets of cell type and cell state markers in the list of
#'   \code{assays}.
#' @import TreeSummarizedExperiment
#' @importFrom dplyr select distinct
#' @importFrom methods is
#' @importFrom stats median
#' @export
#' @return A TreeSummarizedExperiment object
#' @examples 
#' # For a complete workflow example demonstrating each step, please see the
#' # vignette of 'diffcyt'
#' \dontrun{
#' library(diffcyt)
#' 
#' # Function to create random data (one sample)
#' d_random <- function(n = 20000, mean = 0, sd = 1, ncol = 20, cofactor = 5) {
#'     d <- sinh(matrix(rnorm(n, mean, sd), ncol = ncol)) * cofactor
#'     colnames(d) <- paste0("marker", sprintf("%02d", 1:ncol))
#'     d
#' }
#' 
#' # Create random data (without differential signal)
#' set.seed(123)
#' d_input <- list(
#'     sample1 = d_random(), 
#'     sample2 = d_random(), 
#'     sample3 = d_random(), 
#'     sample4 = d_random()
#' )
#' 
#' experiment_info <- data.frame(
#'     sample_id = factor(paste0("sample", 1:4)), 
#'     group_id = factor(c("group1", "group1", "group2", "group2")), 
#'     stringsAsFactors = FALSE
#' )
#' 
#' marker_info <- data.frame(
#'     channel_name = paste0("channel", sprintf("%03d", 1:20)), 
#'     marker_name = paste0("marker", sprintf("%02d", 1:20)), 
#'     marker_class = factor(c(rep("type", 10), rep("state", 10)), 
#'                           levels = c("type", "state", "none")), 
#'     stringsAsFactors = FALSE
#' )
#' 
#' # Prepare data
#' d_se <- prepareData(d_input, experiment_info, marker_info)
#' 
#' # Transform data
#' d_se <- transformData(d_se)
#' 
#' # Generate clusters
#' d_se <- generateClusters(d_se)
#' 
#' # build a tree
#' tr <- buildTree(d_se)
#' 
#' # calculate medians on nodes of a tree
#' d_medians_tree <- calcTreeMedians(d_se = d_se, tree = tr)
#' }

calcTreeMedians <- function(d_se, tree, message = FALSE) {
    
    if (!("cluster_id" %in% (colnames(rowData(d_se))))) {
        stop("Data object does not contain cluster labels. 
             Run 'diffcyt::generateClusters' to generate cluster labels.")
    }
    
    if (!is(tree, "phylo")) {
        stop("tree is not a phylo object.
             Run 'buildTree(d_se)' to generate the tree")
    }
    
    ## create the tse object
    rlab <- as.character(rowData(d_se)$cluster_id)
    d_lse <- TreeSummarizedExperiment(assays = assays(d_se),
                                      rowData = rowData(d_se),
                                      rowTree = tree, 
                                      rowNodeLab = rlab,
                                      colData = colData(d_se),
                                      metadata = metadata(d_se))
    d_lse <- d_lse[, colData(d_lse)$marker_class %in% c("type", "state")]
    
    ## split by patient_id
    rd <- rowData(d_lse)
    sid <- unique(rd$sample_id)
    nodes <- sort(unique(as.vector(tree$edge)))
    labs <- convertNode(tree = tree, 
                        node = nodes, use.alias = TRUE)
    lc <- lapply(seq_along(sid), FUN = function(i) {
        if (message) {
            message("Working on ", i, " out of ", length(sid), " samples.")
        }
        x <- sid[i]
        sel <- rd$sample_id == x
        xx <- d_lse[sel, ]
        ax <- aggTSE(x = xx, rowLevel = nodes, 
                       rowFun = function(x){
                           median(x, na.rm = TRUE)
                       }, message = message)
        
        # counts
        cx <- assays(ax)[[1]]
        
        # for missing nodes
        nam <- setdiff(labs, rownames(cx))
        if (length(nam)) {
            mx <- matrix(NA, nrow = length(nam), ncol = ncol(cx))
            rownames(mx) <- nam
            colnames(mx) <- colnames(cx)
        } else {
            mx <- NULL
        }
        
        cx <- rbind(cx, mx)[labs, ]
        return(cx)
    })
    
    ## assays : reform data to split by markers
    rc <- do.call(rbind, lc)
    lcm <- lapply(seq_len(ncol(rc)), FUN = function(x) {
        xx <- matrix(rc[, x], ncol = length(sid), byrow = FALSE)
        colnames(xx) <- sid
        rownames(xx) <- labs
        return(xx)
    })
    names(lcm) <- colnames(d_lse)
    
    # rowdata
    rd <- cbind.data.frame(cluster_id = factor(labs, levels = labs))
    rownames(rd) <- labs
    
    # column data
    group_id <- sample_id <- NULL
    cd <- rowData(d_se) %>%
        data.frame() %>%
        select(group_id, sample_id) %>%
        distinct() 
    cd <- cd[match(cd$sample_id, sid), ]
    
    # metadata
    tM <- colData(d_lse)$marker_name[colData(d_lse)$marker_class == "type"]
    sM <- colData(d_lse)$marker_name[colData(d_lse)$marker_class == "state"]
    md <- list(id_type_markers = colnames(d_lse) %in% tM,
               id_state_markers = colnames(d_lse) %in% sM)
    
    # row node label
    rnl <- convertNode(tree = rowTree(d_lse), node = nodes, use.alias = TRUE)
    # output
    out <- TreeSummarizedExperiment(assays = lcm,
                                    rowData  = rd,
                                    colData = cd,
                                    metadata = md,
                                    rowTree = rowTree(d_lse),
                                    rowNodeLab = rnl)
    return(out)
}
