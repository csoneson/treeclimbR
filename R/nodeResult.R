#' Table of the top differentially abundant entities
#'
#' Extracts the most differentially abundant entities from a test object, ranked
#' either by p-value or by absolute log-fold-change.
#'
#' @param object The output from \link{runDA} or \link{runDS}.
#' @param type "DA" (\strong{object} from \link{runDA}) or "DS" (\strong{object}
#'   from \link{runDS})
#' @param n A integer, maximum number of entities to return.
#' @param adjust_method A character string specifying the method used to adjust
#'   p-values for multiple testing. See \code{\link[stats]{p.adjust}} for
#'   possible values.
#' @param sort_by A character string specifying the sort method. Possibilities
#'   are "PValue" for p-value, "logFC" for absolute log-fold change or "none"
#'   for no sorting.
#' @param p_value A numeric cutoff value for adjusted p-values. Only entities
#'   with adjusted p-values equal or lower than specified are returned.
#'
#' @importFrom edgeR topTags
#' @importFrom dplyr arrange slice filter '%>%'
#' @importFrom TreeSummarizedExperiment convertNode
#' @importFrom stats p.adjust
#' @importFrom data.table rbindlist
#' @export
#' @return A data frame. Columns including \strong{logFC}, \strong{logCPM},
#'   \strong{PValue}, \strong{FDR}, \strong{F} (or \strong{LR}) are from (the
#'   output table of) \code{\link[edgeR]{topTags}}. The \strong{node} column
#'   stores the node number for each entities. Note: \strong{FDR} is corrected
#'   over all features and nodes when the specified \code{type = "DS"}.
#' @examples
#' library(TreeSummarizedExperiment)
#' library(treeclimbR)
#' set.seed(1)
#' count <- matrix(rnbinom(300,size=1,mu=10),nrow=10)
#' colnames(count) <- paste(rep(LETTERS[1:3], each = 10),
#'                          rep(1:10, 3), sep = "_")
#' rownames(count) <- tinyTree$tip.label
#' count[1, ] <- 0
#' rowInf <- DataFrame(var1 = sample(letters[1:3], 10, replace = TRUE),
#'                     var2 = sample(c(TRUE, FALSE), 10, replace = TRUE))
#' colInf <- DataFrame(gg = factor(sample(1:3, 30, replace = TRUE)),
#'                     group = rep(LETTERS[1:3], each = 10))
#' lse <- TreeSummarizedExperiment(assays = list(count),
#'                                 rowData = rowInf,
#'                                 colData = colInf,
#'                                 rowTree = tinyTree)
#' nodes <- showNode(tree = tinyTree, only.leaf = FALSE)
#' tse <- aggTSE(x = lse, rowLevel = nodes)
#'
#' dd <- model.matrix( ~ group, data = colInf)
#' out <- runDA(TSE = tse, feature_on_row = TRUE,
#'              assay = 1, option = "glmQL",
#'              design = dd, contrast = NULL,
#'              normalize = TRUE,
#'              group_column = "group")
#'
#' topOut <- nodeResult(out, n = 10)
#'

nodeResult <- function(object, n = 10,
                     type = c("DA", "DS"),
                     adjust_method = "BH",
                     sort_by = "PValue",
                     p_value = 1) {
    type <- match.arg(type)
    if (type == "DA") {
        tt <- topTags(object = object$edgeR_results, n = n,
                      adjust.method = adjust_method,
                      sort.by = sort_by,
                      p.value = p_value)$table
        # add nodes
        nod <- convertNode(tree = object$tree, node = rownames(tt))
        ct <- cbind(node = nod, tt)
    }

    if (type == "DS") {
        res <- object$edgeR_results
        tt <- lapply(seq_along(res), FUN = function(x) {
            xx <- topTags(object = res[[x]], n = Inf,
                          adjust.method = adjust_method,
                          sort.by = sort_by,
                          p.value = p_value)$table
            # add nodes
            nod <- convertNode(tree = object$tree,
                             node = names(res)[x])
            cx <- cbind(xx, node = nod, feature = rownames(xx),
                        stringsAsFactors = FALSE)
            return(cx)
        })
        ct <- rbindlist(tt)

        # correct FDR
        ct$FDR <- p.adjust(p = ct$PValue, method = adjust_method)
        n <- min(nrow(ct), n)
        ct <- ct %>%
            arrange(!!as.symbol(sort_by)) %>%
            filter(FDR <= p_value) %>%
            slice(seq_len(n))

    }

    return(ct)
}
