% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/topNodes.R
\name{topNodes}
\alias{topNodes}
\title{Generate a table of top-ranked entities (nodes)}
\usage{
topNodes(
  object,
  n = 10,
  sort_by = NULL,
  sort_decreasing = FALSE,
  sort_by_absolute = FALSE,
  p_value = 1
)
}
\arguments{
\item{object}{An output object from \link{evalCand}.}

\item{n}{An integer, the maximum number of entities to return.}

\item{sort_by}{A character string specifying the column of
\code{object$output} to sort by. Set to \code{NULL} to return without
sorting.}

\item{sort_decreasing}{A logical value indicating whether to sort by
decreasing value of the \code{sort_by} column.}

\item{sort_by_absolute}{A logical value indicating whether to take the
absolute value of the \code{sort_by} column before sorting.}

\item{p_value}{A numeric cutoff value for adjusted p-values. Only entities
with adjusted p-values equal or lower than specified are returned.}
}
\value{
A \code{data.frame} with test results. The \strong{node}
column stores the node number for each entity.
}
\description{
Generate a table of top-ranked nodes from the optimal resolution candidate
of entities on a tree.
}
\examples{
suppressPackageStartupMessages({
    library(TreeSummarizedExperiment)
    library(ggtree)
})

data(tinyTree)
ggtree(tinyTree, branch.length = "none") +
   geom_text2(aes(label = node)) +
   geom_hilight(node = 13, fill = "blue", alpha = 0.3) +
   geom_hilight(node = 18, fill = "orange", alpha = 0.3)
set.seed(1)
pv <- runif(19, 0, 1)
pv[c(seq_len(5), 13, 14, 18)] <- runif(8, 0, 0.001)

fc <- sample(c(-1, 1), 19, replace = TRUE)
fc[c(seq_len(3), 13, 14)] <- 1
fc[c(4, 5, 18)] <- -1
df <- data.frame(node = seq_len(19),
                 pvalue = pv,
                 logFoldChange = fc)
ll <- getCand(tree = tinyTree, score_data = df,
               node_column = "node",
               p_column = "pvalue",
               sign_column = "logFoldChange")
cc <- evalCand(tree = tinyTree, levels = ll$candidate_list,
               score_data = df, node_column = "node",
               p_column = "pvalue", sign_column = "logFoldChange",
               limit_rej = 0.05)

## Unsorted result table
topNodes(cc)

## Sort by p-value in increasing order
topNodes(cc, sort_by = "pvalue")

}
\author{
Ruizhu Huang, Charlotte Soneson
}
