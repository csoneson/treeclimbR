% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nodeResult.R
\name{nodeResult}
\alias{nodeResult}
\title{Extract table with node-level DA/DS results}
\usage{
nodeResult(
  object,
  n = 10,
  type = c("DA", "DS"),
  adjust_method = "BH",
  sort_by = "PValue",
  p_value = 1
)
}
\arguments{
\item{object}{The output from \code{\link{runDA}} or \code{\link{runDS}}.}

\item{n}{An integer indicating the maximum number of entities to return.}

\item{type}{Either "DA" (for \strong{object} from \code{\link{runDA}}) or
"DS" (for \strong{object} from \code{\link{runDS}}).}

\item{adjust_method}{A character string specifying the method used to adjust
p-values for multiple testing. See \code{\link[stats]{p.adjust}} for
possible values.}

\item{sort_by}{A character string specifying the sorting method. This will
be passed to \code{\link[edgeR]{topTags}}.
Possibilities are "PValue" for p-value, "logFC" for absolute log-fold
change or "none" for no sorting.}

\item{p_value}{A numeric cutoff value for adjusted p-values. This will
be passed to \code{\link[edgeR]{topTags}}. Only entities
with adjusted p-values equal or lower than specified are returned.}
}
\value{
A data frame with results for all nodes passing the imposed
    thresholds. The columns \strong{logFC}, \strong{logCPM},
    \strong{PValue}, \strong{FDR}, \strong{F} (or \strong{LR}) are from (the
    output table of) \code{\link[edgeR]{topTags}}. The \strong{node} column
    stores the node number for each entity. Note: \strong{FDR} is corrected
    over all features and nodes when the specified \code{type = "DS"}.
}
\description{
Extract a table with the top-ranked nodes from a DA/DS analysis output
(generated by \code{\link{runDA}} or \code{\link{runDS}}).
}
\examples{
suppressPackageStartupMessages({
    library(TreeSummarizedExperiment)
})

lse <- readRDS(system.file("extdata", "da_sim_100_30_18de.rds",
                           package = "treeclimbR"))
tse <- aggTSE(lse, rowLevel = showNode(tree = rowTree(lse),
                                       only.leaf = FALSE))
dd <- model.matrix( ~ group, data = colData(tse))
out <- runDA(TSE = tse, feature_on_row = TRUE,
             assay = "counts", option = "glmQL",
             design = dd, contrast = NULL,
             normalize = TRUE)

## Top 10 nodes with DA
nodeResult(out, n = 10)

}
\author{
Ruizhu Huang, Charlotte Soneson
}
