% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tpr.R
\name{tpr}
\alias{tpr}
\title{Calculate true positive rate (TPR) on a tree structure}
\usage{
tpr(tree, truth, found, only.leaf = TRUE)
}
\arguments{
\item{tree}{A phylo object}

\item{truth}{Nodes that have signals (eg. differentally abundant at different
experimental conditions.). \strong{Note:} when TPR is required at leaf
level (\code{only.leaf = TRUE}), the descendant leaves of the given nodes
will be found out and the TPR is calculated on the leaf level;}

\item{found}{Nodes that have been found to have signal (eg. differentally
abundant at different experimental conditions). \strong{Note:} when TPR is
required at leaf level (\code{only.leaf = TRUE}), the descendant leaves of
the given nodes will be found out and the TPR is calculated on the leaf
level;}

\item{only.leaf}{A logical value, TRUE or FALSE. If TRUE, true positive rate
is calculated at the leaf (tip) level; otherwise it is calculated at node
level. The default is TRUE}
}
\value{
A true positive rate
}
\description{
\code{tpr} calculates the true positive rate (TPR) on a tree structure at
leaf or node level.
}
\examples{
library(ggtree)
library(TreeSummarizedExperiment)
data("tinyTree")
ggtree(tinyTree) + 
   geom_text2(aes(label = node)) + 
   geom_hilight(node = 16, fill = "orange", alpha = 0.3) +
   geom_hilight(node = 13, fill = "blue", alpha = 0.3)

# Truth: two branches have differential
# abundance under different conditions.

# Found: branches with node 17 and node 14
 # TPR at the tip level
tpr1 <- tpr(tree = tinyTree, truth = c(16, 13),
            found = c(17, 14), only.leaf = TRUE)
 # TPR at the node level
tpr2 <- tpr(tree = tinyTree, truth = c(16, 13),
            found = c(15, 14), only.leaf = FALSE)


}
\author{
Ruizhu Huang
}