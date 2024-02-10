library(TreeSummarizedExperiment)
library(ggtree)
library(ggplot2)
library(scales)

data(tinyTree)

## Simulate some data to display
p1 <- c(rep(0.1/3, 3), rep(0.4/2, 2), rep(0.1, 5))
p2 <- c(rep(0.4/3, 3), rep(0.1/2, 2), rep(0.1, 5))
set.seed(1)
ct0 <- cbind(rmultinom(n = 5, size = 50, prob = p1),
             rmultinom(n = 5, size = 50, prob = p2))
colnames(ct0) <- paste("S", seq_len(10), sep = "")
rownames(ct0) <- convertNode(tree = tinyTree, node = seq_len(10))
oo <- sample(seq_len(10))
ct0 <- ct0[, oo] ## count matrix on the leaf level

## Create TSE
tse <- TreeSummarizedExperiment(
    assays = list(counts = ct0),
    rowTree = tinyTree,
    colData = DataFrame(group = rep(c("A", "B"), each = 5)[oo]))

saveRDS(tse, file = "inst/extdata/tinytree_counts.rds")

