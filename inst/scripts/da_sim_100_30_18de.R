suppressPackageStartupMessages({
    library(dplyr)
    library(treeclimbR)
    library(TreeSummarizedExperiment)
    library(ape)
})

## Generate a tree with 100 leaves
set.seed(2020)
n <- 100
tr <- rtree(n)
tr

## Generate a random probability vector for leaves
p <- runif(n = n, 0, 1)
p <- p/sum(p)
names(p) <- tr$tip.label

## Introduce differential signal
## First select an internal node - all leaves with signal will be in the
## subtree defined by this node (to simplify visualization later)
## The node must have between 20 and 30 leaves
df <- selNode(pr = p, tree = tr, all = TRUE)
nd <- df |>
    filter(numTip > 20 & numTip < 30) |>
    top_n(1) |>
    select(nodeNum) |>
    unlist()
nd

## Next, randomly select 18 leaves from the branch
m <- 18
lf <- unlist(findDescendant(tree = tr, node = nd, only.leaf = TRUE))
lfs <- sample(lf, size = m, replace = FALSE)
lfs <- convertNode(tree = tr, node = lfs)
lfs

## Simulate samples in two groups
nSam <- c(15, 15)
gr <- rep(LETTERS[seq_len(2)], nSam)

## Fold change
fc <- 2

## Simulate counts
count <- rmultinom(n = sum(nSam), size = 500, prob = p)
rownames(count) <- names(p)

## Multiply counts of selected leaves with fc in the first group
count[lfs, seq_len(nSam[1])] <- count[lfs, nSam[1] + seq_len(nSam[1])] * fc
colnames(count) <- paste(gr, seq_len(sum(nSam)), sep = "_")

## Build TSE
lse <- TreeSummarizedExperiment(assays = list(counts = count),
                                colData = data.frame(group = gr),
                                rowTree = tr)

## Add information about signal leaves
rowData(lse)$Signal <- FALSE
rowData(lse)$Signal[rownames(lse) %in% lfs] <- TRUE

## Add information about the node number under which the signal leaves
## were generated
metadata(lse)$parentNodeForSignal <- nd

saveRDS(lse, file = "inst/extdata/da_sim_100_30_18de.rds")
