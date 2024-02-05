library(TreeSummarizedExperiment)
library(ggtree)

data(tinyTree)

set.seed(1L)

## Simulate data
ncell <- 500
ngene <- 20
dat <- matrix(rpois(ngene*ncell, lambda = 20), nrow = ngene)

## Cells belong to different clusters
clus <- sample(tinyTree$tip.label, 500, replace = TRUE)
sid <- sample(seq_len(4), 500, replace = TRUE)
gid <- ifelse(sid %in% seq_len(2), "A", "B")
colD <- data.frame(cluster_id = clus,
                   sample_id = sid,
                   group = gid,
                   stringsAsFactors = FALSE)

## Induce signal
## Genes that are DE between groups A and B for all cell types
allde_up <- c(1, 3)
dat[allde_up, colD$group == "B"] <- 3 * dat[allde_up, colD$group == "B"]
allde_dn <- c(5, 6)
dat[allde_dn, colD$group == "A"] <- 3 * dat[allde_dn, colD$group == "A"]

## Genes that are DE for cell types t6, t7 and t2
somede_up <- c(2)
dat[somede_up, colD$group == "B" & colD$cluster_id %in% c("t2", "t6", "t7")] <-
    3 * dat[somede_up, colD$group == "B" & colD$cluster_id %in% c("t2", "t6", "t7")]
somede_dn <- c(4)
dat[somede_dn, colD$group == "A" & colD$cluster_id %in% c("t2", "t6", "t7")] <-
    3 * dat[somede_dn, colD$group == "A" & colD$cluster_id %in% c("t2", "t6", "t7")]

## Genes that are DE for cell type t4 and t5 (not in the same subtree)
onede_up <- c(7)
dat[onede_up, colD$group == "B" & colD$cluster_id %in% c("t4", "t5")] <-
    3 * dat[onede_up, colD$group == "B" & colD$cluster_id %in% c("t4", "t5")]
onede_dn <- c(8)
dat[onede_dn, colD$group == "A" & colD$cluster_id %in% c("t4", "t5")] <-
    3 * dat[onede_dn, colD$group == "A" & colD$cluster_id %in% c("t4", "t5")]

## TreeSummarizedExperiment object
tse <- TreeSummarizedExperiment(
    assays = list(counts = dat),
    colData = colD,
    colTree = tinyTree,
    colNodeLab = clus)

rowData(tse)$Signal <- "no"
rowData(tse)$Signal[c(allde_up, allde_dn)] <- "all"
rowData(tse)$Signal[c(somede_up, somede_dn)] <- "t2, t6, t7"
rowData(tse)$Signal[c(onede_up, onede_dn)] <- "t4, t5"

saveRDS(tse, file = "inst/extdata/ds_sim_20_500_8de.rds")
