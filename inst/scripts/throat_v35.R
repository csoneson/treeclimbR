## The throat_v35 dataset consists of a subset of a data set from HMP16SData.
## It contains 16S OTU counts from throat samples.
## It is provided as a TreeSummarizedExperiment object with 956 OTUs and 153 samples.

library(dirmult)
library(TreeSummarizedExperiment)
library(ape)
library(HMP16SData)
library(treeclimbR)
library(dplyr)

w <- V35()

xt <- table(colData(w)$HMP_BODY_SUBSITE)
names(xt)
# [1] "Anterior Nares"               "Attached Keratinized Gingiva" "Buccal Mucosa"
# [4] "Hard Palate"                  "Left Antecubital Fossa"       "Left Retroauricular Crease"
# [7] "Mid Vagina"                   "Palatine Tonsils"             "Posterior Fornix"
# [10] "Right Antecubital Fossa"      "Right Retroauricular Crease"  "Saliva"
# [13] "Stool"                        "Subgingival Plaque"           "Supragingival Plaque"
# [16] "Throat"                       "Tongue Dorsum"                "Vaginal Introitus"
names(xt)[16]
# [1] "Throat"

throat <- w %>%
    subset(select = HMP_BODY_SUBSITE == names(xt)[16])
throat
# class: SummarizedExperiment
# dim: 45383 307
# metadata(2): experimentData phylogeneticTree
# assays(1): 16SrRNA
# rownames(45383): OTU_97.1 OTU_97.10 ... OTU_97.9998 OTU_97.9999
# rowData names(7): CONSENSUS_LINEAGE SUPERKINGDOM ... FAMILY GENUS
# colnames(307): 700014424 700014520 ... 700114337 700114379
# colData names(7): RSID VISITNO ... HMP_BODY_SUBSITE SRS_SAMPLE_ID

# keep OTUs: non-zero count in more than 25% samples
count <- assays(throat)[[1]]
lib <- apply(count, 2, sum)
qt <- quantile(lib, c(0.25, 0.75))

throat <- throat[, lib >= qt[1] & lib <= qt[2]]
count <- assays(throat)[[1]]
nz_samp <- apply(count, 1, FUN = function(x){sum(x != 0)/length(x)})
throat <- throat[nz_samp > 0.25, ]

# tree
phyTree <- metadata(w)$phylogeneticTree
leaf <- phyTree$tip.label

treeR <- drop.tip(phy = phyTree,
                  tip = setdiff(leaf, rownames(throat)),
                  trim.internal = TRUE, collapse.singles = TRUE)

throat <- TreeSummarizedExperiment(assays = assays(throat),
                                   rowData = rowData(throat),
                                   colData = colData(throat),
                                   rowTree = treeR)
throat_v35 <- parEstimate(obj = throat)
metadata(throat_v35)

saveRDS(throat_v35, "inst/extdata/throat_v35.rds")
