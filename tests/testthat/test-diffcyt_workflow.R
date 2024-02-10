test_that("diffcyt workflow works", {
    ## Generate some data and run through diffcyt workflow
    ## -------------------------------------------------------------------------
    library(diffcyt)
    ## Helper function to create random data (one sample)
    d_random <- function(n = 20000, mean = 0, sd = 1, ncol = 20, cofactor = 5) {
        d <- sinh(matrix(rnorm(n, mean, sd), ncol = ncol)) * cofactor
        colnames(d) <- paste0("marker", sprintf("%02d", seq_len(ncol)))
        d
    }
    ## Create random data (without differential signal)
    set.seed(123)
    d_input <- list(
        sample1 = d_random(), sample2 = d_random(),
        sample3 = d_random(), sample4 = d_random()
    )
    experiment_info <- data.frame(
        sample_id = factor(paste0("sample", seq_len(4))),
        group_id = factor(c("group1", "group1", "group2", "group2"))
    )
    marker_info <- data.frame(
        channel_name = paste0("channel", sprintf("%03d", seq_len(20))),
        marker_name = paste0("marker", sprintf("%02d", seq_len(20))),
        marker_class = factor(c(rep("type", 10), rep("state", 10)),
                              levels = c("type", "state", "none"))
    )
    d_se <- diffcyt::prepareData(d_input, experiment_info, marker_info)
    d_se <- diffcyt::transformData(d_se)
    d_se <- diffcyt::generateClusters(d_se)

    ## buildTree - check that function errors if provided wrong input
    ## -------------------------------------------------------------------------
    expect_error(buildTree(d_se = 1, dist_method = "euclidean",
                           hclust_method = "average"),
                 "'d_se' must be of class 'SummarizedExperiment'")
    expect_error(buildTree(d_se = d_se, dist_method = 1,
                           hclust_method = "average"),
                 "'dist_method' must be of class 'character'")
    expect_error(buildTree(d_se = d_se, dist_method = c("euclidean", "cosine"),
                           hclust_method = "average"),
                 "'dist_method' must have length 1")
    expect_error(buildTree(d_se = d_se, dist_method = "euclidean",
                           hclust_method = 1),
                 "'hclust_method' must be of class 'character'")
    expect_error(buildTree(d_se = d_se, dist_method = "missing",
                           hclust_method = "average"),
                 "invalid distance method")
    expect_error(buildTree(d_se = d_se, dist_method = "euclidean",
                           hclust_method = "missing"),
                 "invalid clustering method")
    expect_error(buildTree(d_se = d_se, dist_method = "euclidean",
                           hclust_method = c("average", "complete")),
                 "'hclust_method' must have length 1")

    ## buildTree - check that function works with correct input
    ## ------------------------------------------------------------------------
    tr <- buildTree(d_se, dist_method = "euclidean", hclust_method = "average")
    expect_s3_class(tr, "phylo")
    expect_equal(tr$tip.label, as.character(seq_len(100)))
    expect_equal(tr$node.label, paste0("Node_", seq(from = 101, to = 199,
                                                    by = 1)))
    expect_equal(tr$Nnode, 99)

    tr <- buildTree(d_se, dist_method = "euclidean", hclust_method = "complete")
    expect_s3_class(tr, "phylo")
    expect_equal(tr$tip.label, as.character(seq_len(100)))
    expect_equal(tr$node.label, paste0("Node_", seq(from = 101, to = 199,
                                                    by = 1)))
    expect_equal(tr$Nnode, 99)

    ## Build a tree to use in downstream functions
    ## -------------------------------------------------------------------------
    tr <- buildTree(d_se, dist_method = "euclidean", hclust_method = "average")

    ## calcTreeMedians - check that function errors if provided wrong input
    ## -------------------------------------------------------------------------
    expect_error(calcTreeMedians(d_se = 1, tree = tr, message = FALSE),
                 "'d_se' must be of class 'SummarizedExperiment'")
    expect_error(calcTreeMedians(d_se = d_se, tree = 1, message = FALSE),
                 "tree is not a phylo object")
    expect_error(calcTreeMedians(d_se = d_se, tree = tr, message = 1),
                 "'message' must be of class 'logical'")
    expect_error(calcTreeMedians(d_se = d_se, tree = tr,
                                 message = c(TRUE, FALSE)),
                 "'message' must have length 1")
    tmp <- d_se
    SummarizedExperiment::rowData(tmp)$cluster_id <- NULL
    expect_error(calcTreeMedians(d_se = tmp, tree = tr, message = TRUE),
                 "Data object does not contain cluster labels")
    tmp <- d_se
    SummarizedExperiment::colData(tmp)$marker_class <- "unknown"
    expect_warning({
        expect_warning({
            expect_warning({
                expect_error(calcTreeMedians(d_se = tmp, tree = tr,
                                             message = FALSE),
                             "No type or state markers found in the object")
            }, "Multiple nodes are found to have the same label")
        }, "Multiple nodes are found to have the same label")
    }, "Multiple nodes are found to have the same label")

    ## calcTreeMedians - check that function works with correct input
    ## ------------------------------------------------------------------------
    expect_warning({
        expect_warning({
            expect_warning({
                out <- calcTreeMedians(d_se = d_se, tree = tr, message = FALSE)
            }, "Multiple nodes are found to have the same label")
        }, "Multiple nodes are found to have the same label")
    }, "Multiple nodes are found to have the same label")
    expect_s4_class(out, "TreeSummarizedExperiment")
    expect_equal(length(SummarizedExperiment::assayNames(out)), 20) ## nmarkers
    expect_equal(nrow(out), 199) ## nclusters
    expect_equal(ncol(out), 4) ## nsamples
    ## Spot check some values
    for (l in list(list(node = 127, sample = "sample1", marker = "marker17"),
                   list(node = 45, sample = "sample2", marker = "marker06"))) {
        leaves <- TreeSummarizedExperiment::findDescendant(
            TreeSummarizedExperiment::rowTree(out), node = l$node,
            only.leaf = TRUE, self.include = TRUE)[[1]]
        idx <- which(rowData(d_se)$cluster_id %in% leaves &
                         rowData(d_se)$sample_id == l$sample)
        expect_equal(
            SummarizedExperiment::assay(out, l$marker)[paste0("alias_", l$node),
                                                       l$sample],
            stats::median(
                SummarizedExperiment::assay(d_se, "exprs")[idx, l$marker]))
    }

    ## Run with message = TRUE
    expect_warning({
        expect_warning({
            expect_warning({
                out2 <- calcTreeMedians(d_se = d_se, tree = tr, message = TRUE)
            }, "Multiple nodes are found to have the same label")
        }, "Multiple nodes are found to have the same label")
    }, "Multiple nodes are found to have the same label")
    expect_identical(out, out2)

    ## With missing nodes
    set.seed(1)
    kp <- sample(seq_len(nrow(d_se)), 3000)
    tmp <- d_se[kp, ]
    tmp2 <- d_se
    SummarizedExperiment::assay(tmp2, "exprs")[
        setdiff(seq_len(nrow(d_se)), kp), ] <- NA
    expect_warning({
        expect_warning({
            expect_warning({
                expect_warning({
                    expect_warning({
                        out <- calcTreeMedians(d_se = tmp, tree = tr,
                                               message = FALSE)
                    }, "Multiple nodes are found to have the same label")
                }, "Multiple nodes are found to have the same label")
            }, "Multiple nodes are found to have the same label")
        }, "leaves couldn't be found")
    }, "Missing leaves")
    expect_s4_class(out, "TreeSummarizedExperiment")
    expect_equal(length(SummarizedExperiment::assayNames(out)), 20) ## nmarkers
    expect_equal(nrow(out), 199) ## nclusters
    expect_equal(ncol(out), 4) ## nsamples
    ## Spot check some values
    for (l in list(list(node = 127, sample = "sample1", marker = "marker17"),
                   list(node = 45, sample = "sample2", marker = "marker06"))) {
        leaves <- TreeSummarizedExperiment::findDescendant(
            TreeSummarizedExperiment::rowTree(out), node = l$node,
            only.leaf = TRUE, self.include = TRUE)[[1]]
        idx <- which(rowData(tmp2)$cluster_id %in% leaves &
                         rowData(tmp2)$sample_id == l$sample)
        expect_equal(
            SummarizedExperiment::assay(out, l$marker)[paste0("alias_", l$node),
                                                       l$sample],
            stats::median(
                SummarizedExperiment::assay(tmp2, "exprs")[idx, l$marker],
                na.rm = TRUE))
    }

    ## calcTreeCounts - check that function errors if provided wrong input
    ## -------------------------------------------------------------------------
    expect_error(calcTreeCounts(d_se = 1, tree = tr),
                 "'d_se' must be of class 'SummarizedExperiment'")
    expect_error(calcTreeCounts(d_se = d_se, tree = 1),
                 "tree is not a phylo object")
    tmp <- d_se
    SummarizedExperiment::rowData(tmp)$cluster_id <- NULL
    expect_error(calcTreeCounts(d_se = tmp, tree = tr),
                 "Data object does not contain cluster labels")

    ## calcTreeCounts - check that function works with correct input
    ## ------------------------------------------------------------------------
    out <- calcTreeCounts(d_se = d_se, tree = tr)
    expect_s4_class(out, "TreeSummarizedExperiment")
    expect_equal(nrow(out), 199) ## nclusters
    expect_equal(ncol(out), 4) ## nsamples
    expect_equal(SummarizedExperiment::assayNames(out), "counts")
    ## Spot check some values
    for (l in list(list(node = 127, sample = "sample1"),
                   list(node = 45, sample = "sample2"))) {
        leaves <- TreeSummarizedExperiment::findDescendant(
            TreeSummarizedExperiment::rowTree(out), node = l$node,
            only.leaf = TRUE, self.include = TRUE)[[1]]
        idx <- which(rowData(d_se)$cluster_id %in% leaves &
                         rowData(d_se)$sample_id == l$sample)
        expect_equal(
            SummarizedExperiment::assay(out, "counts")[paste0("alias_", l$node),
                                                       l$sample],
            length(idx))
    }

    ## calcMediansByTreeMarker - check that function errors if provided wrong
    ## input
    ## -------------------------------------------------------------------------
    expect_error(calcMediansByTreeMarker(d_se = 1, tree = tr),
                 "'d_se' must be of class 'SummarizedExperiment'")
    expect_error(calcMediansByTreeMarker(d_se = d_se, tree = 1),
                 "tree is not a phylo object")
    tmp <- d_se
    SummarizedExperiment::rowData(tmp)$cluster_id <- NULL
    expect_error(calcMediansByTreeMarker(d_se = tmp, tree = tr),
                 "Data object does not contain cluster labels")
    tmp <- d_se
    SummarizedExperiment::colData(tmp)$marker_class <- "unknown"
    expect_warning({
        expect_warning({
            expect_warning({
                expect_error(calcMediansByTreeMarker(d_se = tmp, tree = tr),
                             "No type or state markers found in the object")
            }, "Multiple nodes are found to have the same label")
        }, "Multiple nodes are found to have the same label")
    }, "Multiple nodes are found to have the same label")

    ## calcMediansByTreeMarker - check that function works with correct input
    ## ------------------------------------------------------------------------
    expect_warning({
        expect_warning({
            expect_warning({
                out <- calcMediansByTreeMarker(d_se = d_se, tree = tr)
            }, "Multiple nodes are found to have the same label")
        }, "Multiple nodes are found to have the same label")
    }, "Multiple nodes are found to have the same label")
    expect_s4_class(out, "TreeSummarizedExperiment")
    expect_equal(nrow(out), 199) ## nclusters
    expect_equal(ncol(out), 20) ## nmarkers
    expect_equal(SummarizedExperiment::assayNames(out), "exprs")
    ## Spot check some values
    for (l in list(list(node = 127, marker = "marker17"),
                   list(node = 45, marker = "marker06"))) {
        leaves <- TreeSummarizedExperiment::findDescendant(
            TreeSummarizedExperiment::rowTree(out), node = l$node,
            only.leaf = TRUE, self.include = TRUE)[[1]]
        idx <- which(rowData(d_se)$cluster_id %in% leaves)
        expect_equal(
            SummarizedExperiment::assay(out, "exprs")[l$node, l$marker],
            stats::median(
                SummarizedExperiment::assay(d_se, "exprs")[idx, l$marker]))
    }
})
