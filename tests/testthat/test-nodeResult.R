test_that("nodeResult works", {
    library(TreeSummarizedExperiment)
    ## Generate some example data - DA
    da_lse <- readRDS(system.file("extdata", "da_sim_100_30_18de.rds",
                                  package = "treeclimbR"))
    da_tse <- aggTSE(da_lse, rowLevel = showNode(tree = rowTree(da_lse),
                                                 only.leaf = FALSE))
    da_res <- runDA(TSE = da_tse, feature_on_row = TRUE,
                    assay = "counts", option = "glmQL",
                    design = model.matrix( ~ group, data = colData(da_tse)),
                    contrast = NULL, normalize = TRUE)

    ## Generate some example data - DS
    ds_tse <- readRDS(system.file("extdata", "ds_sim_20_500_8de.rds",
                                  package = "treeclimbR"))
    expect_warning({
        ds_se <- aggDS(TSE = ds_tse, assay = "counts", sample_id = "sample_id",
                       group_id = "group", cluster_id = "cluster_id", FUN = sum)
    })
    expect_warning({
        ds_res <- runDS(SE = ds_se, tree = colTree(ds_tse), option = "glmQL",
                        group_column = "group", contrast = c(0, 1),
                        filter_min_count = 0, filter_min_total_count = 1,
                        design = model.matrix(~ group, data = colData(ds_se)),
                        filter_min_prop = 0, min_cells = 5, message = FALSE)
    })

    ## Check that function returns error with invalid input
    ## -------------------------------------------------------------------------
    expect_error(nodeResult(object = 1, n = 10, type = "DA",
                            adjust_method = "BH", sort_by = "PValue",
                            p_value = 1),
                 "'object' must be of class 'list'")
    expect_error(nodeResult(object = da_res, n = "10", type = "DA",
                            adjust_method = "BH", sort_by = "PValue",
                            p_value = 1),
                 "'n' must be of class 'numeric'")
    expect_error(nodeResult(object = da_res, n = c(1, 2), type = "DA",
                            adjust_method = "BH", sort_by = "PValue",
                            p_value = 1),
                 "'n' must have length 1")
    expect_error(nodeResult(object = da_res, n = 10, type = 1,
                            adjust_method = "BH", sort_by = "PValue",
                            p_value = 1),
                 "'arg' must be NULL or a character vector")
    expect_error(nodeResult(object = da_res, n = 10, type = "missing",
                            adjust_method = "BH", sort_by = "PValue",
                            p_value = 1),
                 "'arg' should be one of ")
    expect_error(nodeResult(object = da_res, n = 10, type = "DA",
                            adjust_method = 1, sort_by = "PValue",
                            p_value = 1),
                 "'adjust_method' must be of class 'character'")
    expect_error(nodeResult(object = da_res, n = 10, type = "DS",
                            adjust_method = c("BH", "holm"), sort_by = "PValue",
                            p_value = 1),
                 "'adjust_method' must have length 1")
    expect_error(nodeResult(object = da_res, n = 10, type = "DA",
                            adjust_method = "BH", sort_by = 1,
                            p_value = 1),
                 "'sort_by' must be of class 'character'")
    expect_error(nodeResult(object = da_res, n = 10, type = "DA",
                            adjust_method = "BH",
                            sort_by = c("PValue", "logFC"), p_value = 1),
                 "'sort_by' must have length 1")
    expect_error(nodeResult(object = da_res, n = 10, type = "DS",
                            adjust_method = "BH", sort_by = "missing",
                            p_value = 1),
                 "All values in 'sort_by' must be one of")
    expect_error(nodeResult(object = da_res, n = 10, type = "DA",
                            adjust_method = "BH", sort_by = "PValue",
                            p_value = TRUE),
                 "'p_value' must be of class 'numeric'")
    expect_error(nodeResult(object = da_res, n = 10, type = "DA",
                            adjust_method = "BH", sort_by = "PValue",
                            p_value = c(0.5, 1)),
                 "'p_value' must have length 1")

    ## Check that function works as expected with valid input - DA
    ## -------------------------------------------------------------------------
    ## Sort by p-value
    out <- nodeResult(object = da_res, n = 8, type = "DA",
                      adjust_method = "BH", sort_by = "PValue",
                      p_value = 1)
    expect_s3_class(out, "data.frame")
    expect_equal(nrow(out), 8)
    expect_named(out, c("node", "logFC", "logCPM", "F", "PValue", "FDR"))
    expect_true(all(diff(out$PValue) >= 0))
    expect_equal(out$node, c(102, 114, 115, 103, 116, 118, 110, 101))

    ## Sort by logFC
    out <- nodeResult(object = da_res, n = 8, type = "DA",
                      adjust_method = "BH", sort_by = "logFC",
                      p_value = 1)
    expect_s3_class(out, "data.frame")
    expect_equal(nrow(out), 8)
    expect_named(out, c("node", "logFC", "logCPM", "F", "PValue", "FDR"))
    expect_true(all(diff(abs(out$logFC)) <= 0))
    expect_equal(out$node, c(120, 112, 122, 9, 108, 10, 7, 121))

    ## No sorting, set FDR threshold
    out <- nodeResult(object = da_res, n = Inf, type = "DA",
                      adjust_method = "BH", sort_by = "none",
                      p_value = 1e-5)
    expect_s3_class(out, "data.frame")
    expect_equal(nrow(out), 11)
    expect_named(out, c("node", "logFC", "logCPM", "F", "PValue", "FDR"))
    expect_false(all(diff(abs(out$logFC)) <= 0))
    expect_equal(out$node, c(101, 102, 103, 104, 110, 112, 114, 115,
                             116, 118, 120))

    ## No sorting, no FDR threshold
    out <- nodeResult(object = da_res, n = Inf, type = "DA",
                      adjust_method = "BH", sort_by = "none",
                      p_value = 1)
    expect_s3_class(out, "data.frame")
    expect_equal(nrow(out), nrow(da_res$edgeR_results))
    expect_named(out, c("node", "logFC", "logCPM", "F", "PValue", "FDR"))
    expect_false(all(diff(abs(out$logFC)) <= 0))
    expect_equal(rownames(out), rownames(da_res$edgeR_results))

    ## Check that function works as expected with valid input - DS
    ## -------------------------------------------------------------------------
    ## Sort by p-value
    out <- nodeResult(object = ds_res, n = 8, type = "DS",
                      adjust_method = "BH", sort_by = "PValue",
                      p_value = 1)
    expect_s3_class(out, "data.frame")
    expect_equal(nrow(out), 8)
    expect_named(out, c("logFC", "logCPM", "F", "PValue", "FDR", "node",
                        "feature"))
    expect_true(all(diff(out$PValue) >= 0))
    expect_equal(out$node, c(11, 13, 13, 11, 11, 11, 13, 12))
    expect_equal(out$feature, as.character(c(3, 5, 6, 1, 6, 5, 2, 3)))

    ## Sort by logFC
    out <- nodeResult(object = ds_res, n = 8, type = "DS",
                      adjust_method = "BH", sort_by = "logFC",
                      p_value = 1)
    expect_s3_class(out, "data.frame")
    expect_equal(nrow(out), 8)
    expect_named(out, c("logFC", "logCPM", "F", "PValue", "FDR", "node",
                        "feature"))
    expect_true(all(diff(abs(out$logFC)) <= 0))
    expect_equal(out$node, c(9, 7, 4, 1, 8, 9, 19, 2))
    expect_equal(out$feature, as.character(c(3, 1, 1, 3, 3, 8, 1, 6)))

    ## No sorting, set FDR threshold
    out <- nodeResult(object = ds_res, n = 8, type = "DS",
                      adjust_method = "BH", sort_by = "none",
                      p_value = 1e-5)
    expect_s3_class(out, "data.frame")
    expect_equal(nrow(out), 8)
    expect_named(out, c("logFC", "logCPM", "F", "PValue", "FDR", "node",
                        "feature"))
    expect_false(all(diff(abs(out$logFC)) <= 0))
    expect_equal(out$node, c(1, 1, 1, 1, 1, 1, 3, 3))
    expect_equal(out$feature, as.character(c(1, 2, 3, 4, 5, 6, 1, 2)))
})
