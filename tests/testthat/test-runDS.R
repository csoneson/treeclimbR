test_that("runDS works", {
    ## Load example data
    ds_tse <- readRDS(system.file("extdata", "ds_sim_20_500_8de.rds",
                                  package = "treeclimbR"))
    expect_warning({
        ds_se <- aggDS(TSE = ds_tse, assay = "counts", sample_id = "sample_id",
                       group_id = "group", cluster_id = "cluster_id", FUN = sum)
    })
    des <- model.matrix(~ group, data = SummarizedExperiment::colData(ds_se))

    ## Check that function returns error for invalid input
    ## -------------------------------------------------------------------------
    .args <- list(SE = ds_se, tree = colTree(ds_tse),
                  option = "glm", design = des, contrast = c(0, 1),
                  filter_min_count = 1, filter_min_total_count = 15,
                  filter_large_n = 10, filter_min_prop = 1,
                  min_cells = 10, normalize = TRUE, normalize_method = "TMM",
                  group_column = "group", design_terms = "group",
                  message = TRUE)

    args <- .args
    args$SE <- 1
    expect_error(do.call(runDS, args),
                 "'SE' must be of class 'SummarizedExperiment'")

    args <- .args
    args$tree <- 1
    expect_error(do.call(runDS, args),
                 "'tree' must be of class 'phylo'")

    args <- .args
    args$option <- 1
    expect_error(do.call(runDS, args),
                 "'arg' must be NULL or a character vector")
    args$option <- "missing"
    expect_error(do.call(runDS, args),
                 "'arg' should be one of")

    args <- .args
    args$design <- "x"
    expect_error(do.call(runDS, args),
                 "'design' must be of class 'matrix'")

    args <- .args
    args$contrast <- c("x", "y")
    expect_error(do.call(runDS, args),
                 "'contrast' must be of class 'numeric'")

    args <- .args
    args$filter_min_count <- "x"
    expect_error(do.call(runDS, args),
                 "'filter_min_count' must be of class 'numeric'")
    args$filter_min_count <- c(0.1, 0.2)
    expect_error(do.call(runDS, args),
                 "'filter_min_count' must have length 1")

    args <- .args
    args$filter_min_total_count <- "x"
    expect_error(do.call(runDS, args),
                 "'filter_min_total_count' must be of class 'numeric'")
    args$filter_min_total_count <- c(0.1, 0.2)
    expect_error(do.call(runDS, args),
                 "'filter_min_total_count' must have length 1")

    args <- .args
    args$filter_large_n <- "x"
    expect_error(do.call(runDS, args),
                 "'filter_large_n' must be of class 'numeric'")
    args$filter_large_n <- c(0.1, 0.2)
    expect_error(do.call(runDS, args),
                 "'filter_large_n' must have length 1")

    args <- .args
    args$filter_min_prop <- "x"
    expect_error(do.call(runDS, args),
                 "'filter_min_prop' must be of class 'numeric'")
    args$filter_min_prop <- c(0.1, 0.2)
    expect_error(do.call(runDS, args),
                 "'filter_min_prop' must have length 1")

    args <- .args
    args$min_cells <- "x"
    expect_error(do.call(runDS, args),
                 "'min_cells' must be of class 'numeric'")
    args$min_cells <- c(0.1, 0.2)
    expect_error(do.call(runDS, args),
                 "'min_cells' must have length 1")

    args <- .args
    args$normalize <- "x"
    expect_error(do.call(runDS, args),
                 "'normalize' must be of class 'logical'")
    args$normalize <- c(TRUE, FALSE)
    expect_error(do.call(runDS, args),
                 "'normalize' must have length 1")

    args <- .args
    args$normalize_method <- 1
    expect_error(do.call(runDS, args),
                 "'normalize_method' must be of class 'character'")
    args$normalize_method <- c("none", "TMM")
    expect_error(do.call(runDS, args),
                 "'normalize_method' must have length 1")

    args <- .args
    args$group_column <- 1
    expect_error(do.call(runDS, args),
                 "'group_column' must be of class 'character'")

    args <- .args
    args$design_terms <- 1
    expect_error(do.call(runDS, args),
                 "'design_terms' must be of class 'character'")

    args <- .args
    args$message <- 1
    expect_error(do.call(runDS, args),
                 "'message' must be of class 'logical'")
    args$message <- c(TRUE, FALSE)
    expect_error(do.call(runDS, args),
                 "'message' must have length 1")

    ## Check that function works as expected for valid input
    ## -------------------------------------------------------------------------
    ## glmQL
    ds_res <- runDS(SE = ds_se, tree = colTree(ds_tse), option = "glmQL",
                    group_column = "group", contrast = c(0, 1),
                    filter_min_count = 0, filter_min_total_count = 1,
                    design = model.matrix(~ group, data = colData(ds_se)),
                    filter_min_prop = 0, min_cells = 5, message = FALSE,
                    legacy = FALSE)
    expect_type(ds_res, "list")
    expect_named(ds_res, c("edgeR_results", "tree", "nodes_drop"))
    expect_equal(ds_res$nodes_drop, character(0))
    expect_type(ds_res$edgeR_results, "list")
    expect_length(ds_res$edgeR_results, 19)
    expect_s4_class(ds_res$edgeR_results[[1]], "DGELRT")
    ## Check that truly DE features show up where they should
    expect_equal(sort(rownames(edgeR::topTags(ds_res$edgeR_results$alias_11,
                                              n = Inf, p.value = 1e-4))),
                 as.character(seq_len(8)))
    expect_equal(sort(rownames(edgeR::topTags(ds_res$edgeR_results$alias_13,
                                              n = Inf, p.value = 1e-4))),
                 as.character(seq_len(6)))
    expect_equal(sort(rownames(edgeR::topTags(ds_res$edgeR_results$alias_14,
                                              n = Inf, p.value = 1e-4))),
                 as.character(seq_len(6)))
    expect_equal(sort(rownames(edgeR::topTags(ds_res$edgeR_results$alias_2,
                                              n = Inf, p.value = 1e-2))),
                 as.character(seq_len(6)))
    expect_equal(sort(rownames(edgeR::topTags(ds_res$edgeR_results$alias_7,
                                              n = Inf, p.value = 1e-2))),
                 as.character(c(1, 3, 5, 6)))
    expect_equal(sort(rownames(edgeR::topTags(ds_res$edgeR_results$alias_19,
                                              n = Inf, p.value = 1e-2))),
                 as.character(c(1, 3, 5, 6)))
    expect_equal(sort(rownames(edgeR::topTags(ds_res$edgeR_results$alias_5,
                                              n = Inf, p.value = 1e-4))),
                 as.character(c(1, 3, 5, 6, 7, 8)))

    ## glm
    ds_res <- runDS(SE = ds_se, tree = colTree(ds_tse), option = "glm",
                    group_column = "group", contrast = c(0, 1),
                    filter_min_count = 0, filter_min_total_count = 1,
                    design = model.matrix(~ group, data = colData(ds_se)),
                    filter_min_prop = 0, min_cells = 5, message = FALSE)
    expect_type(ds_res, "list")
    expect_named(ds_res, c("edgeR_results", "tree", "nodes_drop"))
    expect_equal(ds_res$nodes_drop, character(0))
    expect_type(ds_res$edgeR_results, "list")
    expect_length(ds_res$edgeR_results, 19)
    expect_s4_class(ds_res$edgeR_results[[1]], "DGELRT")
    ## Check that truly DE features show up where they should
    expect_equal(sort(rownames(edgeR::topTags(ds_res$edgeR_results$alias_11,
                                              n = Inf, p.value = 1e-4))),
                 as.character(seq_len(8)))
    expect_equal(sort(rownames(edgeR::topTags(ds_res$edgeR_results$alias_13,
                                              n = Inf, p.value = 1e-4))),
                 as.character(seq_len(6)))
    expect_equal(sort(rownames(edgeR::topTags(ds_res$edgeR_results$alias_14,
                                              n = Inf, p.value = 1e-4))),
                 as.character(seq_len(6)))
    expect_equal(sort(rownames(edgeR::topTags(ds_res$edgeR_results$alias_2,
                                              n = Inf, p.value = 1e-2))),
                 as.character(seq_len(6)))
    expect_equal(sort(rownames(edgeR::topTags(ds_res$edgeR_results$alias_7,
                                              n = Inf, p.value = 1e-2))),
                 as.character(c(1, 3, 5, 6)))
    expect_equal(sort(rownames(edgeR::topTags(ds_res$edgeR_results$alias_19,
                                              n = Inf, p.value = 1e-2))),
                 as.character(c(1, 3, 5, 6)))
    expect_equal(sort(rownames(edgeR::topTags(ds_res$edgeR_results$alias_5,
                                              n = Inf, p.value = 1e-4))),
                 as.character(c(1, 3, 5, 6, 7, 8)))

    ## glmQL, filter by min_cells
    ds_res <- runDS(SE = ds_se, tree = colTree(ds_tse), option = "glmQL",
                    group_column = "group", contrast = c(0, 1),
                    filter_min_count = 0, filter_min_total_count = 1,
                    design = model.matrix(~ group, data = colData(ds_se)),
                    filter_min_prop = 0, min_cells = 13, message = FALSE,
                    legacy = FALSE)
    expect_type(ds_res, "list")
    expect_named(ds_res, c("edgeR_results", "tree", "nodes_drop"))
    expect_equal(ds_res$nodes_drop, c("alias_1", "alias_3", "alias_10"))
    expect_type(ds_res$edgeR_results, "list")
    expect_length(ds_res$edgeR_results, 16)
    expect_s4_class(ds_res$edgeR_results[[1]], "DGELRT")
    ## Check that truly DE features show up where they should
    expect_equal(sort(rownames(edgeR::topTags(ds_res$edgeR_results$alias_11,
                                              n = Inf, p.value = 1e-4))),
                 as.character(seq_len(8)))
    expect_equal(sort(rownames(edgeR::topTags(ds_res$edgeR_results$alias_13,
                                              n = Inf, p.value = 1e-4))),
                 as.character(seq_len(6)))
    expect_equal(sort(rownames(edgeR::topTags(ds_res$edgeR_results$alias_14,
                                              n = Inf, p.value = 1e-4))),
                 as.character(seq_len(6)))
    expect_equal(sort(rownames(edgeR::topTags(ds_res$edgeR_results$alias_2,
                                              n = Inf, p.value = 1e-2))),
                 as.character(seq_len(6)))
    expect_equal(sort(rownames(edgeR::topTags(ds_res$edgeR_results$alias_7,
                                              n = Inf, p.value = 1e-2))),
                 as.character(c(1, 3, 5, 6)))
    expect_equal(sort(rownames(edgeR::topTags(ds_res$edgeR_results$alias_19,
                                              n = Inf, p.value = 1e-2))),
                 as.character(c(1, 3, 5, 6)))
    expect_equal(sort(rownames(edgeR::topTags(ds_res$edgeR_results$alias_5,
                                              n = Inf, p.value = 1e-4))),
                 as.character(c(1, 3, 5, 6, 7, 8)))

    ## glmQL, some nodes don't have samples from all groups
    ds2 <- ds_se
    SummarizedExperiment::assay(ds2, "alias_5")[, c("1", "2")] <- 0
    SummarizedExperiment::assay(ds2, "alias_13")[, c("3", "4")] <- 0
    ds_res <- runDS(SE = ds2, tree = colTree(ds_tse), option = "glmQL",
                    group_column = "group", contrast = c(0, 1),
                    filter_min_count = 0, filter_min_total_count = 1,
                    design = model.matrix(~ group, data = colData(ds_se)),
                    filter_min_prop = 0, min_cells = 1, message = FALSE,
                    legacy = FALSE)
    expect_type(ds_res, "list")
    expect_named(ds_res, c("edgeR_results", "tree", "nodes_drop"))
    expect_equal(ds_res$nodes_drop, c("alias_5", "alias_13"),
                 ignore_attr = TRUE)
    expect_type(ds_res$edgeR_results, "list")
    expect_length(ds_res$edgeR_results, 17)
    expect_s4_class(ds_res$edgeR_results[[1]], "DGELRT")
    ## Check that truly DE features show up where they should
    expect_equal(sort(rownames(edgeR::topTags(ds_res$edgeR_results$alias_11,
                                              n = Inf, p.value = 1e-4))),
                 as.character(seq_len(8)))
    expect_null(ds_res$edgeR_results$alias_13)
    expect_equal(sort(rownames(edgeR::topTags(ds_res$edgeR_results$alias_14,
                                              n = Inf, p.value = 1e-4))),
                 as.character(seq_len(6)))
    expect_equal(sort(rownames(edgeR::topTags(ds_res$edgeR_results$alias_2,
                                              n = Inf, p.value = 1e-2))),
                 as.character(seq_len(6)))
    expect_equal(sort(rownames(edgeR::topTags(ds_res$edgeR_results$alias_7,
                                              n = Inf, p.value = 1e-2))),
                 as.character(c(1, 3, 5, 6)))
    expect_equal(sort(rownames(edgeR::topTags(ds_res$edgeR_results$alias_19,
                                              n = Inf, p.value = 1e-2))),
                 as.character(c(1, 3, 5, 6)))
    expect_null(ds_res$edgeR_results$alias_5)

    ## glmQL, some samples have zero counts in some nodes,
    ## specify design_terms
    ds2 <- ds_se
    SummarizedExperiment::assay(ds2, "alias_5")[, c("1")] <- 0
    SummarizedExperiment::assay(ds2, "alias_13")[, c("4")] <- 0
    ds_res <- runDS(SE = ds2, tree = colTree(ds_tse), option = "glmQL",
                    group_column = "group", contrast = c(0, 1),
                    filter_min_count = 0, filter_min_total_count = 1,
                    design = NULL, design_terms = "group",
                    filter_min_prop = 0, min_cells = 1, message = FALSE,
                    legacy = FALSE)
    expect_type(ds_res, "list")
    expect_named(ds_res, c("edgeR_results", "tree", "nodes_drop"))
    expect_equal(ds_res$nodes_drop, character(0),
                 ignore_attr = TRUE)
    expect_type(ds_res$edgeR_results, "list")
    expect_length(ds_res$edgeR_results, 19)
    expect_s4_class(ds_res$edgeR_results[[1]], "DGELRT")
    ## Check that truly DE features show up where they should
    expect_equal(sort(rownames(edgeR::topTags(ds_res$edgeR_results$alias_11,
                                              n = Inf, p.value = 1e-4))),
                 as.character(seq_len(8)))
    expect_equal(sort(rownames(edgeR::topTags(ds_res$edgeR_results$alias_14,
                                              n = Inf, p.value = 1e-4))),
                 as.character(seq_len(6)))
    expect_equal(sort(rownames(edgeR::topTags(ds_res$edgeR_results$alias_2,
                                              n = Inf, p.value = 1e-2))),
                 as.character(seq_len(6)))
    expect_equal(sort(rownames(edgeR::topTags(ds_res$edgeR_results$alias_7,
                                              n = Inf, p.value = 1e-2))),
                 as.character(c(1, 3, 5, 6)))
    expect_equal(sort(rownames(edgeR::topTags(ds_res$edgeR_results$alias_19,
                                              n = Inf, p.value = 1e-2))),
                 as.character(c(1, 3, 5, 6)))
    expect_equal(sort(rownames(edgeR::topTags(ds_res$edgeR_results$alias_13,
                                              n = Inf, p.value = 1e-4))),
                 as.character(seq_len(6)))
    expect_equal(sort(rownames(edgeR::topTags(ds_res$edgeR_results$alias_5,
                                              n = Inf, p.value = 1e-4))),
                 as.character(c(1, 3, 5, 6, 7, 8)))
})
