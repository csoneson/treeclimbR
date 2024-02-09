test_that("runDA works", {
    x <- readRDS(system.file("extdata/da_sim_100_30_18de.rds",
                             package = "treeclimbR"))
    des <- model.matrix(~ group, data = SummarizedExperiment::colData(x))

    ## Check that function fails with incorrect input
    ## -------------------------------------------------------------------------
    .args <- list(TSE = x, feature_on_row = TRUE, assay = "counts",
                  option = "glm", design = des, contrast = c(0, 1),
                  filter_min_count = 10, filter_min_total_count = 15,
                  filter_large_n = 10, filter_min_prop = 0.7,
                  normalize = TRUE, normalize_method = "TMM",
                  group_column = "group", design_terms = "group")

    args <- .args
    args$TSE <- 1
    expect_error(do.call(runDA, args),
                 "'TSE' must be of class 'TreeSummarizedExperiment'")

    args <- .args
    args$feature_on_row <- 1
    expect_error(do.call(runDA, args),
                 "'feature_on_row' must be of class 'logical'")
    args$feature_on_row <- c(TRUE, FALSE)
    expect_error(do.call(runDA, args),
                 "'feature_on_row' must have length 1")

    args <- .args
    args$assay <- TRUE
    expect_error(do.call(runDA, args),
                 "'assay' must be of class 'numeric' or 'character'")
    args$assay <- "missing"
    expect_error(do.call(runDA, args),
                 "assay %in% SummarizedExperiment::assayNames(TSE)",
                 fixed = TRUE)
    args$assay <- 100
    expect_error(do.call(runDA, args),
                 "assay <= length(SummarizedExperiment::assays(TSE))",
                 fixed = TRUE)

    args <- .args
    args$option <- 1
    expect_error(do.call(runDA, args),
                 "'arg' must be NULL or a character vector")
    args$option <- "missing"
    expect_error(do.call(runDA, args),
                 "'arg' should be one of")

    args <- .args
    args$design <- "x"
    expect_error(do.call(runDA, args),
                 "'design' must be of class 'matrix'")

    args <- .args
    args$contrast <- c("x", "y")
    expect_error(do.call(runDA, args),
                 "'contrast' must be of class 'numeric'")

    args <- .args
    args$filter_min_count <- "x"
    expect_error(do.call(runDA, args),
                 "'filter_min_count' must be of class 'numeric'")
    args$filter_min_count <- c(0.1, 0.2)
    expect_error(do.call(runDA, args),
                 "'filter_min_count' must have length 1")

    args <- .args
    args$filter_min_total_count <- "x"
    expect_error(do.call(runDA, args),
                 "'filter_min_total_count' must be of class 'numeric'")
    args$filter_min_total_count <- c(0.1, 0.2)
    expect_error(do.call(runDA, args),
                 "'filter_min_total_count' must have length 1")

    args <- .args
    args$filter_large_n <- "x"
    expect_error(do.call(runDA, args),
                 "'filter_large_n' must be of class 'numeric'")
    args$filter_large_n <- c(0.1, 0.2)
    expect_error(do.call(runDA, args),
                 "'filter_large_n' must have length 1")

    args <- .args
    args$filter_min_prop <- "x"
    expect_error(do.call(runDA, args),
                 "'filter_min_prop' must be of class 'numeric'")
    args$filter_min_prop <- c(0.1, 0.2)
    expect_error(do.call(runDA, args),
                 "'filter_min_prop' must have length 1")

    args <- .args
    args$normalize <- "x"
    expect_error(do.call(runDA, args),
                 "'normalize' must be of class 'logical'")
    args$normalize <- c(TRUE, FALSE)
    expect_error(do.call(runDA, args),
                 "'normalize' must have length 1")

    args <- .args
    args$normalize_method <- 1
    expect_error(do.call(runDA, args),
                 "'normalize_method' must be of class 'character'")
    args$normalize_method <- c("none", "TMM")
    expect_error(do.call(runDA, args),
                 "'normalize_method' must have length 1")

    args <- .args
    args$design_terms <- 1
    expect_error(do.call(runDA, args),
                 "'design_terms' must be of class 'character'")

    ## Check that function works as expected with valid input
    ## -------------------------------------------------------------------------
    truede <- TreeSummarizedExperiment::rowLinks(x)$nodeLab_alias[
        SummarizedExperiment::rowData(x)$Signal
    ]

    out <- runDA(TSE = x, feature_on_row = TRUE, assay = "counts",
                 option = "glm", design = des, contrast = c(0, 1),
                 filter_min_count = 0, filter_min_total_count = 0,
                 filter_large_n = 0, filter_min_prop = 0,
                 normalize = TRUE, normalize_method = "TMM",
                 group_column = "group", design_terms = "group")
    ## Should get identical results without specifying contrast
    ## (only the 'comparison' string will be different)
    out2 <- runDA(TSE = x, feature_on_row = TRUE, assay = "counts",
                 option = "glm", design = des, contrast = NULL,
                 filter_min_count = 0, filter_min_total_count = 0,
                 filter_large_n = 0, filter_min_prop = 0,
                 normalize = TRUE, normalize_method = "TMM",
                 group_column = "group", design_terms = "group")
    expect_type(out, "list")
    expect_named(out, c("edgeR_results", "nodes_drop", "tree"))
    expect_s4_class(out$edgeR_results, "DGELRT")
    expect_identical(out$edgeR_results$table, out2$edgeR_results$table)
    expect_identical(out$edgeR_results$AveLogCPM, out2$edgeR_results$AveLogCPM)
    expect_identical(out$edgeR_results$dispersion,
                     out2$edgeR_results$dispersion)
    expect_identical(out$edgeR_results$coefficients,
                     out2$edgeR_results$coefficients)
    ## Check that truly differential features are high in the results list
    ## (these should be the same as in the tests for edgerWrp, as in this
    ## case we're only testing the tips)
    pos <- sort(match(truede, rownames(edgeR::topTags(out$edgeR_results,
                                                      n = Inf))))
    expect_equal(pos, c(1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 12, 13, 14, 15,
                        16, 17, 30, 47))

    ## Without normalization factors
    out <- runDA(TSE = x, feature_on_row = TRUE, assay = "counts",
                 option = "glm", design = des, contrast = c(0, 1),
                 filter_min_count = 0, filter_min_total_count = 0,
                 filter_large_n = 0, filter_min_prop = 0,
                 normalize = FALSE, normalize_method = "TMM",
                 group_column = "group", design_terms = "group")
    expect_s4_class(out$edgeR_results, "DGELRT")
    ## Check that truly differential features are high in the results list
    pos <- sort(match(truede, rownames(edgeR::topTags(out$edgeR_results,
                                                      n = Inf))))
    expect_equal(pos, c(2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 15, 17,
                        18, 21, 40, 62))

    ## With the QL framework, with assay = NULL, design = NULL
    out <- runDA(TSE = x, feature_on_row = TRUE, assay = NULL,
                 option = "glmQL", design = NULL, contrast = c(0, 1),
                 filter_min_count = 0, filter_min_total_count = 0,
                 filter_large_n = 0, filter_min_prop = 0,
                 normalize = TRUE, normalize_method = "TMM",
                 group_column = "group", design_terms = "group",
                 legacy = TRUE)
    expect_s4_class(out$edgeR_results, "DGELRT")
    ## Check that truly differential features are high in the results list
    pos <- sort(match(truede, rownames(edgeR::topTags(out$edgeR_results,
                                                      n = Inf))))
    expect_equal(pos, c(1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 13, 14, 15,
                        16, 17, 29, 45))

    ## With some filtering
    out <- runDA(TSE = x, feature_on_row = TRUE, assay = NULL,
                 option = "glmQL", design = NULL, contrast = c(0, 1),
                 filter_min_count = 0, filter_min_total_count = 10,
                 filter_large_n = 0, filter_min_prop = 0,
                 normalize = TRUE, normalize_method = "TMM",
                 group_column = "group", design_terms = "group",
                 legacy = TRUE)
    expect_s4_class(out$edgeR_results, "DGELRT")
    expect_type(out$nodes_drop, "character")
    expect_length(out$nodes_drop, 5)
    ## Check that truly differential features are high in the results list
    pos <- sort(match(truede, rownames(edgeR::topTags(out$edgeR_results,
                                                      n = Inf))))
    expect_equal(pos, c(1, 2, 3, 4, 5, 6, 7, 9, 10, 11, 12, 13, 14, 15,
                        16, 17, 30, 48))

    ## Transposed
    tmp <- TreeSummarizedExperiment::TreeSummarizedExperiment(
        assays = list(counts = t(SummarizedExperiment::assay(x, "counts"))),
        rowData = SummarizedExperiment::colData(x),
        colData = SummarizedExperiment::rowData(x),
        colTree = TreeSummarizedExperiment::rowTree(x)
    )
    out <- runDA(TSE = tmp, feature_on_row = FALSE, assay = NULL,
                 option = "glmQL", design = des, contrast = c(0, 1),
                 filter_min_count = 0, filter_min_total_count = 0,
                 filter_large_n = 0, filter_min_prop = 0,
                 normalize = TRUE, normalize_method = "TMM",
                 group_column = "group", design_terms = "group",
                 legacy = TRUE)
    expect_s4_class(out$edgeR_results, "DGELRT")
    ## Check that truly differential features are high in the results list
    pos <- sort(match(truede, rownames(edgeR::topTags(out$edgeR_results,
                                                      n = Inf))))
    expect_equal(pos, c(1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 13, 14, 15,
                        16, 17, 29, 45))
})
