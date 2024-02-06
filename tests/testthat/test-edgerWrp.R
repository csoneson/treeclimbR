test_that("edgerWrp works", {
    x <- readRDS(system.file("extdata/da_sim_100_30_18de.rds",
                             package = "treeclimbR"))
    des <- model.matrix(~ group, data = SummarizedExperiment::colData(x))

    ## Check that function fails with incorrect input
    ## -------------------------------------------------------------------------
    expect_error(edgerWrp(count = 1, lib_size = NULL, option = "glm",
                          design = des, contrast = c(0, 1), normalize = TRUE,
                          normalize_method = "TMM"),
                 "'count' must be of class 'matrix'")
    expect_error(edgerWrp(count = SummarizedExperiment::assay(x),
                          lib_size = c("a", "b"), option = "glm",
                          design = des, contrast = c(0, 1), normalize = TRUE,
                          normalize_method = "TMM"),
                 "'lib_size' must be of class 'numeric'")
    expect_error(edgerWrp(count = SummarizedExperiment::assay(x),
                          lib_size = c(100, 200), option = "glm",
                          design = des, contrast = c(0, 1), normalize = TRUE,
                          normalize_method = "TMM"),
                 "'lib_size' must have length 30")
    expect_error(edgerWrp(count = SummarizedExperiment::assay(x),
                          lib_size = NULL, option = 1,
                          design = des, contrast = c(0, 1), normalize = TRUE,
                          normalize_method = "TMM"),
                 "'arg' must be NULL or a character vector")
    expect_error(edgerWrp(count = SummarizedExperiment::assay(x),
                          lib_size = NULL, option = "missing",
                          design = des, contrast = c(0, 1), normalize = TRUE,
                          normalize_method = "TMM"),
                 "'arg' should be one of")
    expect_error(edgerWrp(count = SummarizedExperiment::assay(x),
                          lib_size = NULL, option = "glm",
                          design = NULL, contrast = c(0, 1), normalize = TRUE,
                          normalize_method = "TMM"),
                 "'design' must not be NULL")
    expect_error(edgerWrp(count = SummarizedExperiment::assay(x),
                          lib_size = NULL, option = "glm",
                          design = "NULL", contrast = c(0, 1), normalize = TRUE,
                          normalize_method = "TMM"),
                 "'design' must be of class 'matrix'")
    expect_error(edgerWrp(count = SummarizedExperiment::assay(x),
                          lib_size = NULL, option = "glm",
                          design = des, contrast = "2", normalize = TRUE,
                          normalize_method = "TMM"),
                 "'contrast' must be of class 'numeric'")
    expect_error(edgerWrp(count = SummarizedExperiment::assay(x),
                          lib_size = NULL, option = "glm",
                          design = des, contrast = c(1, 2, 3), normalize = TRUE,
                          normalize_method = "TMM"),
                 "'contrast' must have length 2")
    expect_error(edgerWrp(count = SummarizedExperiment::assay(x),
                          lib_size = NULL, option = "glm",
                          design = des, contrast = c(0, 1), normalize = 1,
                          normalize_method = "TMM"),
                 "'normalize' must be of class 'logical'")
    expect_error(edgerWrp(count = SummarizedExperiment::assay(x),
                          lib_size = NULL, option = "glm",
                          design = des, contrast = c(0, 1),
                          normalize = c(TRUE, FALSE),
                          normalize_method = "TMM"),
                 "'normalize' must have length 1")
    expect_error(edgerWrp(count = SummarizedExperiment::assay(x),
                          lib_size = NULL, option = "glm",
                          design = des, contrast = c(0, 1), normalize = TRUE,
                          normalize_method = 1),
                 "'normalize_method' must be of class 'character'")
    expect_error(edgerWrp(count = SummarizedExperiment::assay(x),
                          lib_size = NULL, option = "glm",
                          design = des, contrast = c(0, 1), normalize = TRUE,
                          normalize_method = "missing"),
                 "'arg' should be one of")

    ## Check that the function works with correct inputs
    ## -------------------------------------------------------------------------
    mat <- SummarizedExperiment::assay(x)
    truede <- rownames(x)[SummarizedExperiment::rowData(x)$Signal]

    out <- edgerWrp(count = mat, lib_size = colSums(mat),
                    option = "glm", design = des, contrast = c(0, 1),
                    normalize = TRUE, normalize_method = "TMM")
    ## Should get identical results without specifying lib_size and contrast
    ## (only the 'comparison' string will be different)
    out2 <- edgerWrp(count = mat, lib_size = NULL,
                     option = "glm", design = des, contrast = NULL,
                     normalize = TRUE, normalize_method = "TMM")
    expect_s4_class(out, "DGELRT")
    expect_identical(out$table, out2$table)
    expect_identical(out$AveLogCPM, out2$AveLogCPM)
    expect_identical(out$dispersion, out2$dispersion)
    expect_identical(out$coefficients, out2$coefficients)
    ## Check that truly differential features are high in the results list
    pos <- sort(match(truede, rownames(edgeR::topTags(out, n = Inf))))
    expect_equal(pos, c(1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 12, 13, 14, 15,
                        16, 17, 30, 47))

    ## Without normalization factors
    out <- edgerWrp(count = mat, lib_size = colSums(mat),
                    option = "glm", design = des, contrast = c(0, 1),
                    normalize = FALSE, normalize_method = "TMM")
    expect_s4_class(out, "DGELRT")
    ## Check that truly differential features are high in the results list
    pos <- sort(match(truede, rownames(edgeR::topTags(out, n = Inf))))
    expect_equal(pos, c(2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 15, 17,
                        18, 21, 40, 62))

    ## With the QL framework
    out <- edgerWrp(count = mat, lib_size = colSums(mat),
                    option = "glmQL", design = des, contrast = c(0, 1),
                    normalize = TRUE, normalize_method = "TMM")
    expect_s4_class(out, "DGELRT")
    ## Check that truly differential features are high in the results list
    pos <- sort(match(truede, rownames(edgeR::topTags(out, n = Inf))))
    expect_equal(pos, c(1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 13, 14, 15,
                        16, 17, 29, 45))
})
