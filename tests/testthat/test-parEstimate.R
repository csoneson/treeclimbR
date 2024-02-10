test_that("Dirichlet multinomial parameter estimation works", {
    library(TreeSummarizedExperiment)
    library(dirmult)
    set.seed(1)
    y <- matrix(rnbinom(200, size = 1, mu = 10), nrow = 10)
    colnames(y) <- paste("S", seq_len(20), sep = "")
    rownames(y) <- tinyTree$tip.label
    toy_tse <- TreeSummarizedExperiment(rowTree = tinyTree,
                                        assays = list(y))

    ## Test that function fails with invalid input
    ## -------------------------------------------------------------------------
    expect_error(.estimateA(obj = list(a = 1)),
                 "should contain pi and theta")
    expect_error(.estimateA(obj = matrix(seq_len(20), nrow = 5)),
                 "must have rownames")
    expect_error(.estimateB(obj = 1),
                 "'obj' must be of class 'list'")
    expect_error(.estimateB(obj = list(a = 1)),
                 "should contain pi and theta")
    expect_error(.estimateC(obj = 1),
                 "'obj' must be of class 'TreeSummarizedExperiment'")
    expect_error(parEstimate(obj = 1),
                 "is not TRUE")
    expect_error(parEstimate(obj = list(1), assay = c(1, 2)),
                 "length(assay) == 1 || is.null(assay) is not TRUE",
                 fixed = TRUE)
    expect_error(parEstimate(obj = toy_tse, assay = "missing"),
                 "assay %in% SummarizedExperiment::assayNames(obj) is not TRUE",
                 fixed = TRUE)

    ## Test that function works with correct input
    ## -------------------------------------------------------------------------
    set.seed(1L)
    res <- parEstimate(obj = toy_tse, assay = 1)
    set.seed(1L)
    res2 <- parEstimate(obj = toy_tse, assay = NULL)
    expect_equal(res, res2)
    expect_s4_class(res, "TreeSummarizedExperiment")
    expect_equal(dim(res), dim(toy_tse))
    expect_type(S4Vectors::metadata(res)$assays.par, "list")
    p <- S4Vectors::metadata(res)$assays.par
    expect_length(p, 2)
    expect_named(p, c("pi", "theta"))
    expect_type(p$pi, "double")
    expect_length(p$pi, nrow(res))
    expect_type(p$theta, "double")
    expect_length(p$theta, 1)
    expect_equal(p$theta, 0.1016139, tolerance = 1e-6)
    expect_equal(p$pi[c(1, 4, 9)], c(0.10811579, 0.09701295, 0.10676714),
                 tolerance = 1e-6, ignore_attr = TRUE)
    expect_named(p$pi, c("t2", "t7", "t6", "t9", "t4", "t8", "t10", "t1",
                         "t5", "t3"))
    expect_named(p$pi, TreeSummarizedExperiment::rowTree(toy_tse)$tip.label)

    ## Test with simulated input with known parameters
    ## -------------------------------------------------------------------------
    set.seed(1L)
    simdatL <- dirmult::simPop(J = 10, K = 5, n = 1000,
                               pi = c(0.05, 0.1, 0.15, 0.2, 0.5), theta = 0.25)
    simdat <- simdatL$data
    colnames(simdat) <- paste0("X", seq_len(ncol(simdat)))
    res <- parEstimate(t(simdat))
    expect_type(res, "list")
    expect_named(res, c("pi", "theta"))
    expect_length(res$pi, 5)
    expect_true(max(abs(res$pi - c(0.05, 0.1, 0.15, 0.2, 0.5))) < 0.08)

    ## Test helper functions explicitly
    ## -------------------------------------------------------------------------
    expect_equal(.estimateA(obj = list(pi = c(0.4, 0.6), theta = 0.2)),
                 list(pi = c(0.4, 0.6), theta = 0.2))
    expect_equal(.estimateB(obj = list(pi = c(0.4, 0.6), theta = 0.2)),
                 list(pi = c(0.4, 0.6), theta = 0.2))
    expect_equal(parEstimate(obj = list(pi = c(0.4, 0.6), theta = 0.2)),
                 list(pi = c(0.4, 0.6), theta = 0.2))
})
