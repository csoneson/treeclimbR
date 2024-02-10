test_that("fdr and tpr calculations works", {
    suppressPackageStartupMessages({
        library(ggtree)
        library(TreeSummarizedExperiment)
    })
    data("tinyTree")

    ## Check that functions fail if input is incorrect
    ## -------------------------------------------------------------------------
    expect_error(fdr(tree = 1, truth = c(1, 2),
                     found = c(1, 2), only.leaf = TRUE),
                 "'tree' must be of class 'phylo'")
    expect_error(tpr(tree = 1, truth = c(1, 2),
                     found = c(1, 2), only.leaf = TRUE),
                 "'tree' must be of class 'phylo'")

    expect_error(fdr(tree = tinyTree, truth = TRUE,
                     found = c(1, 2), only.leaf = TRUE),
                 "'truth' should be either a character vector or a numeric")
    expect_error(tpr(tree = tinyTree, truth = FALSE,
                     found = c(1, 2), only.leaf = TRUE),
                 "'truth' should be either a character vector or a numeric")

    expect_error(fdr(tree = tinyTree, truth = c(1, 2),
                     found = TRUE, only.leaf = TRUE),
                 "'found' should be either a character vector or a numeric")
    expect_error(tpr(tree = tinyTree, truth = c(1, 2),
                     found = FALSE, only.leaf = TRUE),
                 "'found' should be either a character vector or a numeric")

    expect_error(fdr(tree = tinyTree, truth = c(1, 2),
                     found = c(1, 2), only.leaf = 1),
                 "'only.leaf' must be of class 'logical'")
    expect_error(tpr(tree = tinyTree, truth = c(1, 2),
                     found = c(1, 2), only.leaf = 2),
                 "'only.leaf' must be of class 'logical'")


    ## Check that output is correct if input is valid
    ## -------------------------------------------------------------------------
    expect_equal(fdr(tree = tinyTree, truth = c(16, 13),
                     found = c(13), only.leaf = TRUE),
                 c(fdr = 0))
    expect_equal(tpr(tree = tinyTree, truth = c(16, 13),
                     found = c(13), only.leaf = TRUE),
                 c(tpr = 3/8))
    expect_equal(fdr(tree = tinyTree, truth = c(16, 13),
                     found = c(13), only.leaf = FALSE),
                 c(fdr = 0))
    expect_equal(tpr(tree = tinyTree, truth = c(16, 13),
                     found = c(13), only.leaf = FALSE),
                 c(tpr = 5/14))

    expect_equal(fdr(tree = tinyTree, truth = c(16, 13),
                     found = c(13, 16), only.leaf = TRUE),
                 c(fdr = 0))
    expect_equal(tpr(tree = tinyTree, truth = c(16, 13),
                     found = c(13, 16), only.leaf = TRUE),
                 c(tpr = 1))
    expect_equal(fdr(tree = tinyTree, truth = c(16, 13),
                     found = c(13, 16), only.leaf = FALSE),
                 c(fdr = 0))
    expect_equal(tpr(tree = tinyTree, truth = c(16, 13),
                     found = c(13, 16), only.leaf = FALSE),
                 c(tpr = 1))

    expect_equal(fdr(tree = tinyTree, truth = c(16, 13),
                     found = NULL, only.leaf = TRUE),
                 c(fdr = 0))
    expect_equal(tpr(tree = tinyTree, truth = c(16, 13),
                     found = NULL, only.leaf = TRUE),
                 c(tpr = 0))
    expect_equal(fdr(tree = tinyTree, truth = c(16, 13),
                     found = NULL, only.leaf = FALSE),
                 c(fdr = 0))
    expect_equal(tpr(tree = tinyTree, truth = c(16, 13),
                     found = NULL, only.leaf = FALSE),
                 c(tpr = 0))

    expect_equal(fdr(tree = tinyTree, truth = NULL,
                     found = c(13), only.leaf = TRUE),
                 c(fdr = 1))
    expect_equal(tpr(tree = tinyTree, truth = NULL,
                     found = c(13), only.leaf = TRUE),
                 c(tpr = 1))
    expect_equal(fdr(tree = tinyTree, truth = NULL,
                     found = c(13), only.leaf = FALSE),
                 c(fdr = 1))
    expect_equal(tpr(tree = tinyTree, truth = NULL,
                     found = c(13), only.leaf = FALSE),
                 c(tpr = 1))

    expect_equal(fdr(tree = tinyTree, truth = c(16, 13),
                     found = c(14, 15), only.leaf = TRUE),
                 c(fdr = 1/8))
    expect_equal(tpr(tree = tinyTree, truth = c(16, 13),
                     found = c(14, 15), only.leaf = TRUE),
                 c(tpr = 7/8))
    expect_equal(fdr(tree = tinyTree, truth = c(16, 13),
                     found = c(14, 15), only.leaf = FALSE),
                 c(fdr = 2/14))
    expect_equal(tpr(tree = tinyTree, truth = c(16, 13),
                     found = c(14, 15), only.leaf = FALSE),
                 c(tpr = 12/14))

    expect_equal(fdr(tree = tinyTree, truth = c(16, 13),
                     found = convertNode(tinyTree, c(14, 15)), only.leaf = TRUE),
                 c(fdr = 1/8))
    expect_equal(tpr(tree = tinyTree, truth = c(16, 13),
                     found = convertNode(tinyTree, c(14, 15)), only.leaf = TRUE),
                 c(tpr = 7/8))
    expect_equal(fdr(tree = tinyTree, truth = c(16, 13),
                     found = convertNode(tinyTree, c(14, 15)), only.leaf = FALSE),
                 c(fdr = 2/14))
    expect_equal(tpr(tree = tinyTree, truth = c(16, 13),
                     found = convertNode(tinyTree, c(14, 15)), only.leaf = FALSE),
                 c(tpr = 12/14))

    expect_equal(fdr(tree = tinyTree, truth = convertNode(tinyTree, c(16, 13)),
                     found = convertNode(tinyTree, c(14, 15)), only.leaf = TRUE),
                 c(fdr = 1/8))
    expect_equal(tpr(tree = tinyTree, truth = convertNode(tinyTree, c(16, 13)),
                     found = convertNode(tinyTree, c(14, 15)), only.leaf = TRUE),
                 c(tpr = 7/8))
    expect_equal(fdr(tree = tinyTree, truth = convertNode(tinyTree, c(16, 13)),
                     found = convertNode(tinyTree, c(14, 15)), only.leaf = FALSE),
                 c(fdr = 2/14))
    expect_equal(tpr(tree = tinyTree, truth = convertNode(tinyTree, c(16, 13)),
                     found = convertNode(tinyTree, c(14, 15)), only.leaf = FALSE),
                 c(tpr = 12/14))

    expect_equal(fdr(tree = tinyTree, truth = convertNode(tinyTree, c(16, 13)),
                     found = convertNode(tinyTree, c(17, 3, 9)), only.leaf = TRUE),
                 c(fdr = 1/5))
    expect_equal(tpr(tree = tinyTree, truth = convertNode(tinyTree, c(16, 13)),
                     found = convertNode(tinyTree, c(17, 3, 9)), only.leaf = TRUE),
                 c(tpr = 4/8))
    expect_equal(fdr(tree = tinyTree, truth = convertNode(tinyTree, c(16, 13)),
                     found = convertNode(tinyTree, c(17, 3, 9)), only.leaf = FALSE),
                 c(fdr = 1/7))
    expect_equal(tpr(tree = tinyTree, truth = convertNode(tinyTree, c(16, 13)),
                     found = convertNode(tinyTree, c(17, 3, 9)), only.leaf = FALSE),
                 c(tpr = 6/14))

    expect_equal(fdr(tree = tinyTree, truth = convertNode(tinyTree, c(16, 13)),
                     found = c(17, 3, 9), only.leaf = TRUE),
                 c(fdr = 1/5))
    expect_equal(tpr(tree = tinyTree, truth = convertNode(tinyTree, c(16, 13)),
                     found = c(17, 3, 9), only.leaf = TRUE),
                 c(tpr = 4/8))
    expect_equal(fdr(tree = tinyTree, truth = convertNode(tinyTree, c(16, 13)),
                     found = c(17, 3, 9), only.leaf = FALSE),
                 c(fdr = 1/7))
    expect_equal(tpr(tree = tinyTree, truth = convertNode(tinyTree, c(16, 13)),
                     found = c(17, 3, 9), only.leaf = FALSE),
                 c(tpr = 6/14))

})
