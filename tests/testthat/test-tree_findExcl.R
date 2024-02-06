test_that("findExcl works", {
    library(ggtree)
    data(tinyTree)

    ## Check that the function fails with wrongly formatted input
    ## -------------------------------------------------------------------------
    expect_error(findExcl(tree = 1, node = c(12, 17), use.alias = FALSE),
                 "'tree' must be of class 'phylo'")
    expect_error(findExcl(tree = tinyTree, node = TRUE, use.alias = FALSE),
                 "'node' should be either a character vector or")
    expect_error(findExcl(tree = tinyTree, node = c(), use.alias = FALSE),
                 "'node' should be either a character vector or")
    expect_error(findExcl(tree = tinyTree, node = c(12, 17), use.alias = 1),
                 "'use.alias' must be of class 'logical'")
    expect_error(findExcl(tree = tinyTree, node = 120, use.alias = FALSE),
                 "Node 120 can't be found in the tree")
    expect_error(findExcl(tree = tinyTree, node = "Node_120",
                          use.alias = FALSE),
                 "1 nodes mismatch with the tree")

    ## Check that the function works with the correct input
    ## -------------------------------------------------------------------------
    expect_equal(findExcl(tree = tinyTree, node = c(12, 17),
                          use.alias = FALSE),
                 c(t3 = 10))
    expect_equal(findExcl(tree = tinyTree, node = 11,
                          use.alias = FALSE),
                 integer(0), ignore_attr = TRUE)
    expect_equal(findExcl(tree = tinyTree, node = c(18, 19),
                          use.alias = FALSE),
                 c(Node_13 = 13, t8 = 6, t5 = 9, t3 = 10), ignore_attr = TRUE)

    expect_equal(findExcl(tree = tinyTree, node = c("Node_12", "Node_17"),
                          use.alias = FALSE),
                 c(t3 = 10))
    expect_equal(findExcl(tree = tinyTree, node = "Node_11",
                          use.alias = FALSE),
                 integer(0), ignore_attr = TRUE)
    expect_equal(findExcl(tree = tinyTree, node = c("Node_18", "Node_19"),
                          use.alias = FALSE),
                 c(Node_13 = 13, t8 = 6, t5 = 9, t3 = 10), ignore_attr = TRUE)

    expect_equal(findExcl(tree = tinyTree, node = c(12, 17),
                          use.alias = TRUE),
                 c(alias_10 = 10))
    expect_equal(findExcl(tree = tinyTree, node = 11,
                          use.alias = TRUE),
                 integer(0), ignore_attr = TRUE)
    expect_equal(findExcl(tree = tinyTree, node = c(18, 19),
                          use.alias = TRUE),
                 c(alias_13 = 13, alias_6 = 6, alias_9 = 9, alias_10 = 10),
                 ignore_attr = TRUE)
})
