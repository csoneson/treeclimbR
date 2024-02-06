test_that("findChild works", {
    library(ggtree)
    data(tinyTree)

    ## Check that the function fails with wrongly formatted input
    ## -------------------------------------------------------------------------
    expect_error(findChild(tree = 1, node = c(12, 17), use.alias = FALSE),
                 "'tree' must be of class 'phylo'")
    expect_error(findChild(tree = tinyTree, node = TRUE, use.alias = FALSE),
                 "'node' should be either a character vector or")
    expect_error(findChild(tree = tinyTree, node = c(12, 17), use.alias = 1),
                 "'use.alias' must be of class 'logical'")
    expect_error(findChild(tree = tinyTree, node = 120, use.alias = FALSE),
                 "Node 120 can't be found in tinyTree")
    expect_error(findChild(tree = tinyTree, node = "Node_120",
                           use.alias = FALSE),
                 "1 nodes mismatch with the tree")

    ## Check that the function works with the correct input
    ## -------------------------------------------------------------------------
    expect_equal(findChild(tree = tinyTree, node = c(12, 17),
                           use.alias = FALSE),
                 list(Node_12 = c(13L, 15L), Node_17 = c(6L, 18L)))
    expect_equal(findChild(tree = tinyTree, node = c(12, 17),
                           use.alias = TRUE),
                 list(alias_12 = c(13L, 15L), alias_17 = c(6L, 18L)))
    expect_equal(findChild(tree = tinyTree, node = c("Node_12", "Node_17"),
                           use.alias = FALSE),
                 list(Node_12 = c(13L, 15L), Node_17 = c(6L, 18L)))
    expect_equal(findChild(tree = tinyTree, node = c("Node_12", "Node_17"),
                           use.alias = TRUE),
                 list(alias_12 = c(13L, 15L), alias_17 = c(6L, 18L)))

    expect_equal(findChild(tree = tinyTree, node = c(11, 10),
                           use.alias = FALSE),
                 list(Node_11 = c(10L, 12L), t3 = integer(0)))
    expect_equal(findChild(tree = tinyTree, node = c(11, 10),
                           use.alias = TRUE),
                 list(alias_11 = c(10L, 12L), alias_10 = integer(0)))
    expect_equal(findChild(tree = tinyTree, node = c("Node_11", "t3"),
                           use.alias = FALSE),
                 list(Node_11 = c(10L, 12L), t3 = integer(0)))
    expect_equal(findChild(tree = tinyTree, node = c("Node_11", "t3"),
                           use.alias = TRUE),
                 list(alias_11 = c(10L, 12L), alias_10 = integer(0)))
})
