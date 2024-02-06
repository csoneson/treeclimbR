test_that("isConnect works", {
    library(ggtree)
    data(tinyTree)

    ## Check that the function fails with wrongly formatted input
    ## -------------------------------------------------------------------------
    expect_error(isConnect(tree = 1, node_a = c(5, 12),
                           node_b = c(18, 19), connect = "any"),
                 "'tree' must be of class 'phylo'")
    expect_error(isConnect(tree = tinyTree, node_a = TRUE, node_b = 18,
                           connect = "any"),
                 "'node_a' should be either a character vector or")
    expect_error(isConnect(tree = tinyTree, node_a = c(), node_b = 18,
                           connect = "any"),
                 "'node_a' should be either a character vector or")
    expect_error(isConnect(tree = tinyTree, node_a = 18, node_b = TRUE,
                           connect = "any"),
                 "'node_b' should be either a character vector or")
    expect_error(isConnect(tree = tinyTree, node_a = 18, node_b = c(),
                           connect = "any"),
                 "'node_b' should be either a character vector or")
    expect_error(isConnect(tree = tinyTree, node_a = c(12, 17),
                           node_b = 18, connect = 1),
                 "'connect' must be of class 'character'")
    expect_error(isConnect(tree = tinyTree, node_a = c(12, 17),
                           node_b = 18, connect = c("any", "direct")),
                 "'connect' must have length 1")
    expect_error(isConnect(tree = tinyTree, node_a = c(12, 17),
                           node_b = 18, connect = "missing"),
                 "All values in 'connect' must be one of")

    ## Check that the function works with correct input
    ## -------------------------------------------------------------------------
    ## Choose pairs to cover all situations
    node_a <- c(5, 5, 18, 5, 15, 10, 18, 18, 5)
    node_b <- c(18, 4, 17, 17, 17, 11, 14, 2, 15)

    expect_equal(isConnect(tree = tinyTree, node_a = node_a, node_b = node_b,
                           connect = "direct"),
                 c(TRUE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE))
    expect_equal(isConnect(tree = tinyTree, node_a = node_a, node_b = node_b,
                           connect = "indirect"),
                 c(FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE))
    expect_equal(isConnect(tree = tinyTree, node_a = node_a, node_b = node_b,
                           connect = "any"),
                 c(TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE))

    ## Check that it works with a single value as well
    for (i in seq_along(node_a)) {
        expect_equal(isConnect(tree = tinyTree, node_a = node_a[i],
                               node_b = node_b[i], connect = "direct"),
                     c(TRUE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE,
                       FALSE, FALSE)[i])
        expect_equal(isConnect(tree = tinyTree, node_a = node_a[i],
                               node_b = node_b[i], connect = "indirect"),
                     c(FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE,
                       FALSE, TRUE)[i])
        expect_equal(isConnect(tree = tinyTree, node_a = node_a[i],
                               node_b = node_b[i], connect = "any"),
                     c(TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE,
                       FALSE, TRUE)[i])
    }

    ## Same with node labels
    node_a <- c("t4", "t4", "Node_18", "t4", "Node_15", "t3", "Node_18",
                "Node_18", "t4")
    node_b <- c("Node_18", "t9", "Node_17", "Node_17", "Node_17",
                "Node_11", "Node_14", "t7", "Node_15")

    expect_warning(
        expect_warning(
            expect_equal(isConnect(tree = tinyTree, node_a = node_a,
                                   node_b = node_b, connect = "direct"),
                         c(TRUE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE,
                           FALSE, FALSE),
                         ignore_attr = TRUE),
            "Multiple nodes are found to have the same label"),
        "Multiple nodes are found to have the same label")
    expect_warning(
        expect_warning(
            expect_equal(isConnect(tree = tinyTree, node_a = node_a,
                                   node_b = node_b, connect = "indirect"),
                         c(FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE,
                           FALSE, TRUE),
                         ignore_attr = TRUE),
            "Multiple nodes are found to have the same label"),
        "Multiple nodes are found to have the same label")
    expect_warning(
        expect_warning(
            expect_equal(isConnect(tree = tinyTree, node_a = node_a,
                                   node_b = node_b, connect = "any"),
                         c(TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE,
                           FALSE, TRUE),
                         ignore_attr = TRUE),
            "Multiple nodes are found to have the same label"),
        "Multiple nodes are found to have the same label")

    ## Same with node aliases
    node_a <- paste0("alias_", c(5, 5, 18, 5, 15, 10, 18, 18, 5))
    node_b <- c(18, 4, 17, 17, 17, 11, 14, 2, 15)

    expect_equal(isConnect(tree = tinyTree, node_a = node_a, node_b = node_b,
                           connect = "direct"),
                 c(TRUE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE),
                 ignore_attr = TRUE)
    expect_equal(isConnect(tree = tinyTree, node_a = node_a, node_b = node_b,
                           connect = "indirect"),
                 c(FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE),
                 ignore_attr = TRUE)
    expect_equal(isConnect(tree = tinyTree, node_a = node_a, node_b = node_b,
                           connect = "any"),
                 c(TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE),
                 ignore_attr = TRUE)

    ## Same with node aliases for node_b
    node_a <- c(5, 5, 18, 5, 15, 10, 18, 18, 5)
    node_b <- paste0("alias_", c(18, 4, 17, 17, 17, 11, 14, 2, 15))

    expect_equal(isConnect(tree = tinyTree, node_a = node_a, node_b = node_b,
                           connect = "direct"),
                 c(TRUE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE))
    expect_equal(isConnect(tree = tinyTree, node_a = node_a, node_b = node_b,
                           connect = "indirect"),
                 c(FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE))
    expect_equal(isConnect(tree = tinyTree, node_a = node_a, node_b = node_b,
                           connect = "any"),
                 c(TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE))

})
