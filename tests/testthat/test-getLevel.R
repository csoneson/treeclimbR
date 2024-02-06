test_that("getLevel works", {
    ## Generate some data
    library(TreeSummarizedExperiment)
    library(ggtree)
    data(tinyTree)
    set.seed(1)
    pv <- runif(19, min = 0.09, max = 0.11)
    pv[c(16, 13, 17)] <- c(0.01, 0.04, 0.005)
    out <- data.frame(node = 1:19, pvalue = pv)

    ## Check that function returns errors for invalid input
    ## -------------------------------------------------------------------------
    expect_error(getLevel(tree = 1, score_data = out, drop = pvalue > 0.05,
                          score_column = "pvalue", node_column = "node",
                          get_max = FALSE, parent_first = TRUE,
                          message = FALSE),
                 "'tree' must be of class 'phylo'")
    expect_error(getLevel(tree = tinyTree, score_data = 1,
                          drop = pvalue > 0.05,
                          score_column = "pvalue", node_column = "node",
                          get_max = FALSE, parent_first = TRUE,
                          message = FALSE),
                 "score_data should be a data.frame")
    expect_error(getLevel(tree = tinyTree, score_data = out,
                          drop = pvalue > 0.05,
                          score_column = 1, node_column = "node",
                          get_max = FALSE, parent_first = TRUE,
                          message = FALSE),
                 "'score_column' must be of class 'character'")
    expect_error(getLevel(tree = tinyTree, score_data = out,
                          drop = pvalue > 0.05,
                          score_column = c("node", "pvalue"),
                          node_column = "node",
                          get_max = FALSE, parent_first = TRUE,
                          message = FALSE),
                 "'score_column' must have length 1")
    expect_error(getLevel(tree = tinyTree, score_data = out,
                          drop = pvalue > 0.05,
                          score_column = "missing",
                          node_column = "node",
                          get_max = FALSE, parent_first = TRUE,
                          message = FALSE),
                 "All values in 'score_column' must be one of")
    expect_error(getLevel(tree = tinyTree, score_data = out,
                          drop = pvalue > 0.05,
                          score_column = "pvalue", node_column = 1,
                          get_max = FALSE, parent_first = TRUE,
                          message = FALSE),
                 "'node_column' must be of class 'character'")
    expect_error(getLevel(tree = tinyTree, score_data = out,
                          drop = pvalue > 0.05,
                          score_column = "pvalue",
                          node_column = c("node", "pvalue"),
                          get_max = FALSE, parent_first = TRUE,
                          message = FALSE),
                 "'node_column' must have length 1")
    expect_error(getLevel(tree = tinyTree, score_data = out,
                          drop = pvalue > 0.05,
                          score_column = "pvalue",
                          node_column = "missing",
                          get_max = FALSE, parent_first = TRUE,
                          message = FALSE),
                 "All values in 'node_column' must be one of")
    expect_error(getLevel(tree = tinyTree, score_data = out,
                          drop = pvalue > 0.05,
                          score_column = "pvalue", node_column = "node",
                          get_max = 1, parent_first = TRUE,
                          message = FALSE),
                 "'get_max' must be of class 'logical'")
    expect_error(getLevel(tree = tinyTree, score_data = out,
                          drop = pvalue > 0.05,
                          score_column = "pvalue", node_column = "node",
                          get_max = c(TRUE, FALSE), parent_first = TRUE,
                          message = FALSE),
                 "'get_max' must have length 1")
    expect_error(getLevel(tree = tinyTree, score_data = out,
                          drop = pvalue > 0.05,
                          score_column = "pvalue", node_column = "node",
                          get_max = FALSE, parent_first = 1,
                          message = FALSE),
                 "'parent_first' must be of class 'logical'")
    expect_error(getLevel(tree = tinyTree, score_data = out,
                          drop = pvalue > 0.05,
                          score_column = "pvalue", node_column = "node",
                          get_max = FALSE, parent_first = c(TRUE, FALSE),
                          message = FALSE),
                 "'parent_first' must have length 1")
    expect_error(getLevel(tree = tinyTree, score_data = out,
                          drop = pvalue > 0.05,
                          score_column = "pvalue", node_column = "node",
                          get_max = FALSE, parent_first = TRUE,
                          message = 1),
                 "'message' must be of class 'logical'")
    expect_error(getLevel(tree = tinyTree, score_data = out,
                          drop = pvalue > 0.05,
                          score_column = "pvalue", node_column = "node",
                          get_max = FALSE, parent_first = TRUE,
                          message = c(TRUE, FALSE)),
                 "'message' must have length 1")

    tmp <- out
    tmp$keep <- 1
    expect_error(getLevel(tree = tinyTree, score_data = tmp,
                          drop = pvalue > 0.05,
                          score_column = "pvalue", node_column = "node",
                          get_max = FALSE, parent_first = TRUE,
                          message = TRUE),
                 "The result will be output in the 'keep' column")

    expect_error(getLevel(tree = tinyTree, score_data = out,
                          drop = "pvalue > 0.05",
                          score_column = "pvalue", node_column = "node",
                          get_max = FALSE, parent_first = TRUE,
                          message = TRUE),
                 "'drop' must be or evaluate to logical")

    ## TODO: Add tests where some p-values are NA

    ## Check that function works as expected for valid input
    ## -------------------------------------------------------------------------
    final <- getLevel(tree = tinyTree, score_data = out,
                      drop =  pvalue > 0.05, score_column = "pvalue",
                      node_column = "node", get_max = FALSE,
                      parent_first = TRUE, message = FALSE)
    expect_s3_class(final, "data.frame")
    expect_equal(nrow(final), nrow(out))
    expect_equal(final$node, out$node)
    expect_equal(final$pvalue, out$pvalue)
    expect_equal(final$node[final$keep], c(13, 17))

    ## parent_first = FALSE
    final <- getLevel(tree = tinyTree, score_data = DataFrame(out),
                      drop =  pvalue > 0.05, score_column = "pvalue",
                      node_column = "node", get_max = FALSE,
                      parent_first = FALSE, message = FALSE)
    expect_s3_class(final, "data.frame")
    expect_equal(nrow(final), nrow(out))
    expect_equal(final$node, out$node)
    expect_equal(final$pvalue, out$pvalue)
    expect_equal(final$node[final$keep], c(13, 17))

    ## Make node 16 lower than its descendant (node 17)
    out2 <- out
    out2$pvalue[out2$node == 16] <- 0.001
    final <- getLevel(tree = tinyTree, score_data = out2,
                      drop =  pvalue > 0.05, score_column = "pvalue",
                      node_column = "node", get_max = FALSE,
                      parent_first = TRUE, message = FALSE)
    expect_s3_class(final, "data.frame")
    expect_equal(nrow(final), nrow(out2))
    expect_equal(final$node, out2$node)
    expect_equal(final$pvalue, out2$pvalue)
    expect_equal(final$node[final$keep], c(13, 16))

    ## Search for the highest value instead
    final <- getLevel(tree = tinyTree, score_data = out,
                      drop =  pvalue > 0.05, score_column = "pvalue",
                      node_column = "node", get_max = TRUE,
                      parent_first = TRUE, message = FALSE)
    expect_s3_class(final, "data.frame")
    expect_equal(nrow(final), nrow(out))
    expect_equal(final$node, out$node)
    expect_equal(final$pvalue, out$pvalue)
    expect_equal(final$node[final$keep], integer(0))

    ## Search for the highest value, don't filter
    final <- getLevel(tree = tinyTree, score_data = out,
                      drop =  pvalue > 1, score_column = "pvalue",
                      node_column = "node", get_max = TRUE,
                      parent_first = TRUE, message = FALSE)
    expect_s3_class(final, "data.frame")
    expect_equal(nrow(final), nrow(out))
    expect_equal(final$node, out$node)
    expect_equal(final$pvalue, out$pvalue)
    expect_equal(final$node[final$keep], c(1, 2, 3, 6, 7, 8, 9, 10, 18))

    ## Search for the highest value, don't filter, parent_first = FALSE
    final <- getLevel(tree = tinyTree, score_data = out,
                      drop =  pvalue > 1, score_column = "pvalue",
                      node_column = "node", get_max = TRUE,
                      parent_first = FALSE, message = TRUE)
    expect_s3_class(final, "data.frame")
    expect_equal(nrow(final), nrow(out))
    expect_equal(final$node, out$node)
    expect_equal(final$pvalue, out$pvalue)
    expect_equal(final$node[final$keep], c(1, 2, 3, 6, 7, 8, 9, 10, 18))

    ## Search for the highest value, don't filter, parent_first = FALSE
    final <- getLevel(tree = tinyTree, score_data = out,
                      score_column = "pvalue",
                      node_column = "node", get_max = TRUE,
                      parent_first = FALSE, message = TRUE)
    expect_s3_class(final, "data.frame")
    expect_equal(nrow(final), nrow(out))
    expect_equal(final$node, out$node)
    expect_equal(final$pvalue, out$pvalue)
    expect_equal(final$node[final$keep], c(1, 2, 3, 6, 7, 8, 9, 10, 18))
})
