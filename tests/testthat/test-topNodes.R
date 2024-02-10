test_that("topNodes works", {
    library(TreeSummarizedExperiment)
    ## Generate example data
    data(tinyTree)
    set.seed(2L)
    pv <- runif(19, 0, 1)
    pv[c(seq_len(5), 13, 14, 18)] <- runif(8, 0, 0.001)
    fc <- sample(c(-1, 1), 19, replace = TRUE) *
        runif(n = 19, min = 0.9, max = 1.1)
    fc[c(seq_len(3), 13, 14)] <- runif(n = 5, min = 0.9, max = 1.1)
    fc[c(4, 5, 18)] <- runif(n = 3, min = -1.1, max = -0.9)
    df <- data.frame(node = seq_len(19),
                     pvalue = pv,
                     foldChange = fc)
    ll <- getCand(tree = tinyTree, score_data = df,
                  node_column = "node",
                  p_column = "pvalue",
                  sign_column = "foldChange")
    cc <- evalCand(tree = tinyTree, levels = ll$candidate_list,
                   score_data = df, node_column = "node",
                   p_column = "pvalue", sign_column = "foldChange",
                   limit_rej = 0.05)

    ## Test that function returns error for misspecified input
    ## -------------------------------------------------------------------------
    expect_error(topNodes(object = 1, n = 10, sort_by = NULL,
                          sort_decreasing = FALSE,
                          sort_by_absolute = FALSE, p_value = 1),
                 "'object' must be of class 'list'")
    expect_error(topNodes(object = list(x = 3), n = 10, sort_by = NULL,
                          sort_decreasing = FALSE,
                          sort_by_absolute = FALSE, p_value = 1),
                 '"output" %in% names(object) is not TRUE', fixed = TRUE)
    expect_error(topNodes(object = list(output = 3), n = 10, sort_by = NULL,
                          sort_decreasing = FALSE,
                          sort_by_absolute = FALSE, p_value = 1),
                 "'$objectoutput' must be of class 'data.frame'", fixed = TRUE)
    expect_error(topNodes(object = list(output = data.frame(x = 3)),
                          n = 10, sort_by = NULL,
                          sort_decreasing = FALSE,
                          sort_by_absolute = FALSE, p_value = 1),
                 '"adj.p" %in% colnames(object$output) is not TRUE',
                 fixed = TRUE)
    expect_error(topNodes(object = list(output = data.frame(adj.p = 3)),
                          n = 10, sort_by = "x",
                          sort_decreasing = FALSE,
                          sort_by_absolute = FALSE, p_value = 1),
                 'sort_by %in% colnames(object$output) is not TRUE',
                 fixed = TRUE)
    expect_error(topNodes(object = cc, n = "x", sort_by = NULL,
                          sort_decreasing = FALSE,
                          sort_by_absolute = FALSE, p_value = 1),
                 "'n' must be of class 'numeric'")
    expect_error(topNodes(object = cc, n = c(1, 2), sort_by = NULL,
                          sort_decreasing = FALSE,
                          sort_by_absolute = FALSE, p_value = 1),
                 "'n' must have length 1")
    expect_error(topNodes(object = cc, n = 10, sort_by = 1,
                          sort_decreasing = FALSE,
                          sort_by_absolute = FALSE, p_value = 1),
                 "'sort_by' must be of class 'character'")
    expect_error(topNodes(object = cc, n = 10, sort_by = c("adj.p", "pvalue"),
                          sort_decreasing = FALSE,
                          sort_by_absolute = FALSE, p_value = 1),
                 "'sort_by' must have length 1")
    expect_error(topNodes(object = cc, n = 10, sort_by = NULL,
                          sort_decreasing = 1,
                          sort_by_absolute = FALSE, p_value = 1),
                 "'sort_decreasing' must be of class 'logical'")
    expect_error(topNodes(object = cc, n = 10, sort_by = NULL,
                          sort_decreasing = c(TRUE, FALSE),
                          sort_by_absolute = FALSE, p_value = 1),
                 "'sort_decreasing' must have length 1")
    expect_error(topNodes(object = cc, n = 10, sort_by = NULL,
                          sort_decreasing = FALSE,
                          sort_by_absolute = 1, p_value = 1),
                 "'sort_by_absolute' must be of class 'logical'")
    expect_error(topNodes(object = cc, n = 10, sort_by = NULL,
                          sort_decreasing = FALSE,
                          sort_by_absolute = c(TRUE, FALSE), p_value = 1),
                 "'sort_by_absolute' must have length 1")
    expect_error(topNodes(object = cc, n = 10, sort_by = NULL,
                          sort_decreasing = FALSE,
                          sort_by_absolute = FALSE, p_value = TRUE),
                 "'p_value' must be of class 'numeric'")
    expect_error(topNodes(object = cc, n = 10, sort_by = NULL,
                          sort_decreasing = FALSE,
                          sort_by_absolute = FALSE, p_value = c(0.5, 0.7)),
                 "'p_value' must have length 1")
    expect_error(topNodes(object = cc, n = 10, sort_by = NULL,
                          sort_decreasing = FALSE,
                          sort_by_absolute = FALSE, p_value = 2),
                 "'p_value' must be within [0,1]", fixed = TRUE)


    ## Test that topNodes works
    ## -------------------------------------------------------------------------
    ## No sorting, no filtering
    out <- topNodes(object = cc, n = 10, sort_by = NULL,
                    sort_decreasing = FALSE, sort_by_absolute = FALSE,
                    p_value = 1)
    expect_s3_class(out, "data.frame")
    expect_equal(nrow(out), 7)
    expect_named(out, c("node", "pvalue", "foldChange", "adj.p", "signal.node"))
    expect_equal(out$node, c(6, 7, 8, 9, 10, 13, 18))
    expect_equal(out$signal.node, rep(c(FALSE, TRUE), c(5, 2)))

    ## No sorting, filtering
    out <- topNodes(object = cc, n = 10, sort_by = NULL,
                    sort_decreasing = FALSE, sort_by_absolute = FALSE,
                    p_value = 0.35)
    expect_s3_class(out, "data.frame")
    expect_equal(nrow(out), 3)
    expect_named(out, c("node", "pvalue", "foldChange", "adj.p", "signal.node"))
    expect_equal(out$node, c(7, 13, 18))
    expect_equal(out$signal.node, rep(c(FALSE, TRUE), c(1, 2)))

    ## Sorting, filtering
    out <- topNodes(object = cc, n = 10, sort_by = "pvalue",
                    sort_decreasing = FALSE, sort_by_absolute = FALSE,
                    p_value = 0.35)
    expect_s3_class(out, "data.frame")
    expect_equal(nrow(out), 3)
    expect_named(out, c("node", "pvalue", "foldChange", "adj.p", "signal.node"))
    expect_equal(out$node, c(18, 13, 7))
    expect_equal(out$signal.node, rep(c(TRUE, FALSE), c(2, 1)))

    ## Sorting decreasing, filtering
    out <- topNodes(object = cc, n = 10, sort_by = "pvalue",
                    sort_decreasing = TRUE, sort_by_absolute = FALSE,
                    p_value = 0.35)
    expect_s3_class(out, "data.frame")
    expect_equal(nrow(out), 3)
    expect_named(out, c("node", "pvalue", "foldChange", "adj.p", "signal.node"))
    expect_equal(out$node, c(7, 13, 18))
    expect_equal(out$signal.node, rep(c(FALSE, TRUE), c(1, 2)))

    ## Top 2 only
    out <- topNodes(object = cc, n = 2, sort_by = "pvalue",
                    sort_decreasing = TRUE, sort_by_absolute = FALSE,
                    p_value = 0.35)
    expect_s3_class(out, "data.frame")
    expect_equal(nrow(out), 2)
    expect_named(out, c("node", "pvalue", "foldChange", "adj.p", "signal.node"))
    expect_equal(out$node, c(7, 13))
    expect_equal(out$signal.node, rep(c(FALSE, TRUE), c(1, 1)))

    ## Sort by absolute
    out <- topNodes(object = cc, n = 10, sort_by = "foldChange",
                    sort_decreasing = TRUE, sort_by_absolute = TRUE,
                    p_value = 1)
    expect_s3_class(out, "data.frame")
    expect_equal(nrow(out), 7)
    expect_named(out, c("node", "pvalue", "foldChange", "adj.p", "signal.node"))
    expect_equal(out$node, c(18, 8, 10, 7, 9, 13, 6))
    expect_equal(out$signal.node, rep(c(TRUE, FALSE, TRUE, FALSE),
                                      c(1, 4, 1, 1)))

    ## Sort by absolute, increasing
    out <- topNodes(object = cc, n = 10, sort_by = "foldChange",
                    sort_decreasing = FALSE, sort_by_absolute = TRUE,
                    p_value = 1)
    expect_s3_class(out, "data.frame")
    expect_equal(nrow(out), 7)
    expect_named(out, c("node", "pvalue", "foldChange", "adj.p", "signal.node"))
    expect_equal(out$node, c(6, 13, 9, 7, 10, 8, 18))
    expect_equal(out$signal.node, rep(c(FALSE, TRUE, FALSE, TRUE),
                                      c(1, 1, 4, 1)))
})
