test_that("infoCand works", {
    library(TreeSummarizedExperiment)
    library(ggtree)
    ## Generate example data
    data(tinyTree)
    set.seed(2L)
    pv <- runif(19, 0, 1)
    pv[c(seq_len(5), 13, 14, 18)] <- runif(8, 0, 0.001)
    fc <- sample(c(-1, 1), 19, replace = TRUE)
    fc[c(seq_len(3), 13, 14)] <- 1
    fc[c(4, 5, 18)] <- -1
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
    expect_error(infoCand(object = 1),
                 "'object' must be of class 'list'")
    obj <- cc
    obj$level_info <- NULL
    expect_error(infoCand(object = obj),
                 "object needs to have a 'level_info' slot")

    ## Test that infoCand works
    ## -------------------------------------------------------------------------
    out <- infoCand(object = cc)
    expect_s3_class(out, "data.frame")
    expect_equal(nrow(out), 25L)
    expect_equal(ncol(out), 9L)
    expect_named(out, c("t", "upper_t", "is_valid", "method", "limit_rej",
                        "level_name", "best", "rej_leaf", "rej_node"))
    expect_equal(out$t, c(seq(0, 0.05, by = 0.01), seq(0.1, 1, by = 0.05)))
    expect_equal(sum(out$best), 6)
    expect_equal(out$rej_leaf, rep(5, 25))
    expect_equal(out$rej_node, c(5, rep(2, 24)))
    expect_equal(out$is_valid, rep(c(TRUE, FALSE), c(7, 18)))
})
