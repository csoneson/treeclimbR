test_that("getCand works", {
    ## Generate some data
    library(TreeSummarizedExperiment)
    data(tinyTree)
    set.seed(1)
    pv <- runif(19, 0, 1)
    pv[c(seq_len(5), 13, 14, 18)] <- runif(8, 0, 0.001)
    fc <- sample(c(-1, 1), 19, replace = TRUE)
    fc[c(seq_len(3), 13, 14)] <- 1
    fc[c(4, 5, 18)] <- -1
    df <- data.frame(node = seq_len(19),
                     pvalue = pv,
                     foldChange = fc)

    ## Check that function returns error for invalid inputs
    ## -------------------------------------------------------------------------
    .args <- list(tree = tinyTree, t = NULL, score_data = df,
                  node_column = "node", p_column = "pvalue",
                  sign_column = "foldChange", threshold = 0.05, pct_na = 0.5,
                  message = FALSE)

    args <- .args
    args$tree <- 1
    expect_error(do.call(getCand, args),
                 "'tree' must be of class 'phylo'")

    args <- .args
    args$t <- "x"
    expect_error(do.call(getCand, args),
                 "'t' must be of class 'numeric'")

    args <- .args
    args$score_data <- "x"
    expect_error(do.call(getCand, args),
                 "'score_data' must be of class 'data.frame'")

    args <- .args
    args$node_column <- 1
    expect_error(do.call(getCand, args),
                 "'node_column' must be of class 'character'")
    args$node_column <- c("node", "pvalue")
    expect_error(do.call(getCand, args),
                 "'node_column' must have length 1")
    args$node_column <- "missing"
    expect_error(do.call(getCand, args),
                 "all(c(node_column, p_column, sign_column) %in%",
                 fixed = TRUE)

    args <- .args
    args$p_column <- 1
    expect_error(do.call(getCand, args),
                 "'p_column' must be of class 'character'")
    args$p_column <- c("node", "pvalue")
    expect_error(do.call(getCand, args),
                 "'p_column' must have length 1")
    args$p_column <- "missing"
    expect_error(do.call(getCand, args),
                 "all(c(node_column, p_column, sign_column) %in%",
                 fixed = TRUE)

    args <- .args
    args$sign_column <- 1
    expect_error(do.call(getCand, args),
                 "'sign_column' must be of class 'character'")
    args$sign_column <- c("node", "pvalue")
    expect_error(do.call(getCand, args),
                 "'sign_column' must have length 1")
    args$sign_column <- "missing"
    expect_error(do.call(getCand, args),
                 "all(c(node_column, p_column, sign_column) %in%",
                 fixed = TRUE)

    args <- .args
    args$threshold <- "x"
    expect_error(do.call(getCand, args),
                 "'threshold' must be of class 'numeric'")
    args$threshold <- c(0.1, 0.2)
    expect_error(do.call(getCand, args),
                 "'threshold' must have length 1")
    args$threshold <- 2
    expect_error(do.call(getCand, args),
                 "'threshold' must be within [0,1]", fixed = TRUE)

    args <- .args
    args$pct_na <- "x"
    expect_error(do.call(getCand, args),
                 "'pct_na' must be of class 'numeric'")
    args$pct_na <- c(0.1, 0.2)
    expect_error(do.call(getCand, args),
                 "'pct_na' must have length 1")
    args$pct_na <- 2
    expect_error(do.call(getCand, args),
                 "'pct_na' must be within [0,1]", fixed = TRUE)

    args <- .args
    args$message <- "x"
    expect_error(do.call(getCand, args),
                 "'message' must be of class 'logical'")
    args$message <- c(TRUE, FALSE)
    expect_error(do.call(getCand, args),
                 "'message' must have length 1")

    ## Check that function performs as expected for valid inputs
    ## -------------------------------------------------------------------------
    ll <- getCand(tree = tinyTree, score_data = df,
                  t = c(0.01, 0.05, 0.1, 0.25, 0.75),
                  node_column = "node", p_column = "pvalue",
                  sign_column = "foldChange", threshold = 0.05,
                  pct_na = 0.5, message = FALSE)
    expect_type(ll, "list")
    expect_named(ll, c("candidate_list", "score_data"))
    expect_type(ll$candidate_list, "list")
    expect_s3_class(ll$score_data, "data.frame")
    expect_equal(df, ll$score_data[, seq_len(ncol(df))])
    expect_length(ll$candidate_list, 5)
    expect_named(ll$candidate_list,
                 as.character(c(0.01, 0.05, 0.1, 0.25, 0.75)))
    for (i in c(0.01, 0.05, 0.1, 0.25, 0.75)) {
        expect_equal(as.numeric(ll$score_data$pvalue <= i) *
                         ll$score_data$foldChange,
                     ll$score_data[[paste0("q_", i)]])
    }
    expect_equal(ll$candidate_list$`0.05`, c(6, 7, 8, 9, 10, 13, 18))

    ## Candidates are the same even with threshold = 1 - issue is with the sign
    ll <- getCand(tree = tinyTree, score_data = df,
                  t = c(0.01, 0.05, 0.1, 0.25, 0.75),
                  node_column = "node", p_column = "pvalue",
                  sign_column = "foldChange", threshold = 1,
                  pct_na = 0.5, message = FALSE)
    expect_type(ll, "list")
    expect_named(ll, c("candidate_list", "score_data"))
    expect_type(ll$candidate_list, "list")
    expect_s3_class(ll$score_data, "data.frame")
    expect_equal(df, ll$score_data[, seq_len(ncol(df))])
    expect_length(ll$candidate_list, 5)
    expect_named(ll$candidate_list,
                 as.character(c(0.01, 0.05, 0.1, 0.25, 0.75)))
    for (i in c(0.01, 0.05, 0.1, 0.25, 0.75)) {
        expect_equal(as.numeric(ll$score_data$pvalue <= i) *
                         ll$score_data$foldChange,
                     ll$score_data[[paste0("q_", i)]])
    }
    expect_equal(ll$candidate_list$`0.75`, c(6, 7, 8, 9, 10, 13, 18))

    ## Lower threshold so that 13/14 are no longer valid
    ll <- getCand(tree = tinyTree, score_data = df,
                  t = c(0.01, 0.05, 0.1, 0.25, 0.75),
                  node_column = "node", p_column = "pvalue",
                  sign_column = "foldChange", threshold = 1e-4,
                  pct_na = 0.5, message = FALSE)
    expect_type(ll, "list")
    expect_named(ll, c("candidate_list", "score_data"))
    expect_type(ll$candidate_list, "list")
    expect_s3_class(ll$score_data, "data.frame")
    expect_equal(df, ll$score_data[, seq_len(ncol(df))])
    expect_length(ll$candidate_list, 5)
    expect_named(ll$candidate_list,
                 as.character(c(0.01, 0.05, 0.1, 0.25, 0.75)))
    for (i in c(0.01, 0.05, 0.1, 0.25, 0.75)) {
        expect_equal(as.numeric(ll$score_data$pvalue <= i) *
                         ll$score_data$foldChange,
                     ll$score_data[[paste0("q_", i)]])
    }
    expect_equal(ll$candidate_list$`0.05`, c(1, 2, 3, 6, 7, 8, 9, 10, 18))

    ## Change signs so that node 16 is a valid node for a high t threshold
    df2 <- df
    df2$foldChange[df2$node %in% c(17, 7, 19)] <- -1
    ll <- getCand(tree = tinyTree, score_data = df2,
                  t = c(0.01, 0.05, 0.1, 0.25, 0.95),
                  node_column = "node", p_column = "pvalue",
                  sign_column = "foldChange", threshold = 1,
                  pct_na = 0.5, message = FALSE)
    expect_type(ll, "list")
    expect_named(ll, c("candidate_list", "score_data"))
    expect_type(ll$candidate_list, "list")
    expect_s3_class(ll$score_data, "data.frame")
    expect_equal(df2, ll$score_data[, seq_len(ncol(df2))])
    expect_length(ll$candidate_list, 5)
    expect_named(ll$candidate_list,
                 as.character(c(0.01, 0.05, 0.1, 0.25, 0.95)))
    for (i in c(0.01, 0.05, 0.1, 0.25, 0.95)) {
        expect_equal(as.numeric(ll$score_data$pvalue <= i) *
                         ll$score_data$foldChange,
                     ll$score_data[[paste0("q_", i)]])
    }
    expect_equal(ll$candidate_list$`0.95`, c(9, 10, 13, 16))

    ## Add some NAs, set a lenient threshold
    df2 <- df
    df2$pvalue[df2$node %in% c(1, 2, 4)] <- NA
    ll <- getCand(tree = tinyTree, score_data = df2,
                  t = c(0.01, 0.05, 0.1, 0.25, 0.75),
                  node_column = "node", p_column = "pvalue",
                  sign_column = "foldChange", threshold = 0.05,
                  pct_na = 0.25, message = FALSE)
    expect_type(ll, "list")
    expect_named(ll, c("candidate_list", "score_data"))
    expect_type(ll$candidate_list, "list")
    expect_s3_class(ll$score_data, "data.frame")
    expect_equal(df2, ll$score_data[, seq_len(ncol(df2))])
    expect_length(ll$candidate_list, 5)
    expect_named(ll$candidate_list,
                 as.character(c(0.01, 0.05, 0.1, 0.25, 0.75)))
    for (i in c(0.01, 0.05, 0.1, 0.25, 0.75)) {
        expect_equal(as.numeric(ll$score_data$pvalue <= i) *
                         ll$score_data$foldChange,
                     ll$score_data[[paste0("q_", i)]])
    }
    expect_equal(ll$candidate_list$`0.05`, c(6, 7, 8, 9, 10, 13, 18))

    ## With NAs, lenient threshold, stringent threshold on node p-value
    ## -> node 13/14 are no longer allowed. Leaves with p=NA will not be
    ## reported individually -> 1/2 are gone
    df2 <- df
    df2$pvalue[df2$node %in% c(1, 2, 4)] <- NA
    ll <- getCand(tree = tinyTree, score_data = df2,
                  t = c(0.01, 0.05, 0.1, 0.25, 0.75),
                  node_column = "node", p_column = "pvalue",
                  sign_column = "foldChange", threshold = 1e-4,
                  pct_na = 0.25, message = FALSE)
    expect_type(ll, "list")
    expect_named(ll, c("candidate_list", "score_data"))
    expect_type(ll$candidate_list, "list")
    expect_s3_class(ll$score_data, "data.frame")
    expect_equal(df2, ll$score_data[, seq_len(ncol(df2))])
    expect_length(ll$candidate_list, 5)
    expect_named(ll$candidate_list,
                 as.character(c(0.01, 0.05, 0.1, 0.25, 0.75)))
    for (i in c(0.01, 0.05, 0.1, 0.25, 0.75)) {
        expect_equal(as.numeric(ll$score_data$pvalue <= i) *
                         ll$score_data$foldChange,
                     ll$score_data[[paste0("q_", i)]])
    }
    expect_equal(sort(ll$candidate_list$`0.05`),
                 sort(c(3, 6, 7, 8, 9, 10, 18)))

    ## Add some NAs also for internal nodes, set a lenient threshold
    ## p=NA for node 13 -> select node 14 instead
    df2 <- df
    df2$pvalue[df2$node %in% c(1, 2, 4, 13)] <- NA
    ll <- getCand(tree = tinyTree, score_data = df2,
                  t = c(0.01, 0.05, 0.1, 0.25, 0.75),
                  node_column = "node", p_column = "pvalue",
                  sign_column = "foldChange", threshold = 0.05,
                  pct_na = 0.25, message = FALSE)
    expect_type(ll, "list")
    expect_named(ll, c("candidate_list", "score_data"))
    expect_type(ll$candidate_list, "list")
    expect_s3_class(ll$score_data, "data.frame")
    expect_equal(df2, ll$score_data[, seq_len(ncol(df2))])
    expect_length(ll$candidate_list, 5)
    expect_named(ll$candidate_list,
                 as.character(c(0.01, 0.05, 0.1, 0.25, 0.75)))
    for (i in c(0.01, 0.05, 0.1, 0.25, 0.75)) {
        expect_equal(as.numeric(ll$score_data$pvalue <= i) *
                         ll$score_data$foldChange,
                     ll$score_data[[paste0("q_", i)]])
    }
    expect_equal(ll$candidate_list$`0.05`, c(6, 7, 8, 9, 10, 14, 18))

    ## With NAs, more stringent threshold
    ## -> Neither 14 nor 18 is eligible since not at least 75% of the direct
    ## children have valid p-values.
    df2 <- df
    df2$pvalue[df2$node %in% c(1, 2, 4, 13)] <- NA
    ll <- getCand(tree = tinyTree, score_data = df2,
                  t = c(0.01, 0.05, 0.1, 0.25, 0.75),
                  node_column = "node", p_column = "pvalue",
                  sign_column = "foldChange", threshold = 0.05,
                  pct_na = 0.75, message = TRUE)
    expect_type(ll, "list")
    expect_named(ll, c("candidate_list", "score_data"))
    expect_type(ll$candidate_list, "list")
    expect_s3_class(ll$score_data, "data.frame")
    expect_equal(df2, ll$score_data[, seq_len(ncol(df2))])
    expect_length(ll$candidate_list, 5)
    expect_named(ll$candidate_list,
                 as.character(c(0.01, 0.05, 0.1, 0.25, 0.75)))
    for (i in c(0.01, 0.05, 0.1, 0.25, 0.75)) {
        expect_equal(as.numeric(ll$score_data$pvalue <= i) *
                         ll$score_data$foldChange,
                     ll$score_data[[paste0("q_", i)]])
    }
    expect_equal(ll$candidate_list$`0.05`, c(3, 5, 6, 7, 8, 9, 10))
})
