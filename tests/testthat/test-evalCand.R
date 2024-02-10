test_that(".pseudoLeaf works", {
    library(TreeSummarizedExperiment)
    library(ggtree)
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

    expect_equal(.pseudoLeaf(tree = tinyTree, score_data = df,
                             node_column = "node", p_column = "pvalue"),
                 seq(1, 10))

    ## 2,3 are NA - replace with next level
    df2 <- df
    df2$pvalue[df2$node %in% c(2, 3)] <- NA
    expect_equal(.pseudoLeaf(tree = tinyTree, score_data = df2,
                             node_column = "node", p_column = "pvalue"),
                 c(1, 14, seq(4, 10)))

    ## only 2 is NA - gets removed
    df2 <- df
    df2$pvalue[df2$node %in% 2] <- NA
    expect_equal(.pseudoLeaf(tree = tinyTree, score_data = df2,
                             node_column = "node", p_column = "pvalue"),
                 c(1, 3, seq(4, 10)))

    ## 4,5,6 are NA - replace 4 and 5 with 18 and remove 6
    df2 <- df
    df2$pvalue[df2$node %in% c(4, 5, 6)] <- NA
    expect_equal(.pseudoLeaf(tree = tinyTree, score_data = df2,
                             node_column = "node", p_column = "pvalue"),
                 c(seq(1, 3), 18, seq(7, 10)))

})

test_that("evalCand works", {
    ## Generate some data
    library(TreeSummarizedExperiment)
    data(tinyTree)
    set.seed(1)
    pv <- runif(19, 0, 1)
    pv[c(seq_len(5), 13, 14, 18)] <- runif(8, 0, 0.001)
    fc <- sample(c(-1, 1), 19, replace = TRUE)
    fc[c(seq_len(3), 13, 14)] <- 1
    fc[c(4, 5, 18)] <- -1
    df <- data.frame(node = seq_len(19), pvalue = pv,
                     foldChange = fc, feature = "gene1")
    ll <- getCand(tree = tinyTree, score_data = df,
                  node_column = "node", p_column = "pvalue",
                  sign_column = "foldChange")

    ## One with NAs
    df2 <- df
    df2$feature <- "gene2"
    df2$pvalue[df2$node %in% c(4, 5, 6)] <- NA
    llna <- getCand(tree = tinyTree, score_data = df2,
                    node_column = "node", p_column = "pvalue",
                    sign_column = "foldChange", threshold = 1e-4,
                    pct_na = 0.25, message = FALSE)

    ## Check that the function returns error for invalid input
    ## -------------------------------------------------------------------------
    .args <- list(tree = tinyTree, type = "single",
                  levels = ll$candidate_list, score_data = ll$score_data,
                  node_column = "node", p_column = "pvalue",
                  sign_column = "foldChange",
                  feature_column = NULL, method = "BH",
                  limit_rej = 0.05, use_pseudo_leaf = FALSE,
                  message = FALSE)

    args <- .args
    args$tree <- 1
    expect_error(do.call(evalCand, args),
                 "'tree' must be of class 'phylo'")

    args <- .args
    args$type <- 1
    expect_error(do.call(evalCand, args),
                 "'arg' must be NULL or a character vector")
    args$type <- "missing"
    expect_error(do.call(evalCand, args),
                 "'arg' should be one of")

    args <- .args
    args$levels <- c(1, 2)
    expect_error(do.call(evalCand, args),
                 "'levels' must be of class 'list'")

    args <- .args
    args$score_data <- c(1, 2)
    expect_error(do.call(evalCand, args),
                 "'score_data' must be of class 'data.frame'")
    args$type <- "multiple"
    expect_error(do.call(evalCand, args),
                 "'score_data' must be of class 'list'")

    args <- .args
    args$node_column <- 1
    expect_error(do.call(evalCand, args),
                 "'node_column' must be of class 'character'")
    args$node_column <- c("node", "pvalue")
    expect_error(do.call(evalCand, args),
                 "'node_column' must have length 1")
    args$node_column <- "missing"
    expect_error(do.call(evalCand, args),
                 "all(c(node_column, p_column, sign_column) %in%",
                 fixed = TRUE)

    args <- .args
    args$p_column <- 1
    expect_error(do.call(evalCand, args),
                 "'p_column' must be of class 'character'")
    args$p_column <- c("node", "pvalue")
    expect_error(do.call(evalCand, args),
                 "'p_column' must have length 1")
    args$p_column <- "missing"
    expect_error(do.call(evalCand, args),
                 "all(c(node_column, p_column, sign_column) %in%",
                 fixed = TRUE)

    args <- .args
    args$sign_column <- 1
    expect_error(do.call(evalCand, args),
                 "'sign_column' must be of class 'character'")
    args$sign_column <- c("node", "pvalue")
    expect_error(do.call(evalCand, args),
                 "'sign_column' must have length 1")
    args$sign_column <- "missing"
    expect_error(do.call(evalCand, args),
                 "all(c(node_column, p_column, sign_column) %in%",
                 fixed = TRUE)

    args <- .args
    args$type <- "multiple"
    args$score_data <- list(args$score_data)
    args$feature_column <- 1
    expect_error(do.call(evalCand, args),
                 "'feature_column' must be of class 'character'")
    args$feature_column <- c("node", "pvalue")
    expect_error(do.call(evalCand, args),
                 "'feature_column' must have length 1")
    args$feature_column <- "missing"
    expect_error(do.call(evalCand, args),
                 "feature_column %in%",
                 fixed = TRUE)

    args <- .args
    args$method <- 1
    expect_error(do.call(evalCand, args),
                 "'method' must be of class 'character'")
    args$method <- c("BH", "holm")
    expect_error(do.call(evalCand, args),
                 "'method' must have length 1")

    args <- .args
    args$limit_rej <- "x"
    expect_error(do.call(evalCand, args),
                 "'limit_rej' must be of class 'numeric'")
    args$limit_rej <- c(0.1, 0.2)
    expect_error(do.call(evalCand, args),
                 "'limit_rej' must have length 1")
    args$limit_rej <- 2
    expect_error(do.call(evalCand, args),
                 "'limit_rej' must be within [0,1]", fixed = TRUE)

    args <- .args
    args$use_pseudo_leaf <- "x"
    expect_error(do.call(evalCand, args),
                 "'use_pseudo_leaf' must be of class 'logical'")
    args$use_pseudo_leaf <- c(TRUE, FALSE)
    expect_error(do.call(evalCand, args),
                 "'use_pseudo_leaf' must have length 1")

    args <- .args
    args$message <- "x"
    expect_error(do.call(evalCand, args),
                 "'message' must be of class 'logical'")
    args$message <- c(TRUE, FALSE)
    expect_error(do.call(evalCand, args),
                 "'message' must have length 1")

    ## Check that the function works as expected for valid input
    ## -------------------------------------------------------------------------
    ## Single
    out <- evalCand(tree = tinyTree, type = "single",
                    levels = ll$candidate_list, score_data = ll$score_data,
                    node_column = "node", p_column = "pvalue",
                    sign_column = "foldChange",
                    feature_column = NULL, method = "BH",
                    limit_rej = 0.05, use_pseudo_leaf = FALSE,
                    message = FALSE)
    expect_type(out, "list")
    expect_named(out, c("candidate_best", "output", "candidate_list",
                        "level_info", "FDR", "method", "column_info"))
    expect_type(out$candidate_best, "list")
    expect_s3_class(out$output, "data.frame")
    expect_type(out$candidate_list, "list")
    expect_s3_class(out$level_info, "data.frame")
    expect_equal(out$FDR, 0.05)
    expect_equal(out$method, "BH")
    expect_equal(out$column_info,
                 list(node_column = "node", p_column = "pvalue",
                      sign_column = "foldChange", feature_column = NULL))
    expect_equal(out$level_info$rej_pseudo_leaf, rep(NA, nrow(out$level_info)))
    expect_equal(out$level_info[2, ],
                 data.frame(t = 0.01, upper_t = 0.15, is_valid = TRUE,
                            method = "BH", limit_rej = 0.05, level_name = "0.01",
                            best = TRUE, rej_leaf = 5, rej_node = 2,
                            rej_pseudo_leaf = NA, rej_pseudo_node = NA),
                 ignore_attr = TRUE)
    expect_equal(nrow(out$output), length(ll$candidate_list$`0.01`))
    expect_equal(as.data.frame(out$output)[, seq_len(7)],
                 ll$score_data[match(out$output$node, ll$score_data$node), seq_len(7)],
                 ignore_attr = TRUE)

    ## Multiple
    out <- evalCand(tree = tinyTree, type = "multiple",
                    levels = list(gene1 = ll$candidate_list,
                                  gene2 = llna$candidate_list),
                    score_data = list(gene1 = ll$score_data,
                                      gene2 = llna$score_data),
                    node_column = "node", p_column = "pvalue",
                    sign_column = "foldChange",
                    feature_column = "feature", method = "BH",
                    limit_rej = 0.05, use_pseudo_leaf = FALSE,
                    message = FALSE)
    expect_type(out, "list")
    expect_named(out, c("candidate_best", "output", "candidate_list",
                        "level_info", "FDR", "method", "column_info"))
    expect_type(out$candidate_best, "list")
    expect_s3_class(out$output, "data.frame")
    expect_type(out$candidate_list, "list")
    expect_s3_class(out$level_info, "data.frame")
    expect_equal(out$FDR, 0.05)
    expect_equal(out$method, "BH")
    expect_equal(out$column_info,
                 list(node_column = "node", p_column = "pvalue",
                      sign_column = "foldChange", feature_column = "feature"))
    expect_equal(out$level_info$rej_pseudo_leaf, rep(NA, nrow(out$level_info)))
    expect_equal(out$level_info[2, ],
                 data.frame(t = 0.01, upper_t = 0.1, is_valid = TRUE,
                            method = "BH", limit_rej = 0.05, level_name = "0.01",
                            best = TRUE, rej_leaf = 10, rej_node = 6,
                            rej_pseudo_leaf = NA, rej_pseudo_node = NA),
                 ignore_attr = TRUE)
    expect_equal(nrow(out$output), length(ll$candidate_list$`0.01`) +
                     length(llna$candidate_list$`0.01`))
    expect_equal(out$output[out$output$feature == "gene1", seq_len(7)],
                 ll$score_data[match(out$output$node[out$output$feature == "gene1"],
                                     ll$score_data$node), seq_len(7)],
                 ignore_attr = TRUE)

    ## Multiple - use pseudo leaf
    out <- evalCand(tree = tinyTree, type = "multiple",
                    levels = list(gene1 = ll$candidate_list,
                                  gene2 = llna$candidate_list),
                    score_data = list(gene1 = ll$score_data,
                                      gene2 = llna$score_data),
                    node_column = "node", p_column = "pvalue",
                    sign_column = "foldChange",
                    feature_column = "feature", method = "BH",
                    limit_rej = 0.05, use_pseudo_leaf = TRUE,
                    message = TRUE)
    expect_type(out, "list")
    expect_named(out, c("candidate_best", "output", "candidate_list",
                        "level_info", "FDR", "method", "column_info"))
    expect_type(out$candidate_best, "list")
    expect_s3_class(out$output, "data.frame")
    expect_type(out$candidate_list, "list")
    expect_s3_class(out$level_info, "data.frame")
    expect_equal(out$FDR, 0.05)
    expect_equal(out$method, "BH")
    expect_equal(out$column_info,
                 list(node_column = "node", p_column = "pvalue",
                      sign_column = "foldChange", feature_column = "feature"))
    expect_equal(out$level_info$rej_pseudo_leaf, rep(9, nrow(out$level_info)))
    expect_equal(out$level_info[2, ],
                 data.frame(t = 0.01, upper_t = 0.08, is_valid = TRUE,
                            method = "BH", limit_rej = 0.05, level_name = "0.01",
                            best = TRUE, rej_leaf = 9, rej_node = 6,
                            rej_pseudo_leaf = 9, rej_pseudo_node = 5),
                 ignore_attr = TRUE)
    expect_equal(nrow(out$output), length(ll$candidate_list$`0.01`) +
                     length(llna$candidate_list$`0.01`))
    expect_equal(out$output[out$output$feature == "gene1", seq_len(7)],
                 ll$score_data[match(out$output$node[out$output$feature == "gene1"],
                                     ll$score_data$node), seq_len(7)],
                 ignore_attr = TRUE)

    ## Multiple - different levels for different features
    ll1 <- ll
    ll1$candidate_list <- ll1$candidate_list[seq_len(3)]
    expect_error(evalCand(tree = tinyTree, type = "multiple",
                          levels = list(gene1 = ll1$candidate_list,
                                        gene2 = llna$candidate_list),
                          score_data = list(gene1 = ll$score_data,
                                            gene2 = llna$score_data),
                          node_column = "node", p_column = "pvalue",
                          sign_column = "foldChange",
                          feature_column = "feature", method = "BH",
                          limit_rej = 0.05, use_pseudo_leaf = FALSE,
                          message = FALSE),
                 "The names of elements in 'levels' are different")

    ## Multiple - no feature column
    expect_warning(evalCand(tree = tinyTree, type = "multiple",
                            levels = list(gene1 = ll$candidate_list,
                                          gene2 = llna$candidate_list),
                            score_data = list(gene1 = ll$score_data,
                                              gene2 = llna$score_data),
                            node_column = "node", p_column = "pvalue",
                            sign_column = "foldChange",
                            feature_column = NULL, method = "BH",
                            limit_rej = 0.05, use_pseudo_leaf = FALSE,
                            message = FALSE),
                 "To distinguish results from different features")
})
