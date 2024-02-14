test_that("getData works", {
    ## Generate some data
    ## -------------------------------------------------------------------------
    library(TreeSummarizedExperiment)
    library(ggtree)
    library(ggplot2)
    library(ggnewscale)
    library(viridis)
    library(dplyr)

    data(tinyTree)

    ## Simulate data with 10 features and 10 samples (5 from each of 2 groups)
    p1 <- c(rep(0.1/3, 3), rep(0.4/2, 2), rep(0.1, 5))
    p2 <- c(rep(0.4/3, 3), rep(0.1/2, 2), rep(0.1, 5))
    set.seed(1)
    ct0 <- cbind(rmultinom(n = 5, size = 50, prob = p1),
                 rmultinom(n = 5, size = 50, prob = p2))
    colnames(ct0) <- paste("S", seq_len(10), sep = "")
    rownames(ct0) <- convertNode(tree = tinyTree, node = seq_len(10))
    oo <- sample(seq_len(10))
    ct0 <- ct0[, oo]
    ct <- rbind(colSums(ct0[seq(1, 3), ]), colSums(ct0[seq(4, 5), ]),
                ct0[seq(6, 10), ])
    colnames(ct) <- paste("S", seq_len(10), sep = "")
    rownames(ct) <- convertNode(tree = tinyTree, node = c(13, 18, seq(6, 10)))
    col_split <- ifelse(colnames(ct) %in% paste0("S", seq_len(5)), "A", "B")
    names(col_split) <- colnames(ct)

    ## Prepare the tree figure
    tree_fig <- ggtree(tinyTree, branch.length = "none",
                       layout = "rectangular", open.angle = 100) +
        geom_hilight(node = 18, fill = "orange", alpha = 0.3) +
        geom_hilight(node = 13, fill = "blue", alpha = 0.3)
    fig <- TreeHeatmap(
        tree = tinyTree, tree_fig = tree_fig, hm_data = ct,
        cluster_column = TRUE, column_split = col_split,
        column_anno = col_split, column_anno_gap = 0.6,
        column_anno_color = c(A = "red", B = "blue"),
        show_colnames = TRUE, colnames_position = "bottom",
        colnames_angle = 90, colnames_size = 2, colnames_offset_y = -0.2,
        show_title = TRUE, title_offset_y = 1.5, title_color = "blue"
    )

    ## Figure with row names
    figrn <- TreeHeatmap(
        tree = tinyTree, tree_fig = tree_fig, hm_data = ct,
        cluster_column = TRUE, column_split = col_split,
        column_anno = col_split, column_anno_gap = 0.6,
        column_anno_color = c(A = "red", B = "blue"),
        show_colnames = TRUE, colnames_position = "bottom",
        colnames_angle = 90, colnames_size = 2, colnames_offset_y = -0.2,
        show_title = TRUE, title_offset_y = 1.5, title_color = "blue",
        show_rownames = TRUE
    )

    ## Check that the function fails with wrong input
    ## -------------------------------------------------------------------------
    expect_error(getData(tree_hm = 1, type = "heatmap"),
                 "'tree_hm' must be of class 'ggtree'")
    expect_error(getData(tree_hm = fig, type = "missing"),
                 "'arg' should be one of")

    ## Check that the function works with correct input
    ## -------------------------------------------------------------------------
    ## Heatmap data
    df_hm <- getData(tree_hm = fig, type = "heatmap")
    expect_s3_class(df_hm, "data.frame")
    expect_equal(nrow(df_hm), 70) ## 7 blocks, 10 samples
    expect_named(df_hm, c("node", "row_label", "x", "y", "height", "width",
                          "variable", "value", "column_order", "split_level"))
    expect_equal(df_hm$row_label, rep(rownames(ct), each = 10))
    expect_equal(df_hm$column_order, rep(c(2, 3, 4, 0, 1, 5, 6, 8, 9, 7), 7))
    expect_equal(df_hm$variable, factor(rep(paste0("S", seq_len(10)), 7),
                                        levels = c("S4", "S5", "S1", "S2", "S3",
                                                   "S6", "S7", "S10", "S8",
                                                   "S9")))
    expect_true(all(diff(df_hm$x[order(df_hm$column_order)]) >= 0))
    expect_length(unique(paste(df_hm$column_order, df_hm$x)), 10)
    expect_equal(unique(df_hm$width), 0.6)
    expect_equal(unique(df_hm$split_level[df_hm$variable == "S1"]), 0)
    expect_equal(unique(df_hm$split_level[df_hm$variable == "S8"]), 1)
    expect_equal(unique(df_hm$height[df_hm$row_label == "Node_13"]), 3)
    expect_equal(unique(df_hm$height[df_hm$row_label == "Node_18"]), 2)
    expect_equal(unique(df_hm$height[df_hm$row_label == "t8"]), 1)
    for (n in rownames(ct)) {
        for (m in colnames(ct)) {
            expect_equal(df_hm$value[df_hm$variable == m &
                                         df_hm$row_label == n],
                         ct[n, m])
        }
    }

    ## ... from figure with row names
    df_hm_rn <- getData(tree_hm = figrn, type = "heatmap")
    expect_equal(df_hm, df_hm_rn)

    ## Row names
    df_rn <- getData(tree_hm = fig, type = "row_name")
    expect_null(df_rn)
    df_rn2 <- getData(tree_hm = figrn, type = "row_name")
    expect_s3_class(df_rn2, "data.frame")
    expect_equal(nrow(df_rn2), 7)
    expect_equal(df_rn2$row_label[order(df_rn2$y)], c("t3", "Node_13", "t5",
                                                      "t10", "t1", "t8", "Node_18"))
    expect_equal(unique(df_rn2$width), 0.6)
    expect_equal(unique(df_rn2$x_right), 12.2)
    expect_equal(unique(df_rn2$x_left), 6)
    expect_equal(unique(df_rn2$x), 12.2)

    ## Column names
    df_cn <- getData(tree_hm = fig, type = "column_name")
    expect_s3_class(df_cn, "data.frame")
    expect_equal(nrow(df_cn), 10)
    expect_named(df_cn, c("variable", "x", "width", "y_top", "y_bottom", "y"))
    expect_equal(diff(sort(df_cn$x)), c(rep(0.6, 4), 0.8, rep(0.6, 4)))
    expect_equal(df_cn$variable, paste0("S", seq_len(10)))
    expect_equal(unique(df_cn$y), 0.5)
    expect_equal(unique(df_cn$y_top), 11.1)
    expect_equal(unique(df_cn$y_bottom), 0.5)
    expect_equal(unique(df_cn$width), 0.6)

    ## Title
    df_title <- getData(tree_hm = fig, type = "title")
    expect_s3_class(df_title, "data.frame")
    expect_equal(df_title, data.frame(x = 9.1, y = 9.5, label = "First heatmap"))

    ## Column annotation
    df_ca <- getData(tree_hm = fig, type = "column_anno")
    expect_s3_class(df_ca, "data.frame")
    expect_equal(nrow(df_ca), 10)
    expect_named(df_ca, c("variable", "x", "width", "xend", "y", "yend",
                          "anno_group", "anno_color"))
    expect_equal(df_ca$variable, paste0("S", seq_len(10)))
    expect_equal(df_ca$x - df_cn$x, rep(-0.3, 10))
    expect_equal(df_ca$x - df_ca$xend, rep(-0.6, 10))
    expect_equal(df_ca$anno_group, rep(c("A", "B"), each = 5),
                 ignore_attr = TRUE)
    expect_equal(df_ca$anno_color, rep(c("red", "blue"), each = 5),
                 ignore_attr = TRUE)

    ## Column order
    df_cord <- getData(tree_hm = fig, type = "column_order")
    expect_type(df_cord, "character")
    ## Samples are clustered
    expect_equal(df_cord, c("S4", "S5", "S1", "S2", "S3",
                            "S6", "S7", "S10", "S8", "S9"))

    ## Column split
    df_spl <- getData(tree_hm = fig, type = "column_split")
    expect_s3_class(df_spl, "tbl_df")
    expect_equal(nrow(df_spl), 2)
    expect_named(df_spl, c("split_level", "x", "y", "column_split"))
    expect_equal(df_spl$column_split, c("A", "B"))

})
