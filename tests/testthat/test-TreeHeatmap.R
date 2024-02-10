test_that("TreeHeatmap works", {
    ## Get some data
    ## -------------------------------------------------------------------------
    library(TreeSummarizedExperiment)
    library(ggtree)
    library(ggplot2)
    library(scales)

    ## Load example data (tiny tree with corresponding count matrix)
    tse <- readRDS(system.file("extdata", "tinytree_counts.rds",
                               package = "treeclimbR"))

    ## Prepare the tree figure
    tree_fig <- ggtree(rowTree(tse), branch.length = "none",
                       layout = "rectangular", open.angle = 100) +
        geom_hilight(node = 18, fill = "orange", alpha = 0.3) +
        geom_hilight(node = 13, fill = "blue", alpha = 0.3)

    ## Aggregate counts for each of the highlighted subtrees
    tseagg <- aggTSE(
        tse,
        rowLevel = c(13, 18,
                     setdiff(showNode(tinyTree, only.leaf = TRUE),
                             unlist(findDescendant(tinyTree, node = c(13, 18),
                                                   only.leaf = TRUE)))))

    ## Check that the function fails if provided with invalid input data
    ## -------------------------------------------------------------------------
    tmp <- SummarizedExperiment::assay(tseagg, "counts")
    rownames(tmp)[rownames(tmp) == "alias_8"] <- "alias_4"
    expect_warning(
        TreeHeatmap(tree = rowTree(tse), tree_fig = tree_fig,
                    hm_data = tmp),
        "Some leaves are contributing to multiple rows")

    col_order <- colnames(tse)
    col_order[1] <- "missing"
    expect_error(
        TreeHeatmap(tree = rowTree(tse), tree_fig = tree_fig,
                    hm_data = SummarizedExperiment::assay(tse, "counts"),
                    column_order = col_order),
        "column_order: Some columns are not present in hm_data")

    col_split <- rep(c("A", "B"), each = 5)
    expect_error(
        TreeHeatmap(tree = rowTree(tse), tree_fig = tree_fig,
                    hm_data = SummarizedExperiment::assay(tse, "counts"),
                    column_split = col_split),
        "column_split: should be a named vector")
    names(col_split) <- colnames(tse)
    names(col_split)[1] <- "missing"
    expect_error(
        TreeHeatmap(tree = rowTree(tse), tree_fig = tree_fig,
                    hm_data = SummarizedExperiment::assay(tse, "counts"),
                    column_split = col_split),
        "column_split: Some columns are not present in hm_data")

    col_anno_color <- c("blue", "red")
    expect_error(
        TreeHeatmap(tree = rowTree(tse), tree_fig = tree_fig,
                    hm_data = SummarizedExperiment::assay(tse, "counts"),
                    column_anno_color = col_anno_color),
        "column_anno_color: should be a named vector")

    row_label <- rep(c("C", "D"), each = 5)
    expect_error(
        TreeHeatmap(tree = rowTree(tse), tree_fig = tree_fig,
                    hm_data = SummarizedExperiment::assay(tse, "counts"),
                    rownames_label = row_label),
        "rownames_label: should be a named vector")
    names(row_label) <- rownames(tse)
    names(row_label)[1] <- "missing"
    expect_error(
        TreeHeatmap(tree = rowTree(tse), tree_fig = tree_fig,
                    hm_data = SummarizedExperiment::assay(tse, "counts"),
                    rownames_label = row_label),
        "rownames_label: Some rows are not present in hm_data")

    ## Check that the function works as expected
    ## -------------------------------------------------------------------------
    ## Simple heatmap with tree
    hm <- TreeHeatmap(tree = rowTree(tse), tree_fig = tree_fig,
                      hm_data = SummarizedExperiment::assay(tse, "counts"))
    expect_s3_class(hm, "ggtree")

    ## Don't show row tree
    hm <- TreeHeatmap(tree = rowTree(tse), tree_fig = tree_fig,
                      hm_data = SummarizedExperiment::assay(tse, "counts"),
                      show_row_tree = FALSE)
    expect_s3_class(hm, "ggplot")

    ## Character matrix
    mat <- SummarizedExperiment::assay(tse, "counts")
    mode(mat) <- "character"
    hm <- TreeHeatmap(tree = rowTree(tse), tree_fig = tree_fig,
                      hm_data = mat)
    expect_s3_class(hm, "ggtree")

    ## Aggregated heatmap with tree
    hm <- TreeHeatmap(tree = rowTree(tseagg), tree_fig = tree_fig,
                      hm_data = SummarizedExperiment::assay(tseagg, "counts"))
    expect_s3_class(hm, "ggtree")

    ## Cluster columns
    hm <- TreeHeatmap(tree = rowTree(tseagg), tree_fig = tree_fig,
                      hm_data = SummarizedExperiment::assay(tseagg, "counts"),
                      cluster_column = TRUE)
    expect_s3_class(hm, "ggtree")

    ## Cluster columns for subset, with split
    col_split <- ifelse(colnames(tseagg) %in% paste0("S", seq_len(5)), "A", "B")
    names(col_split) <- colnames(tseagg)
    idx <- c(1, 3, 4)
    hm <- TreeHeatmap(
        tree = rowTree(tseagg), tree_fig = tree_fig,
        hm_data = SummarizedExperiment::assay(tseagg, "counts")[, idx],
        cluster_column = TRUE, column_split = col_split[idx])
    expect_s3_class(hm, "ggtree")

    ## Split columns
    col_split <- ifelse(colnames(tseagg) %in% paste0("S", seq_len(5)), "A", "B")
    names(col_split) <- colnames(tseagg)
    hm <- TreeHeatmap(tree = rowTree(tseagg), tree_fig = tree_fig,
                      hm_data = SummarizedExperiment::assay(tseagg, "counts"),
                      cluster_column = TRUE, column_split = col_split)
    expect_s3_class(hm, "ggtree")

    ## Split columns, new labels
    col_split <- ifelse(colnames(tseagg) %in% paste0("S", seq_len(5)), "A", "B")
    names(col_split) <- colnames(tseagg)
    col_split_label <- c(A = "hello", B = "goodbye")
    hm <- TreeHeatmap(tree = rowTree(tseagg), tree_fig = tree_fig,
                      hm_data = SummarizedExperiment::assay(tseagg, "counts"),
                      cluster_column = TRUE, column_split = col_split,
                      column_split_label = col_split_label)
    expect_s3_class(hm, "ggtree")

    ## Provide both column split and column order - ignore the latter
    col_split <- ifelse(colnames(tseagg) %in% paste0("S", seq_len(5)), "A", "B")
    names(col_split) <- colnames(tseagg)
    col_order <- sample(colnames(tseagg), 10)
    expect_warning({
        hm <- TreeHeatmap(tree = rowTree(tseagg), tree_fig = tree_fig,
                          hm_data = SummarizedExperiment::assay(tseagg, "counts"),
                          cluster_column = FALSE, column_split = col_split,
                          column_order = col_order)
    }, "column_order is ignored when column_split is given")
    expect_s3_class(hm, "ggtree")

    ## Annotate columns
    col_anno <- col_split
    hm <- TreeHeatmap(tree = rowTree(tseagg), tree_fig = tree_fig,
                      hm_data = SummarizedExperiment::assay(tseagg, "counts"),
                      cluster_column = TRUE, column_split = col_split,
                      column_anno = col_anno, column_anno_gap = 1)
    expect_s3_class(hm, "ggtree")

    ## Change annotation colors
    hm <- TreeHeatmap(tree = rowTree(tseagg), tree_fig = tree_fig,
                      hm_data = SummarizedExperiment::assay(tseagg, "counts"),
                      cluster_column = TRUE, column_split = col_split,
                      column_anno = col_anno, column_anno_gap = 1,
                      column_anno_color = c(A = "red", B = "blue"))
    expect_s3_class(hm, "ggtree")

    ## Add column names
    hm <- TreeHeatmap(tree = rowTree(tseagg), tree_fig = tree_fig,
                      hm_data = SummarizedExperiment::assay(tseagg, "counts"),
                      cluster_column = TRUE, column_split = col_split,
                      column_anno = col_anno, column_anno_gap = 1,
                      column_anno_color = c(A = "red", B = "blue"),
                      show_colnames = TRUE, colnames_position = "bottom",
                      colnames_angle = 90, colnames_size = 2,
                      colnames_offset_y = -0.4)
    expect_s3_class(hm, "ggtree")

    ## Add row labels
    rl <- letters[seq_len(nrow(tseagg))]
    names(rl) <- rownames(tseagg)
    hm <- TreeHeatmap(tree = rowTree(tseagg), tree_fig = tree_fig,
                      hm_data = SummarizedExperiment::assay(tseagg, "counts"),
                      cluster_column = TRUE, column_split = col_split,
                      column_anno = col_anno, column_anno_gap = 1,
                      column_anno_color = c(A = "red", B = "blue"),
                      show_colnames = TRUE, colnames_position = "bottom",
                      colnames_angle = 90, colnames_size = 2,
                      colnames_offset_y = -0.4, rownames_label = rl,
                      show_rownames = TRUE)
    expect_s3_class(hm, "ggtree")

    ## Add title
    hm <- TreeHeatmap(tree = rowTree(tseagg), tree_fig = tree_fig,
                      hm_data = SummarizedExperiment::assay(tseagg, "counts"),
                      cluster_column = TRUE, column_split = col_split,
                      column_anno = col_anno, column_anno_gap = 1,
                      column_anno_color = c(A = "red", B = "blue"),
                      show_colnames = TRUE, colnames_position = "bottom",
                      colnames_angle = 90, colnames_size = 2,
                      colnames_offset_y = -0.4,
                      show_title = TRUE, title_offset_y = 2,
                      title_color = "blue")
    expect_s3_class(hm, "ggtree")

    ## Change colors
    hm <- TreeHeatmap(tree = rowTree(tseagg), tree_fig = tree_fig,
                      hm_data = SummarizedExperiment::assay(tseagg, "counts"),
                      cluster_column = TRUE, column_split = col_split,
                      column_anno = col_anno, column_anno_gap = 1,
                      column_anno_color = c(A = "red", B = "blue"),
                      show_colnames = TRUE, colnames_position = "bottom",
                      colnames_angle = 90, colnames_size = 2,
                      colnames_offset_y = -0.4,
                      show_title = TRUE, title_offset_y = 2,
                      title_color = "blue") +
        scale_fill_gradientn(
            colours = c("blue", "yellow", "red"),
            values = scales::rescale(c(5, 8, 10)),
            guide = "colorbar", limits = c(5, 10))
    expect_s3_class(hm, "ggtree")
})
