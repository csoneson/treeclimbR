#' Generate a heatmap corresponding to an arbitrary level of a tree
#'
#' Generate a heatmap corresponding to an arbitrary level of a tree.
#'
#' @author Ruizhu Huang
#' @export
#'
#' @param tree A \code{phylo} object.
#' @param tree_fig A \code{ggtree} object corresponding to \code{tree}. This
#'     will be used to represent the tree in the resulting figure.
#' @param hm_data A \code{data.frame} with the values to show in the heatmap.
#'     The row names should correspond to the nodes of \strong{tree}.
#' @param tree_hm_gap A numeric scalar specifying the gap between the tree and
#'     the heatmap.
#' @param rel_width A numeric scalar specifying the width of heatmap relative to
#'     the width of the tree. For example, if \code{rel_width = 1}, the width of
#'     the heatmap is the same as the width of the tree.
#' @param cell_line_color A color for the lines separating cells in the
#'     heatmap. The default is NA.
#' @param cell_line_size A numeric scalar specifying the line width for lines
#'     separating cells in the heatmap. The default is 0.
#' @param column_order A character vector specifying the display order of the
#'     columns in the heatmap. Should correspond to the column names of
#'     \code{hm_data}. Ignored when \strong{column_split} is provided.
#' @param column_split A named character vector that provides the grouping
#'     information used to split the columns in the heatmap. The names should
#'     correspond to the column names of \code{hm_data}.
#' @param column_split_label A named character vector to label the column split.
#'     The names should correspond to the values in \code{column_split}.
#' @param column_split_gap A numeric scalar specifying the gap between the
#'     groups of split columns in the heatmap.
#' @param split_label_fontface The fontface of the labels of the column split.
#'     The default is "bold".
#' @param split_label_color The color of the the labels of the column split.
#'     The default is "black".
#' @param split_label_size The size of the the labels of the column split. The
#'     default is 3.
#' @param split_label_angle The angle of the the labels of the column split. The
#'     default is 0.
#' @param split_label_offset_x A numeric value to shift the labels of the column
#'     split along the x-axis. The default is 0.
#' @param split_label_offset_y A numeric value to shift the labels of the column
#'     split along the y-axis. The default is 2.
#' @param split_label_hjust The horizontal justification for the labels of the
#'     column split: 0 (left aligned); 0.5 (centered); 1 (right aligned).
#'     .The default is 0.5
#' @param split_label_vjust Similar to \code{split_label_hjust}, but controls
#'     vertical justification.
#' @param column_anno A named vector to specify labels that are used to
#'     annotate the columns of heatmap.
#' @param column_anno_size A numeric value to specify the size of the annotation
#'     bar.
#' @param column_anno_color A named vector to specify colors that are used to
#'     annotate the columns of the heatmap.
#' @param column_anno_gap A numeric value to specify the gap between the
#'     column annotation bar and the heatmap.
#' @param legend_title_hm The legend title of the heatmap.
#' @param legend_title_column_anno The legend title of the column annotation.
#' @param show_colnames A logical value to specify whether column names should
#'     be displayed. The default is FALSE.
#' @param colnames_position The position of column names, either "top" or
#'     "bottom".
#' @param colnames_angle A numeric scalar specifying the angle of column names.
#' @param colnames_offset_x A numeric value to shift column names on the x-axis.
#'     The default is 0.
#' @param colnames_offset_y A numeric value to shift column names on the y-axis.
#'     The default is 0.
#' @param colnames_size A numeric value to specify the size of column names.
#' @param colnames_hjust The horizontal justification for column names:
#'     0 (left aligned); 0.5 (centered); 1 (right aligned).
#' @param show_rownames A logical value to specify whether row names should
#'     be displayed. The default is FALSE.
#' @param rownames_position The position of the row names, either "right" or
#'     "left".
#' @param rownames_label A named vector to annotate the rows of the heatmap
#'     instead of the row names of \strong{hm_data}.
#' @param rownames_angle A numeric value specifying the angle of row names.
#' @param rownames_offset_x A numeric value to shift row names on the x-axis.
#'     The default is 0.
#' @param rownames_offset_y A numeric value to shift row names on the y-axis.
#'     The default is 0.
#' @param rownames_size A numeric value to specify the size of row names.
#' @param rownames_hjust The horizontal justification for row names:
#'     0 (left aligned); 0.5 (centered); 1 (right aligned).
#' @param show_title A logical value to specify whether the title should
#'     be displayed. The default is FALSE.
#' @param title_hm The title of the heatmap.
#' @param title_fontface The fontface of the title. The default is "bold".
#' @param title_color  The color of the title. The default is "black".
#' @param title_size The size of the title. The default is 3.
#' @param title_angle The angle of the title. The default is 0.
#' @param title_offset_x A numeric value to shift the title along the x-axis.
#'     The default is 0.
#' @param title_offset_y A numeric value to shift the title along the y-axis.
#'     The default is 2.
#' @param title_hjust The horizontal justification for the title:
#'     0 (left aligned); 0.5 (centered); 1 (right aligned). The default is 0.5.
#' @param cluster_column A logical scalar, specifying whether
#'     columns of the heatmap should be clustered by similarity. The default is
#'     \code{TRUE}. This is ignored when \strong{column_order} is given.
#' @param dist_method See \strong{method} in \code{\link[stats]{dist}}. The
#'     distance method used for clustering columns. The default is "euclidean".
#' @param hclust_method See \strong{method} in \code{\link[stats]{hclust}}. The
#'     clustering method used for clustering columns. The default is "ave".
#' @param show_row_tree A logical scalar (default \code{TRUE}). If \code{FALSE},
#'     the figure provided in \code{tree_fig} is not shown.
#'
#' @returns A \code{ggtree} object.
#'
#' @importFrom TreeSummarizedExperiment convertNode findDescendant
#' @importFrom ggtree ggtree
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr mutate select distinct group_by summarise arrange
#' @importFrom ggplot2 geom_tile geom_segment scale_color_manual labs
#'     geom_text scale_fill_viridis_c aes scale_fill_viridis_d theme_void ggplot
#' @importFrom ggnewscale new_scale_color
#' @importFrom viridis viridis
#' @importFrom stats hclust dist
#'
#' @examples
#' library(TreeSummarizedExperiment)
#' library(ggtree)
#' library(ggplot2)
#' library(scales)
#'
#' ## Load example data (tiny tree with corresponding count matrix)
#' tse <- readRDS(system.file("extdata", "tinytree_counts.rds",
#'                            package = "treeclimbR"))
#'
#' ## Prepare the tree figure
#' tree_fig <- ggtree(rowTree(tse), branch.length = "none",
#'                    layout = "rectangular", open.angle = 100) +
#'     geom_hilight(node = 18, fill = "orange", alpha = 0.3) +
#'     geom_hilight(node = 13, fill = "blue", alpha = 0.3)
#' tree_fig
#'
#' ## Simple heatmap with tree
#' TreeHeatmap(tree = rowTree(tse), tree_fig = tree_fig,
#'             hm_data = SummarizedExperiment::assay(tse, "counts"))
#'
#' ## Aggregate counts for each of the highlighted subtrees
#' tseagg <- aggTSE(
#'     tse,
#'     rowLevel = c(13, 18,
#'                  setdiff(showNode(tinyTree, only.leaf = TRUE),
#'                          unlist(findDescendant(tinyTree, node = c(13, 18),
#'                                                only.leaf = TRUE)))))
#'
#' ## Visualize aggregated heatmap with tree
#' TreeHeatmap(tree = rowTree(tseagg), tree_fig = tree_fig,
#'             hm_data = SummarizedExperiment::assay(tseagg, "counts"))
#'
#' ## Cluster columns
#' TreeHeatmap(tree = rowTree(tseagg), tree_fig = tree_fig,
#'             hm_data = SummarizedExperiment::assay(tseagg, "counts"),
#'             cluster_column = TRUE)
#'
#' ## Split columns
#' col_split <- ifelse(colnames(tseagg) %in% paste0("S", seq_len(5)), "A", "B")
#' names(col_split) <- colnames(tseagg)
#' TreeHeatmap(tree = rowTree(tseagg), tree_fig = tree_fig,
#'             hm_data = SummarizedExperiment::assay(tseagg, "counts"),
#'             cluster_column = TRUE, column_split = col_split)
#'
#' ## Annotate columns
#' col_anno <- col_split
#' TreeHeatmap(tree = rowTree(tseagg), tree_fig = tree_fig,
#'             hm_data = SummarizedExperiment::assay(tseagg, "counts"),
#'             cluster_column = TRUE, column_split = col_split,
#'             column_anno = col_anno, column_anno_gap = 1)
#'
#' ## Change annotation colors
#' TreeHeatmap(tree = rowTree(tseagg), tree_fig = tree_fig,
#'             hm_data = SummarizedExperiment::assay(tseagg, "counts"),
#'             cluster_column = TRUE, column_split = col_split,
#'             column_anno = col_anno, column_anno_gap = 1,
#'             column_anno_color = c(A = "red", B = "blue"))
#'
#' ## Add column names
#' TreeHeatmap(tree = rowTree(tseagg), tree_fig = tree_fig,
#'             hm_data = SummarizedExperiment::assay(tseagg, "counts"),
#'             cluster_column = TRUE, column_split = col_split,
#'             column_anno = col_anno, column_anno_gap = 1,
#'             column_anno_color = c(A = "red", B = "blue"),
#'             show_colnames = TRUE, colnames_position = "bottom",
#'             colnames_angle = 90, colnames_size = 2,
#'             colnames_offset_y = -0.4)
#'
#' ## Add title
#' TreeHeatmap(tree = rowTree(tseagg), tree_fig = tree_fig,
#'             hm_data = SummarizedExperiment::assay(tseagg, "counts"),
#'             cluster_column = TRUE, column_split = col_split,
#'             column_anno = col_anno, column_anno_gap = 1,
#'             column_anno_color = c(A = "red", B = "blue"),
#'             show_colnames = TRUE, colnames_position = "bottom",
#'             colnames_angle = 90, colnames_size = 2,
#'             colnames_offset_y = -0.4,
#'             show_title = TRUE, title_offset_y = 2,
#'             title_color = "blue")
#'
#' ## Change colors
#' TreeHeatmap(tree = rowTree(tseagg), tree_fig = tree_fig,
#'             hm_data = SummarizedExperiment::assay(tseagg, "counts"),
#'             cluster_column = TRUE, column_split = col_split,
#'             column_anno = col_anno, column_anno_gap = 1,
#'             column_anno_color = c(A = "red", B = "blue"),
#'             show_colnames = TRUE, colnames_position = "bottom",
#'             colnames_angle = 90, colnames_size = 2,
#'             colnames_offset_y = -0.4,
#'             show_title = TRUE, title_offset_y = 2,
#'             title_color = "blue") +
#'             scale_fill_gradientn(
#'                 colours = c("blue", "yellow", "red"),
#'                 values = scales::rescale(c(5, 8, 10)),
#'                 guide = "colorbar", limits = c(5, 10))
#'
TreeHeatmap <- function(tree, tree_fig, hm_data,
                        tree_hm_gap = 0,
                        rel_width = 1,
                        cell_line_color = NA,
                        cell_line_size = 0,
                        column_order = NULL,
                        column_split = NULL,
                        column_split_gap = 0.2,
                        column_split_label = NULL,
                        split_label_fontface = "bold",
                        split_label_color = "black",
                        split_label_size = 3,
                        split_label_angle = 0,
                        split_label_offset_x = 0,
                        split_label_offset_y = 2,
                        split_label_hjust = 0.5,
                        split_label_vjust = 0,
                        column_anno = NULL,
                        column_anno_size = 1,
                        column_anno_color = NULL,
                        column_anno_gap = 0.1,
                        legend_title_hm = "Expression",
                        legend_title_column_anno = "group",
                        show_colnames = FALSE,
                        colnames_position = "top",
                        colnames_angle = 0,
                        colnames_offset_x = 0,
                        colnames_offset_y = 0,
                        colnames_size = 4,
                        colnames_hjust = 0.5,
                        show_rownames = FALSE,
                        rownames_position = "right",
                        rownames_angle = 0,
                        rownames_offset_x = 0,
                        rownames_offset_y = 0,
                        rownames_size = 4,
                        rownames_hjust = 0.5,
                        rownames_label = NULL,
                        show_title = FALSE,
                        title_hm = "First heatmap",
                        title_fontface = "bold",
                        title_color = "black",
                        title_size = 3,
                        title_angle = 0,
                        title_offset_x = 0,
                        title_offset_y = 2,
                        title_hjust = 0.5,
                        cluster_column = FALSE,
                        dist_method = "euclidean",
                        hclust_method = "ave",
                        show_row_tree = TRUE) {

    ## Check input arguments
    ## -------------------------------------------------------------------------
    ## Check that no leaf is covered multiple times in the matrix, as this
    ## would give a confusing heatmap
    lvs <- unlist(TreeSummarizedExperiment::findDescendant(
        tree, node = rownames(hm_data), only.leaf = TRUE, self.include = TRUE))
    if (any(duplicated(lvs))) {
        warning("Some leaves are contributing to multiple rows in hm_data. ",
                "This is likely unintended and may negatively impact the ",
                "heatmap visualization. Please check your data carefully.")
    }

    if (!is.null(column_order) && !all(column_order %in% colnames(hm_data))) {
        stop("column_order: Some columns are not present in hm_data.")
    }

    if (!is.null(column_split)) {
        if (is.null(names(column_split))) {
            stop("column_split: should be a named vector")
        }
        if (!all(names(column_split) %in% colnames(hm_data))) {
            stop("column_split: Some columns are not present in hm_data.")
        }
    }

    if (!is.null(column_anno_color) && is.null(names(column_anno_color))) {
        stop("column_anno_color: should be a named vector")
    }

    if (!is.null(rownames_label)) {
        if (is.null(names(rownames_label))) {
            stop("rownames_label: should be a named vector")
        }
        if (!all(names(rownames_label) %in% rownames(hm_data))) {
            stop("rownames_label: Some rows are not present in hm_data")
        }
    }

    ## Get data from tree figure
    ## -------------------------------------------------------------------------
    df <- tree_fig$data

    ## Heatmap - build up hm_df
    ## -------------------------------------------------------------------------
    ## Data
    hm_df <- data.frame(hm_data, check.names = FALSE)

    ## Nodes
    hm_df$node <- TreeSummarizedExperiment::convertNode(
        tree = tree, node = rownames(hm_df))

    ## Row labels
    if (is.null(rownames_label)) {
        hm_df$row_label <- rownames(hm_df)
    } else {
        hm_df$row_label <- rownames_label[rownames(hm_df)]
    }

    ## y position for each row
    desd_hm <- TreeSummarizedExperiment::findDescendant(
        tree = tree, node = hm_df$node, only.leaf = FALSE, self.include = TRUE)
    hm_df$y <- vapply(desd_hm, FUN = function(x) {
        xx <- match(x, df$node)
        y <- df$y[xx]
        ## Get the middle point of the range defined by the nodes
        mean(range(y, na.rm = TRUE))
    }, NA_real_)

    ## Height of each row
    hm_df$height <- vapply(desd_hm, FUN = function(x) {
        xx <- match(x, df$node)
        y <- df$y[xx]

        cy <- colnames(df)
        if ("scale" %in% cy) {
            dt <- unique(df$scale[xx])
            if (length(dt) > 1) {
                dt <- max(setdiff(dt, 1), 1)
            }
        } else {
            dt <- 1
        }
        ## Size of range spanned by nodes
        diff(range(y, na.rm = TRUE)) + dt
    }, NA_real_)

    ## Width of each column in the heatmap
    hm_df$width <- rel_width * (df$x |>
                                    range(na.rm = TRUE) |>
                                    diff()) / ncol(hm_data)

    ## Convert to long form
    ## -------------------------------------------------------------------------
    hm_dt <- tidyr::pivot_longer(hm_df, names_to = "variable",
                                 values_to = "value",
                                 cols = -c("node", "row_label", "y",
                                           "height", "width"))

    ## Determine column order
    ## -------------------------------------------------------------------------
    if (!is.null(column_split)) {
        # 1. column_split is given, ignore column order. The order within the
        #    same slice is determined by the input order of column split
        # 2. column_split is given and column_cluster = TRUE, the order within
        #    the same slice is determined by the column similarity
        column_split <- factor(column_split, levels = unique(column_split))
        if (cluster_column) {
            ## Column similarity within the same slice
            ind_split <- lapply(levels(column_split), FUN = function(x) {
                names(column_split)[column_split == x]
            })
            similar <- lapply(ind_split, FUN = function(x) {
                if (length(x) > 2) {
                    st <- hm_data[, x, drop = FALSE]
                    xx <- hclust(dist(t(st), method = dist_method),
                                 method = hclust_method)
                    colnames(st)[xx$order]
                } else {
                    x
                }
            })
            similar_order <- unlist(similar)
            column_split <- column_split[similar_order]
        }
        split_level <- sort(column_split)

        if (!is.null(column_order)) {
            warning("column_order is ignored when column_split is given")
        }
        column_order <- names(split_level)
    } else {
        ## column_split isn't given
        # 1. column_order is given, use column_order and ignore column
        #    similarity.
        # 2. column_order isn't given but allow cluster columns, order columns
        #    by similarity
        # 3. column_order isn't given and cluter columns is not allowed, use the
        #    original column order.
        split_level <- rep(0, ncol(hm_data))
        names(split_level) <- colnames(hm_data)
        split_level <- factor(split_level)
        if (is.null(column_order) && !cluster_column) {
            column_order <- colnames(hm_data)
        }
        if (!is.null(column_order) && cluster_column) {
            ## This needs to be before the next check, since that will
            ## create column_order by clustering and thus this would
            ## always trigger
            warning("cluster_column is ignored because column_order is given")
        }
        if (is.null(column_order) && cluster_column) {
            hc <- hclust(dist(t(hm_data), method = dist_method),
                         method = hclust_method)
            column_order <- colnames(hm_data)[hc$order]
        }
    }

    ## Get x coordinate and prepare for plotting
    ## -------------------------------------------------------------------------
    hm_dt <- hm_dt |>
        dplyr::mutate(variable = factor(.data$variable,
                                        levels = column_order)) |>
        dplyr::mutate(column_order = as.numeric(.data$variable) - 1) |>
        dplyr::mutate(split_level = split_level[.data$variable]) |>
        dplyr::mutate(split_level = as.numeric(.data$split_level) - 1) |>
        dplyr::mutate(x = max(df$x, na.rm = TRUE) +
                          tree_hm_gap + .data$width/2 +
                          .data$column_order * .data$width +
                          .data$split_level * column_split_gap) |>
        dplyr::select("node", "row_label", "x", "y",
                      "height", "width", "variable",
                      "value", "column_order", "split_level")

    ## Plot tree + heatmap
    ## -------------------------------------------------------------------------
    if (!show_row_tree) {
        tree_fig <- ggplot2::ggplot() + ggplot2::theme_void()
        hm_dt <- hm_dt |>
            dplyr::mutate(x = .data$x - min(.data$x))
    }
    p <- tree_fig +
        ggplot2::geom_tile(data = hm_dt,
                           aes(x = x, y = y, height = height,
                               fill = value, width = width),
                           color = cell_line_color,
                           linewidth = cell_line_size,
                           inherit.aes = FALSE) +
        ggplot2::labs(fill = legend_title_hm)

    if (is.numeric(hm_dt$value)) {
        p <- p + ggplot2::scale_fill_viridis_c()
    } else {
        p <- p + ggplot2::scale_fill_viridis_d()
    }

    ## Heatmap annotation
    ## -------------------------------------------------------------------------
    if (!is.null(column_anno)) {
        if (is.null(column_anno_color)) {
            anno_uc <- unique(column_anno)
            column_anno_color <- viridis::viridis(length(anno_uc))
            names(column_anno_color) <- anno_uc
        }
        ## column annotations
        anno_df <- hm_dt |>
            dplyr::select("variable", "x", "width") |>
            dplyr::distinct() |>
            dplyr::mutate(
                variable = as.character(.data$variable),
                x = .data$x - 0.5 * .data$width,
                xend = .data$x + .data$width,
                y = max(df$y, na.rm = TRUE) + column_anno_gap,
                yend = max(df$y, na.rm = TRUE) + column_anno_gap,
                anno_group = column_anno[.data$variable],
                anno_color = column_anno_color[.data$anno_group])
        anno_color <- anno_df$anno_color
        names(anno_color) <- anno_df$anno_group

        p <- p +
            ggnewscale::new_scale_color() +
            ggplot2::geom_segment(
                data = anno_df,
                aes(x = x, y = y, xend = xend, yend = yend, color = anno_group),
                inherit.aes = FALSE,
                linewidth = column_anno_size) +
            ggplot2::scale_color_manual(values = anno_color) +
            ggplot2::labs(color = legend_title_column_anno)
    } else {
        anno_df <- NULL
    }

    ## Heatmap column and row names
    ## -------------------------------------------------------------------------
    if (show_colnames) {
        cn_df <- hm_dt |>
            dplyr::select("variable", "x", "width") |>
            dplyr::distinct() |>
            dplyr::mutate(
                variable = as.character(.data$variable),
                y_top = max(hm_dt$y + 0.5 * hm_dt$height) + column_anno_gap,
                y_bottom = min(hm_dt$y - 0.5 * hm_dt$height),
                y = ifelse(colnames_position == "top",
                           .data$y_top, .data$y_bottom))
        p <- p + ggplot2::geom_text(data = cn_df,
                                    aes(x = x, y = y, label = variable),
                                    size = colnames_size, inherit.aes = FALSE,
                                    angle = colnames_angle,
                                    nudge_x = colnames_offset_x,
                                    nudge_y = colnames_offset_y,
                                    hjust = colnames_hjust)
    } else {
        cn_df <- NULL
    }

    if (show_rownames) {
        rn_df <- hm_dt |>
            dplyr::select("y", "width", "row_label") |>
            dplyr::distinct() |>
            dplyr::mutate(x_right = max(hm_dt$x + 0.5 * hm_dt$width),
                          x_left = min(hm_dt$x - 0.5 * hm_dt$width),
                          x = ifelse(rownames_position == "right",
                                     .data$x_right, .data$x_left))
        p <- p + ggplot2::geom_text(data = rn_df,
                                    aes(x = x, y = y, label = row_label),
                                    size = rownames_size, inherit.aes = FALSE,
                                    angle = rownames_angle,
                                    nudge_x = rownames_offset_x,
                                    nudge_y = rownames_offset_y,
                                    hjust = rownames_hjust)
    } else {
        rn_df <- NULL
    }

    ## Heatmap title
    ## -------------------------------------------------------------------------
    if (show_title) {
        title_df <- data.frame(x = mean(range(hm_dt$x, na.rm = TRUE)),
                               y = max(hm_dt$y),
                               label = title_hm)
        p <- p + ggplot2::geom_text(data = title_df,
                                    aes(x = x, y = y, label = label),
                                    inherit.aes = FALSE,
                                    fontface = title_fontface,
                                    colour = title_color,
                                    size = title_size,
                                    angle = title_angle,
                                    nudge_x = title_offset_x,
                                    nudge_y = title_offset_y,
                                    hjust = title_hjust)
    } else {
        title_df <- NULL
    }

    ## Split title
    ## -------------------------------------------------------------------------
    split_df <- hm_dt |>
        dplyr::select("x", "y", "split_level") |>
        dplyr::group_by(.data$split_level) |>
        dplyr::summarise(x = mean(range(.data$x, na.rm = TRUE)),
                         y = max(range(.data$y, na.rm = TRUE))) |>
        dplyr::arrange(.data$split_level) |>
        dplyr::mutate(column_split = levels(column_split))
    if (!is.null(column_split_label)) {
        split_df <- split_df |>
            dplyr::mutate(
                split_label =
                    column_split_label[as.character(.data$column_split)])
        p <- p + ggplot2::geom_text(data = split_df,
                                    aes(x = x, y = y, label = split_label),
                                    inherit.aes = FALSE,
                                    fontface = split_label_fontface,
                                    colour = split_label_color,
                                    size = split_label_size,
                                    angle = split_label_angle,
                                    nudge_x = split_label_offset_x,
                                    nudge_y = split_label_offset_y,
                                    hjust = split_label_hjust,
                                    vjust = split_label_vjust)

    }

    ## Add data for later possible retrieval
    ## -------------------------------------------------------------------------
    p$temp_data <- list(
        hm_data = hm_dt,
        row_name = rn_df,
        column_name = cn_df,
        hm_title = title_df,
        column_anno = anno_df,
        column_order = column_order,
        column_split = split_df)

    return(p)
}
