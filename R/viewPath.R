#' Visualize the tree paths
#' viewPath shows all paths on the tree. Each connects a leaf to the root. Users
#' could decide what to display on the x and y axis. The data show on x and y
#' axis should be provided as a data frame.
#'
#' @param tree A phylo object.
#' @param data_frame A data frame. It should include at least a column of node
#'   (number), and a column of score values that is to show on the y-axis of the
#'   figure.
#' @param node_column A column name in \code{data_frame}. The column stores the
#'   node number.
#' @param y_axis_column A column name in \code{data_frame}. The column will be
#'   shown on the y axis of the figure.
#' @param x_axis_column A column name in \code{data_frame}. The column will be
#'   shown on the y-axis of the figure. If NULL (default), the number of
#'   descendant leaves is automatically generated and used.
#' @param y_lab A character to be used as the label of y-axis.
#' @param x_lab A character to be used as the label of x-axis.
#' @param node_true A character or numeric vector. It gives the truth: the nodes
#'   with signal (with different levels under different conditions).
#' @param lab_true A character vector. It has the same length as
#'   \code{node_true} and provides the information about the signal direction
#'   (e.g. up or down). 
#' @param col_true A vector of colors. It should be named using the elements of
#'   \code{lab_true}.
#' @param node_found A character or numeric vector. It gives the nodes that are
#'   found.
#' @param legend_title_found A legend title of the found.
#' @param legend_title_true A legend title of the truth.
#' @param only_true_path A logical value, TRUE or FALSE. If TRUE, only paths
#'   with signal are shown. The default is FALSE
#' @param point_size A numeric value to specify the point size
#' @param point_stroke A numeric value to specify the point stroke
#' @param point_shape A numeric value to specify the point shape
#' @param log_x A logic value, TRUE or FALSE. It decides whether x-axis should
#'   be displayed in log scale
#' @param log_y A logic value, TRUE or FALSE. It decides whether y-axis should
#'   be displayed in log scale
#' 
#' @importFrom dplyr mutate "%>%" filter
#' @export
#' @author Ruizhu Huang
#' @examples 
#'   
#' library(TreeSummarizedExperiment)
#' library(ggtree)
#' library(ggplot2)
#' data("tinyTree")
#' ggtree(tinyTree, branch.length = "none") + 
#'     geom_text2(aes(label = node)) +
#'     geom_hilight(node = 13, fill = "orange", alpha = 0.4) +
#'     geom_hilight(node = 18, fill = "blue", alpha = 0.4)
#' 
#' 
#' set.seed(1)
#' sv <- c(0.8, 0.83, 0.86, -0.8, -0.85, 
#'         runif(7, -0.5, 0.5), 0.95, 0.9,
#'         runif(3, -0.7, 0), -0.9, 0.3)
#' pair <- cbind(11, 1:19)
#' dist <- apply(pair, 1, FUN= function(x) {
#'                 distNode(tree = tinyTree, node = x)
#'                 })
#' dat <- data.frame(node = 1:19, score = round(sv, 2),
#'                   dist = dist)
#' 
#' # score in red, node in blue
#' ggtree(tinyTree, branch.length = "none") %<+% dat + 
#'     geom_text2(aes(label = node), color = "blue", vjust = 1) +
#'     geom_text2(aes(label = score), color = "red", vjust = -0.5) +
#'     geom_hilight(node = 13, fill = "orange", alpha = 0.4) +
#'     geom_hilight(node = 18, fill = "blue", alpha = 0.4)
#' 
#' # The score (red textS in figure above) is now in y-axis,
#' # the x-axis is the number of descendant leaves of the corresponding node
#' viewPath(tree = tinyTree, data_frame = dat,
#'          node_column = "node", y_axis_column = "score", 
#'          node_true = c(13, 18), lab_true = c("up", "down"),
#'          col_true = c("up" = "orange", "down" = "blue"), 
#'          node_found = c(13, 18), log_x = FALSE, log_y = FALSE)
#' 
#' # the x-axis is the distance between the node and the root         
#' viewPath(tree = tinyTree, data_frame = dat,
#'          node_column = "node", y_axis_column = "score", 
#'          x_axis_column = "dist",
#'          node_true = c(13, 18), lab_true = c("up", "down"),
#'          col_true = c("up" = "orange", "down" = "blue"), 
#'          node_found = c(13, 18), log_x = FALSE, log_y = FALSE)
#' 

viewPath <- function(tree,
                     data_frame, node_column,
                     y_axis_column, 
                     x_axis_column = NULL,
                     y_lab = NULL, x_lab = NULL,
                     node_true, lab_true, col_true,
                     node_found, 
                     legend_title_found = "Found", 
                     legend_title_true = "Truth",
                     only_true_path = FALSE,
                     point_size = 2,
                     point_stroke = 2,
                     point_shape = 21,
                     log_x = TRUE, log_y = TRUE){
    
    # ================== check inputs ====================
    nam <- names(col_true)
    labU <- unique(lab_true)
    if (!setequal(nam, labU)) {
        stop("The name of col_true should come from lab_true")
    }
    if ("grey" %in% col_true) {
        stop("Please choose colors other than grey or white for col_true")
    }
    
    nodeNum <- data_frame[[node_column]]
    ## the number of descendants
    des <- findOS(tree = tree, node = nodeNum, only.leaf = TRUE)
    len <- lapply(des, length)
    
    ## Categorize the branches: up, down and other.
    trueAll <- findOS(tree = tree, node = node_true, only.leaf = FALSE,
                   self.include = TRUE)
    truth <- rep("other", nrow(data_frame))
    for (i in seq_along(node_true)) {
        truth[nodeNum %in% trueAll[[i]]] <- lab_true[i]
    }
    found <- rep(FALSE, nrow(data_frame))
    found[nodeNum %in% node_found] <- TRUE
    
    ## Add columns to the data
    data_new <- data_frame %>% 
        mutate(truth = truth) %>% 
        mutate(found = found)
    
    if (!is.null(x_axis_column)) {
        data_new <- data_new %>%
            mutate(x_axis = (!!as.name(x_axis_column)))
    } else {
        data_new <- data_new %>%
            mutate(x_axis = unlist(len))
    }
    ## ================= data of edges ========================
    mt <- matTree(tree)
    edg <- data.frame(tree$edge)
    colnames(edg) <- c("parent", "child")
    edg <- edg %>%
        mutate(xstart = rep(NA, nrow(edg))) %>%
        mutate(ystart = rep(NA, nrow(edg))) %>%
        mutate(xend = rep(NA, nrow(edg))) %>%
        mutate(yend = rep(NA, nrow(edg))) %>%
        mutate(truth = rep(NA, nrow(edg))) 
    
    for (i in seq_len(nrow(edg))) {
        ci <- which(data_new[[node_column]] %in% edg[i, "child"])
        pi <- which(data_new[[node_column]] %in% edg[i, "parent"])
        edg[i, c("xstart", "ystart")] <- data_new[ci, c("x_axis", y_axis_column)]
        edg[i, c("xend", "yend")] <- data_new[pi, c("x_axis", y_axis_column)]
        truth.i <- unique(data_new[c(ci, pi), "truth"])
        edg[i, "truth"] <- ifelse(length(truth.i) == 1, 
                                  truth.i, "other")
        edg[i, "found_child"] <- data_new[ci, "found"]
        edg[i, "found_parent"] <- data_new[pi, "found"]
        edg[i, "truth_child"] <- data_new[ci, "truth"]
        edg[i, "truth_parent"] <- data_new[pi, "truth"]
    }
    
    # ================= only show signal path ?=======================
    if (only_true_path) {
        leaf_true <- findOS(tree = tree, node = node_true, 
                        only.leaf = TRUE, self.include = TRUE)
        leaf_true <- unlist(leaf_true)
        sel <- mt[mt[, 1] %in% leaf_true, ]
        sel <- unique(as.vector(sel))
        keep <- sel[!is.na(sel)]
        edg <- edg %>%
            filter(parent %in% keep & child %in% keep)
    }
    
    # =========================  figure ============================ 
    col_true <- c(col_true, "other" = "grey")
    col_fill <- c(col_true, "other" = "white")
    if (is.null(x_lab)) { x_lab <- x_axis_column }
    if (is.null(x_lab)) {x_lab <- "The number of descendant leaves"}
    if (is.null(y_lab)) { y_lab <- y_axis_column }
    p <- ggplot() +
        geom_segment(data = edg, aes(x = xstart, y = ystart, 
                         xend = xend, yend = yend, color = truth)) +
        geom_point(data = edg[edg$found_child, ], 
                   aes(x = xstart, y = ystart, fill = truth_child),
                   color = "black", size = point_size, 
                   stroke = point_stroke, shape = point_shape) +
        geom_point(data = edg[edg$found_parent, ],
                   aes(x = xend, y = yend, fill = truth_parent),
                   color = "black", size = point_size,
                   stroke = point_stroke, shape = point_shape) +
        scale_fill_manual(values = col_fill) +
        scale_color_manual(values = col_true) +
        xlab(x_lab) +
        ylab(y_lab) +
        labs(fill = legend_title_found, color = legend_title_true)
    
    if (log_x) {
        p <- p + scale_x_log10()
    } 
    if (log_y) {
        p <- p + scale_y_log10() 
    }
    print(p)
    
    
}
