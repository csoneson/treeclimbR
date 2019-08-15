#' view tree and its two zoomed branches
#' 
#' \code{viewBranch} shows a full tree and two zommed branches
#' 
#' @param tree A phylo object
#' @param ann_data A data frame to annotate nodes on the tree. Default is NULL,
#'   otherwise, it should include at least two columns: one column about nodes
#'   and named as \code{node} and the other column about the annotation value
#'   and could be named freely.
#' @param ann_column A character vector to specify columns in \code{ann_data}
#'   that provide values to annotate nodes.
#' @param ann_color A character vector to specify the colors to be used for the
#'   annotation. It should have the same length as \code{ann_column}.
#' @param ann_hjust A numeric vector. It is to adjust the label locations
#'   horizontally. Please refer to \code{hjust} in \code{ggtree::geom_text2}
#' @param ann_vjust A numeric vector. It is to adjust the label locations
#'   vertically. Please refer to \code{vjust} in \code{ggtree::geom_text2}
#' @param ann_size A numeric vector. It is to adjust the label size. Please
#'   refer to \code{size} in \code{ggtree::geom_text2}
#' @param hlight_node A vector of node number. The branch to be highlight.
#'   Please refer to \code{node} in \code{ggtree::geom_hilight}
#' @param hlight_alpha A numeric value. Please refer to \code{alpha} in
#'   \code{ggtree::geom_hilight}
#' @param view_node A vector of node number. Branches to be shown in subplots.
#' @param zoom_node A vector of node number. Branches to be zoomed in. Please
#'   refer to \code{node} in \code{ggtree::scaleClade}
#' @param group_leaf A named list of leaf node number. Please refer to
#'   \code{node} in \code{ggtree::groupOTU}
#' @param group_color A named vector provide the color for \code{group_leaf}.
#'   The names should match with the name in \code{group_leaf} and one extra
#'   with name "0". Please see the example.
#' @param zoom_scale A numeric vector. The zoom scale for \code{zoom_node}.
#'   Please refer to \code{scale} in \code{ggtree::scaleClade}
#' @param branch_length Please refer to the argument \code{branch.length} of
#'   \code{\link[ggtree]{ggtree}}. The default is "none".
#' @param layout Please refer to the argument \code{layout} of
#'   \code{\link[ggtree]{ggtree}}. The default is "rectangular".
#' @param edge_size Please refer to the argument \code{size} of
#'   \code{\link[ggtree]{ggtree}}. The default is 1.
#' @param rel_widths The widths of three plots. Please refer to
#'   \code{rel_widths} of \code{\link[cowplot]{plot_grid}}.
#' @param labels The labels of three plots. Please refer to \code{labels} of
#'   \code{\link[cowplot]{plot_grid}}
#' @param ... More arguments accepted by \code{\link[cowplot]{plot_grid}}.
#' 
#' @importFrom ggtree geom_hilight geom_text2 viewClade "%<+%" geom_point2
#' @importFrom cowplot plot_grid 
#' @export
#' 
#' @return a figure
#' 
#' @examples 
#' library(TreeSummarizedExperiment)
#' data("tinyTree")
#' 
#' df <- data.frame(node = printNode(tinyTree, type = "all")$nodeNum,
#'                  score1 = round(runif(19), 2),
#'                  score2 = round(runif(19), 2))   
#'                  
#' grp <- list(A = 1:3, B = 7:8)                  
#' viewBranch(tree = tinyTree, ann_data = df,
#'            ann_column = c("node", "score2"),
#'            ann_color = c("blue", "red"),
#'            ann_hjust = c(0.2, 0.2),
#'            ann_vjust = c(0.6, -0.5),
#'            ann_size = c(3, 3), 
#'            view_node = c(16, 14), nrow = 1,
#'            zoom_node = c(16), zoom_scale = 10, 
#'            group_leaf = grp,
#'            group_color = c(A ="orange", B = "blue", "0" = "grey"))

viewBranch <- function(tree, 
                       ann_data = NULL,
                       ann_column = NULL, 
                       ann_color = "black", 
                       ann_hjust= c(0.2, .2), 
                       ann_vjust = c(.6, -0.5),
                       # ann_hjust, 
                       # ann_vjust,
                       ann_size = 1, 
                       hlight_node = NULL, 
                       hlight_alpha = 0.5,
                       hlight_fill = "blue",
                       point_node = NULL,
                       point_color = "cyan",
                       point_size = 3,
                       view_node = NULL, 
                       zoom_node = NULL,
                       zoom_scale = 2,
                       group_leaf = NULL,
                       group_color = "orange",
                       nrow = nrow,
                       rel_widths = c(0.5, 2, 3),
                       branch_length = "none",
                       layout = "rectangular",
                       edge_size = 1,
                       #                       rel_widths, 
                       ...) {
    
    # ================= check inputs  =================
    # text annotation
    if (!is.null(ann_column)) {
        lcn <- length(ann_column)
        lcl <- length(ann_color)
        lhs <- length(ann_size)
        
        if (lcn != lcl) {
            if (length(ann_color) == 1) {
                ann_color <- rep(ann_color, length(ann_column))
            } else {
                stop("ann_color has different length to ann_column.")
            }
        }
        if (lcn != lhs) {
            if (length(ann_size) == 1) {
                ann_size <- rep(ann_size, length(ann_column))
            } else {
                stop("ann_size has different length to ann_column.")
            }
        }
    }
    
    ## highlight
    if (!is.null(hlight_node)) {
        hn <- length(hlight_node)
        hc <- length(hlight_fill)
        if (hn != hc) {
            if (length(hlight_fill) == 1) {
                hlight_fill <- rep(hlight_fill, length(hlight_node))
            } else {
                stop("hlight_node has different length to hlight_fill.")
            }
        }
    }
    
    # zoom in
    if (!is.null(zoom_node)) {
        zn <- length(zoom_node)
        zs <- length(zoom_scale)
        if (zn != zs) {
            if (zs == 1) {
                zoom_scale <- rep(zoom_scale, length(zoom_node))
            } else {
                stop("zoom_node has different length to zoom_scale.")
            }
        }
    }
    
    # points
    if (!is.null(point_node)) {
        pn <- length(point_node)
        cn <- length(point_color)
        if (pn != cn) {
            if (cn == 1) {
                point_color <- rep(point_color, length(point_node))
            } else {
                stop("point_node has different length to point_color.")
            }
            
        } 
        point_color <- point_color[order(point_node)]
        # if (cn != pn) {
        #     if (cn == 1) {
        #         point_color <- rep(point_color, pn)
        #     } else {
        #         stop("ann_color has different length to ann_column.")
        #     }
        # }
    }
    
    # ===================================================================
    p0 <- ggtree(tree, branch.length = branch_length, layout = layout,
                 size = edge_size) 
    
    # scale branches
    if (!is.null(zoom_node)) {
        for (i in seq_along(zoom_node)) {
            p0 <- scaleClade(p0, zoom_node[i], scale = zoom_scale[i])
        }
    }
    
    # group the branches
    if (!is.null(group_leaf)) {
        p0 <- groupOTU(p0, group_leaf, "grp") + aes(color = grp) +
            scale_color_manual(values = group_color)
    }
    
    # add annotate data
    if (!is.null(ann_data)) {
        p0 <- p0 %<+% ann_data 
    }
    
    # hilight branches
    if (!is.null(hlight_node)) {
        for (i in seq_along(hlight_node)) {
            p0 <- p0 +
                geom_hilight(node = hlight_node[i], 
                             alpha = hlight_alpha,
                             fill = hlight_fill[i])
        }
    }
    
    # # add points
    if (!is.null(point_node)) {
        
        p0 <- p0 +
            geom_point2(aes(subset = (node %in% point_node)),
                        color = point_color, size = point_size)
    }
    # view branches
    if (!is.null(view_node)) {
        plist <- vector("list", length(view_node))
        for (i in seq_along(view_node)) {
            plist[[i]] <- viewClade(p0, node = view_node[i])
            
        }
    } else {
        plist <- NULL
    }
    
    # add texts
    if (!is.null(ann_column)) {
        
        for (i in seq_along(ann_column)) {
            if (is.null(plist)) {
                p0 <- p0 + geom_text2(aes_string(label = ann_column[i]),
                                      hjust = ann_hjust[i],
                                      vjust = ann_vjust[i],
                                      color = ann_color[i],
                                      size = ann_size[i])
            } else {
                plist <- lapply(plist, FUN = function(x) {
                    xx <- x + geom_text2(aes_string(label = ann_column[i]),
                                         hjust = ann_hjust[i],
                                         vjust = ann_vjust[i],
                                         color = ann_color[i],
                                         size = ann_size[i])
                    return(xx)
                })  
            }
            
        }
    }
    
    outlist <- c(list(p0), plist)
    if (missing(view_node)) {
        p <- outlist[[1]]
    } else {
        p <- plot_grid(plotlist = outlist, nrow = 1,
                       rel_widths = rel_widths, ...)
        
    }
    p
}
