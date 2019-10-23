#' visualize simulated scenario
#'
#' \code{viewSim} is to visualize the output from the function \code{simData}.
#'
#' @param obj The output from \code{simData}
#' @param zoom_scale A positive numeric value. If it is above one, branches with
#'   fold change equal to one (non-signal branch) will be zoomed in; If below
#'   one, they will be shrinked. Default is 2
#' @param ref_value Default is 1. The fold change of the non-signal branch. (The
#'   slot \code{FC} in \code{metadata} of \strong(obj))
#' @param legend_position The legend position. see \code{\link[ggplot2]{theme}}
#' @param show_leaf TRUE or FALSE with fold change above or below 1 will be
#'   labelled.
#' @param ... Additional arguments that are accepted by
#'   \code{\link[ggtree]{ggtree}}.
#'
#' @importFrom S4Vectors metadata
#' @importFrom ggtree ggtree scaleClade geom_tiplab groupOTU
#' @import ggplot2
#' @export
#'
#' @return a figure
#' @examples
#' set.seed(1)
#' y <- matrix(rnbinom(100,size=1,mu=10),nrow=10)
#' colnames(y) <- paste("S", 1:10, sep = "")
#' rownames(y) <- tinyTree$tip.label
#'
#'
#' toy_lse <- TreeSummarizedExperiment(rowTree = tinyTree,
#'                                     assays = list(y))
#' res <- parEstimate(obj = toy_lse)
#'
#' set.seed(1122)
#' dat1 <- simData(obj = res)
#' viewSim(obj = dat1, layout = "circular")
#'
#'
#'
viewSim <- function(obj, zoom_scale = 2, legend_position = "left", 
                    show_leaf = FALSE, ref_value = 1, ...){
    
    md <- metadata(obj)
    # tree
    tree <- rowTree(obj)
    
    # branch
    branch <- c(md$branch$A, md$branch$B)
    if (!is.numeric(branch)) {
        branch <- transNode(tree = tree, node = branch, 
                            use.alias = TRUE, message = FALSE)
    }
    
    # scenario
    sc <- md$scenario
    
    # fc
    fc <- md$FC 
    
    # leaves that change abundance
    high <- names(fc)[fc > ref_value]
    low <- names(fc)[fc < ref_value]
    cls <- list(increase = high, decrease = low)
    
    # figure
    tree1 <- groupOTU(tree, cls)
    p1 <- ggtree(tree1, aes(color = group), ...) +
        scale_color_manual(values = c(increase = "orange", 
                                      decrease = "blue", "0" = "grey"),
                           labels = c(increase = "increase", 
                                      decease = "decrease", "0" = "other")) +
        theme(legend.position = legend_position,
              legend.background = element_rect(fill="transparent")) +
        guides(color = guide_legend(title = sc))
    
    p2 <- scaleClade(p1, node = branch[1], scale = zoom_scale)
    p3 <- scaleClade(p2, node = branch[2], scale = zoom_scale)
    if (show_leaf) {
      p3 <- p3 + geom_tiplab(aes(label = label))
    }
    
    fc_data <- data.frame(node = transNode(tree = tree, node = names(fc)),
                          fc = fc)
    fc_data <- fc_data[fc_data$fc != ref_value, ]
    if (sc == "S2") {
        p3 <- p3 %<+% fc_data + geom_point2(aes(size = fc), alpha = 0.3)
    }  
   print(p3) 
    
}

