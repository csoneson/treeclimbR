test_that("aggDS works", {
    ## Generate example data
    library(ape)
    tr <- rtree(3, tip.label = LETTERS[seq_len(3)])
    set.seed(1L)
    cc <- matrix(rpois(60, 10), nrow = 6)
    rownames(cc) <- paste0("gene", seq_len(6))
    colnames(cc) <- paste0("cell", seq_len(10))
    cd <- data.frame(sid = rep(1:2, each = 5),
                     gid = rep(letters[1:2], each = 5),
                     cid = sample(LETTERS[1:3], size = 10, replace = TRUE),
                     stringsAsFactors = FALSE)
    expect_warning({
        tse <- TreeSummarizedExperiment(assays = list(counts = cc),
                                        colTree = tr,
                                        colNodeLab = cd$cid,
                                        colData = cd)},
        "Multiple nodes are found to have the same label"
    )

    ## Test that function returns error for misspecified input
    ## -------------------------------------------------------------------------
    ## TSE
    expect_error(aggDS(TSE = 1, assay = "counts", sample_id = "sid",
                       group_id = "gid", cluster_id = "cid", FUN = sum,
                       message = FALSE),
                 "'TSE' must be of class 'TreeSummarizedExperiment'")
    tse2 <- tse
    TreeSummarizedExperiment::colTree(tse2) <- NULL
    expect_error(aggDS(TSE = tse2, assay = "counts", sample_id = "sid",
                       group_id = "gid", cluster_id = "cid", FUN = sum,
                       message = FALSE),
                 "'TreeSummarizedExperiment::colTreeTSE' must not be NULL")
    tse2 <- tse
    tse2$cid <- paste0(tse2$cid, "_1")
    expect_error(aggDS(TSE = tse2, assay = "counts", sample_id = "sid",
                       group_id = "gid", cluster_id = "cid", FUN = sum,
                       message = FALSE),
                 "[[cluster_id]] %in% TreeSummarizedExperiment::colTree(TSE)$tip.label) is not TRUE",
                 fixed = TRUE)
    ## assay
    expect_error(aggDS(TSE = tse, assay = c("counts", "counts"),
                       sample_id = "sid", group_id = "gid", cluster_id = "cid",
                       FUN = sum, message = FALSE),
                 "length(assay) == 1 is not TRUE", fixed = TRUE)
    expect_error(aggDS(TSE = tse, assay = "missing", sample_id = "sid",
                       group_id = "gid", cluster_id = "cid", FUN = sum,
                       message = FALSE),
                 "assay %in% SummarizedExperiment::assayNames(TSE)",
                 fixed = TRUE)
    expect_error(aggDS(TSE = tse, assay = TRUE, sample_id = "sid",
                       group_id = "gid", cluster_id = "cid", FUN = sum,
                       message = FALSE),
                 "is.numeric(assay)",
                 fixed = TRUE)
    expect_error(aggDS(TSE = tse, assay = 42, sample_id = "sid",
                       group_id = "gid", cluster_id = "cid", FUN = sum,
                       message = FALSE),
                 "assay %in% seq_len(length(SummarizedExperiment::assays(TSE)))",
                 fixed = TRUE)
    ## sample_id
    expect_error(aggDS(TSE = tse, assay = "counts", sample_id = 1,
                       group_id = "gid", cluster_id = "cid", FUN = sum,
                       message = FALSE),
                 "'sample_id' must be of class 'character'")
    expect_error(aggDS(TSE = tse, assay = "counts", sample_id = c("sid", "gid"),
                       group_id = "gid", cluster_id = "cid", FUN = sum,
                       message = FALSE),
                 "'sample_id' must have length 1")
    expect_error(aggDS(TSE = tse, assay = "counts", sample_id = "missing",
                       group_id = "gid", cluster_id = "cid", FUN = sum,
                       message = FALSE),
                 "all(c(sample_id, group_id, cluster_id) %in% colnames(",
                 fixed = TRUE)
    ## group_id
    expect_error(aggDS(TSE = tse, assay = "counts", sample_id = "sid",
                       group_id = 1, cluster_id = "cid", FUN = sum,
                       message = FALSE),
                 "'group_id' must be of class 'character'")
    expect_error(aggDS(TSE = tse, assay = "counts", sample_id = "sid",
                       group_id = c("gid", "cid"), cluster_id = "cid", FUN = sum,
                       message = FALSE),
                 "'group_id' must have length 1")
    expect_error(aggDS(TSE = tse, assay = "counts", sample_id = "sid",
                       group_id = "missing", cluster_id = "cid", FUN = sum,
                       message = FALSE),
                 "all(c(sample_id, group_id, cluster_id) %in% colnames(",
                 fixed = TRUE)
    ## cluster_id
    expect_error(aggDS(TSE = tse, assay = "counts", sample_id = "sid",
                       group_id = "gid", cluster_id = 1, FUN = sum,
                       message = FALSE),
                 "'cluster_id' must be of class 'character'")
    expect_error(aggDS(TSE = tse, assay = "counts", sample_id = "sid",
                       group_id = "gid", cluster_id = c("cid", "sid"), FUN = sum,
                       message = FALSE),
                 "'cluster_id' must have length 1")
    expect_error(aggDS(TSE = tse, assay = "counts", sample_id = "sid",
                       group_id = "gid", cluster_id = "missing", FUN = sum,
                       message = FALSE),
                 "all(c(sample_id, group_id, cluster_id) %in% colnames(",
                 fixed = TRUE)
    ## FUN
    expect_error(aggDS(TSE = tse, assay = "counts", sample_id = "sid",
                       group_id = "gid", cluster_id = "cid", FUN = 1,
                       message = FALSE),
                 "'FUN' must be of class 'function'")
    expect_error(aggDS(TSE = tse, assay = "counts", sample_id = "sid",
                       group_id = "gid", cluster_id = "cid", FUN = "sum",
                       message = FALSE),
                 "'FUN' must be of class 'function'")
    ## message
    expect_error(aggDS(TSE = tse, assay = "counts", sample_id = "sid",
                       group_id = "gid", cluster_id = "cid", FUN = sum,
                       message = 1),
                 "'message' must be of class 'logical'")
    expect_error(aggDS(TSE = tse, assay = "counts", sample_id = "sid",
                       group_id = "gid", cluster_id = "cid", FUN = sum,
                       message = c(TRUE, FALSE)),
                 "'message' must have length 1")

    ## Test that aggregation works - sum
    ## -------------------------------------------------------------------------
    expect_warning({
        out <- aggDS(TSE = tse, assay = "counts", sample_id = "sid",
                     group_id = "gid", cluster_id = "cid", FUN = sum,
                     message = FALSE)},
        "Multiple nodes are found to have the same label"
    )
    expect_s4_class(out, "SummarizedExperiment")
    expect_equal(nrow(out), nrow(tse))
    expect_equal(ncol(out), length(unique(tse$sid)))
    expect_equal(length(SummarizedExperiment::assays(out)), 5L) ## 3 leaves + 2 nodes
    expect_equal(rownames(out), rownames(tse))
    expect_equal(colnames(out), as.character(unique(tse$sid)))
    expect_equal(out$gid, c("a", "b"))
    expect_equal(SummarizedExperiment::assayNames(out),
                 paste0("alias_", seq_len(5)))
    expect_s3_class(S4Vectors::metadata(out)$experiment_info, "data.frame")
    expect_equal(nrow(S4Vectors::metadata(out)$experiment_info),
                 length(unique(tse$sid)))
    expect_equal(ncol(S4Vectors::metadata(out)$experiment_info), 3L)
    expect_equal(colnames(S4Vectors::metadata(out)$experiment_info),
                 c("sample_id", "group_id", "n_cells"))
    expect_equal(S4Vectors::metadata(out)$experiment_info$sample_id,
                 c(1, 2))
    expect_equal(S4Vectors::metadata(out)$experiment_info$group_id,
                 c("a", "b"))
    expect_equal(as.vector(S4Vectors::metadata(out)$experiment_info$n_cells),
                 c(5, 5))
    ## Check aggregated counts
    get_cells <- function(node_alias) {
        descs <- TreeSummarizedExperiment::findDescendant(colTree(tse), node_alias,
                                                          only.leaf = TRUE)[[1]]
        if (is.null(descs)) {
            ## leaf
            cells <- which(colLinks(tse)$nodeLab_alias == node_alias)
        } else {
            ## internal node
            cells <- which(colLinks(tse)$nodeNum %in% descs)
        }
        return(cells)
    }
    expect_equal(SummarizedExperiment::assay(out, "alias_1")[, 1],
                 rowSums(SummarizedExperiment::assay(tse, "counts")[, intersect(
                     get_cells("alias_1"),
                     which(tse$sid == 1)), drop = FALSE]))
    expect_equal(SummarizedExperiment::assay(out, "alias_1")[, 2],
                 rowSums(SummarizedExperiment::assay(tse, "counts")[, intersect(
                     get_cells("alias_1"),
                     which(tse$sid == 2)), drop = FALSE]))
    expect_equal(SummarizedExperiment::assay(out, "alias_2")[, 1],
                 rowSums(SummarizedExperiment::assay(tse, "counts")[, intersect(
                     get_cells("alias_2"),
                     which(tse$sid == 1)), drop = FALSE]))
    expect_equal(SummarizedExperiment::assay(out, "alias_2")[, 2],
                 rowSums(SummarizedExperiment::assay(tse, "counts")[, intersect(
                     get_cells("alias_2"),
                     which(tse$sid == 2)), drop = FALSE]))
    expect_equal(SummarizedExperiment::assay(out, "alias_3")[, 1],
                 rowSums(SummarizedExperiment::assay(tse, "counts")[, intersect(
                     get_cells("alias_3"),
                     which(tse$sid == 1)), drop = FALSE]))
    expect_equal(SummarizedExperiment::assay(out, "alias_3")[, 2],
                 rowSums(SummarizedExperiment::assay(tse, "counts")[, intersect(
                     get_cells("alias_3"),
                     which(tse$sid == 2)), drop = FALSE]))
    expect_equal(SummarizedExperiment::assay(out, "alias_4")[, 1],
                 rowSums(SummarizedExperiment::assay(tse, "counts")[, intersect(
                     get_cells("alias_4"),
                     which(tse$sid == 1)), drop = FALSE]))
    expect_equal(SummarizedExperiment::assay(out, "alias_4")[, 2],
                 rowSums(SummarizedExperiment::assay(tse, "counts")[, intersect(
                     get_cells("alias_4"),
                     which(tse$sid == 2)), drop = FALSE]))
    expect_equal(SummarizedExperiment::assay(out, "alias_5")[, 1],
                 rowSums(SummarizedExperiment::assay(tse, "counts")[, intersect(
                     get_cells("alias_5"),
                     which(tse$sid == 1)), drop = FALSE]))
    expect_equal(SummarizedExperiment::assay(out, "alias_5")[, 2],
                 rowSums(SummarizedExperiment::assay(tse, "counts")[, intersect(
                     get_cells("alias_5"),
                     which(tse$sid == 2)), drop = FALSE]))

    ## Test that aggregation works - mean
    ## -------------------------------------------------------------------------
    expect_warning({
        out <- aggDS(TSE = tse, assay = "counts", sample_id = "sid",
                     group_id = "gid", cluster_id = "cid", FUN = mean,
                     message = TRUE)},
        "Multiple nodes are found to have the same label"
    )
    expect_s4_class(out, "SummarizedExperiment")
    expect_equal(nrow(out), nrow(tse))
    expect_equal(ncol(out), length(unique(tse$sid)))
    expect_equal(length(SummarizedExperiment::assays(out)), 5L) ## 3 leaves + 2 nodes
    expect_equal(rownames(out), rownames(tse))
    expect_equal(colnames(out), as.character(unique(tse$sid)))
    expect_equal(out$gid, c("a", "b"))
    expect_equal(SummarizedExperiment::assayNames(out),
                 paste0("alias_", seq_len(5)))
    expect_s3_class(S4Vectors::metadata(out)$experiment_info, "data.frame")
    expect_equal(nrow(S4Vectors::metadata(out)$experiment_info),
                 length(unique(tse$sid)))
    expect_equal(ncol(S4Vectors::metadata(out)$experiment_info), 3L)
    expect_equal(colnames(S4Vectors::metadata(out)$experiment_info),
                 c("sample_id", "group_id", "n_cells"))
    expect_equal(S4Vectors::metadata(out)$experiment_info$sample_id,
                 c(1, 2))
    expect_equal(S4Vectors::metadata(out)$experiment_info$group_id,
                 c("a", "b"))
    expect_equal(as.vector(S4Vectors::metadata(out)$experiment_info$n_cells),
                 c(5, 5))
    ## Check aggregated counts
    get_cells <- function(node_alias) {
        descs <- TreeSummarizedExperiment::findDescendant(colTree(tse), node_alias,
                                                          only.leaf = TRUE)[[1]]
        if (is.null(descs)) {
            ## leaf
            cells <- which(colLinks(tse)$nodeLab_alias == node_alias)
        } else {
            ## internal node
            cells <- which(colLinks(tse)$nodeNum %in% descs)
        }
        return(cells)
    }
    expect_equal(SummarizedExperiment::assay(out, "alias_1")[, 1],
                 rowMeans(SummarizedExperiment::assay(tse, "counts")[, intersect(
                     get_cells("alias_1"),
                     which(tse$sid == 1)), drop = FALSE]))
    expect_equal(SummarizedExperiment::assay(out, "alias_1")[, 2],
                 rowMeans(SummarizedExperiment::assay(tse, "counts")[, intersect(
                     get_cells("alias_1"),
                     which(tse$sid == 2)), drop = FALSE]))
    expect_equal(SummarizedExperiment::assay(out, "alias_2")[, 1],
                 rowMeans(SummarizedExperiment::assay(tse, "counts")[, intersect(
                     get_cells("alias_2"),
                     which(tse$sid == 1)), drop = FALSE]))
    expect_equal(SummarizedExperiment::assay(out, "alias_2")[, 2],
                 rowMeans(SummarizedExperiment::assay(tse, "counts")[, intersect(
                     get_cells("alias_2"),
                     which(tse$sid == 2)), drop = FALSE]))
    expect_equal(SummarizedExperiment::assay(out, "alias_3")[, 1],
                 rowMeans(SummarizedExperiment::assay(tse, "counts")[, intersect(
                     get_cells("alias_3"),
                     which(tse$sid == 1)), drop = FALSE]))
    expect_equal(SummarizedExperiment::assay(out, "alias_3")[, 2],
                 rowMeans(SummarizedExperiment::assay(tse, "counts")[, intersect(
                     get_cells("alias_3"),
                     which(tse$sid == 2)), drop = FALSE]))
    expect_equal(SummarizedExperiment::assay(out, "alias_4")[, 1],
                 rowMeans(SummarizedExperiment::assay(tse, "counts")[, intersect(
                     get_cells("alias_4"),
                     which(tse$sid == 1)), drop = FALSE]))
    expect_equal(SummarizedExperiment::assay(out, "alias_4")[, 2],
                 rowMeans(SummarizedExperiment::assay(tse, "counts")[, intersect(
                     get_cells("alias_4"),
                     which(tse$sid == 2)), drop = FALSE]))
    expect_equal(SummarizedExperiment::assay(out, "alias_5")[, 1],
                 rowMeans(SummarizedExperiment::assay(tse, "counts")[, intersect(
                     get_cells("alias_5"),
                     which(tse$sid == 1)), drop = FALSE]))
    expect_equal(SummarizedExperiment::assay(out, "alias_5")[, 2],
                 rowMeans(SummarizedExperiment::assay(tse, "counts")[, intersect(
                     get_cells("alias_5"),
                     which(tse$sid == 2)), drop = FALSE]))

})
