test_that("medianByClusterMarker works", {
    ## Generate some data
    library(SummarizedExperiment)
    set.seed(1)
    count <- matrix(rpois(n = 1000, lambda = 10), nrow = 100)
    colnames(count) <- paste0("mk", 1:10)
    rowD <- data.frame("cluster" = sample(seq_len(6), 100, replace = TRUE))
    colD <- data.frame(type_marker = rep(c(FALSE, TRUE), each = 5))

    ## Markers in columns ('transposed' SE)
    d_tr <- SummarizedExperiment(assays = list(counts = count),
                                 rowData = rowD, colData = colD)
    ## Markers in rows ('regular' SE)
    d_reg <- SummarizedExperiment(assays = list(counts = t(count)),
                                  rowData = colD, colData = rowD)

    ## Check that function returns error with invalid input
    ## -------------------------------------------------------------------------
    expect_error(medianByClusterMarker(
        SE = 1, assay = 1, marker_in_column = FALSE, column_cluster = "cluster",
        use_marker = SummarizedExperiment::rowData(d_reg)$type_marker
    ), "'SE' must be of class 'SummarizedExperiment'")
    expect_error(medianByClusterMarker(
        SE = d_reg, assay = c(1, 2), marker_in_column = 1,
        column_cluster = "cluster",
        use_marker = SummarizedExperiment::rowData(d_reg)$type_marker
    ), "length(assay) == 1", fixed = TRUE)
    expect_error(medianByClusterMarker(
        SE = d_reg, assay = TRUE, marker_in_column = 1,
        column_cluster = "cluster",
        use_marker = SummarizedExperiment::rowData(d_reg)$type_marker
    ), "length(assay) == 1", fixed = TRUE)
    expect_error(medianByClusterMarker(
        SE = d_reg, assay = 1, marker_in_column = 1, column_cluster = "cluster",
        use_marker = SummarizedExperiment::rowData(d_reg)$type_marker
    ), "'marker_in_column' must be of class 'logical'")
    expect_error(medianByClusterMarker(
        SE = d_reg, assay = 1, marker_in_column = c(TRUE, FALSE),
        column_cluster = "cluster",
        use_marker = SummarizedExperiment::rowData(d_reg)$type_marker
    ), "'marker_in_column' must have length 1")
    expect_error(medianByClusterMarker(
        SE = d_reg, assay = 1, marker_in_column = FALSE, column_cluster = 1,
        use_marker = SummarizedExperiment::rowData(d_reg)$type_marker
    ), "'column_cluster' must be of class 'character'")
    expect_error(medianByClusterMarker(
        SE = d_reg, assay = 1, marker_in_column = FALSE,
        column_cluster = c("cluster", "cluster"),
        use_marker = SummarizedExperiment::rowData(d_reg)$type_marker
    ), "'column_cluster' must have length 1")

    ## Markers in rows - check that function works with valid input
    ## -------------------------------------------------------------------------
    ## ... 'type' markers
    out <- medianByClusterMarker(
        SE = d_reg, assay = "counts", marker_in_column = FALSE,
        column_cluster = "cluster",
        use_marker = SummarizedExperiment::rowData(d_reg)$type_marker
    )
    expect_s4_class(out, "SummarizedExperiment")
    expect_equal(nrow(out), 5)
    expect_equal(ncol(out), 6)
    expect_equal(rownames(out), paste0("mk", seq(6, 10)))
    expect_equal(colnames(out), as.character(seq_len(6)))
    expect_equal(SummarizedExperiment::assayNames(out), "counts")
    expect_equal(out$cluster, colnames(out))
    for (m in rownames(out)) {
        for (cl in colnames(out)) {
            expect_equal(
                SummarizedExperiment::assay(out, "counts")[m, cl],
                median(SummarizedExperiment::assay(
                    d_reg, "counts")[m, d_reg$cluster == cl]))
        }
    }

    ## ... all markers
    out <- medianByClusterMarker(
        SE = d_reg, assay = "counts", marker_in_column = FALSE,
        column_cluster = "cluster", use_marker = NULL
    )
    expect_s4_class(out, "SummarizedExperiment")
    expect_equal(nrow(out), 10)
    expect_equal(ncol(out), 6)
    expect_equal(rownames(out), paste0("mk", seq(1, 10)))
    expect_equal(colnames(out), as.character(seq_len(6)))
    expect_equal(SummarizedExperiment::assayNames(out), "counts")
    expect_equal(out$cluster, colnames(out))
    for (m in rownames(out)) {
        for (cl in colnames(out)) {
            expect_equal(
                SummarizedExperiment::assay(out, "counts")[m, cl],
                median(SummarizedExperiment::assay(
                    d_reg, "counts")[m, d_reg$cluster == cl]))
        }
    }

    ## ... markers specified by numeric indices
    out <- medianByClusterMarker(
        SE = d_reg, assay = "counts", marker_in_column = FALSE,
        column_cluster = "cluster",
        use_marker = which(SummarizedExperiment::rowData(d_reg)$type_marker)
    )
    expect_s4_class(out, "SummarizedExperiment")
    expect_equal(nrow(out), 5)
    expect_equal(ncol(out), 6)
    expect_equal(rownames(out), paste0("mk", seq(6, 10)))
    expect_equal(colnames(out), as.character(seq_len(6)))
    expect_equal(SummarizedExperiment::assayNames(out), "counts")
    expect_equal(out$cluster, colnames(out))
    for (m in rownames(out)) {
        for (cl in colnames(out)) {
            expect_equal(
                SummarizedExperiment::assay(out, "counts")[m, cl],
                median(SummarizedExperiment::assay(
                    d_reg, "counts")[m, d_reg$cluster == cl]))
        }
    }

    ## ... markers specified by name
    out <- medianByClusterMarker(
        SE = d_reg, assay = "counts", marker_in_column = FALSE,
        column_cluster = "cluster",
        use_marker = paste0("mk", seq(10, 6))
    )
    expect_s4_class(out, "SummarizedExperiment")
    expect_equal(nrow(out), 5)
    expect_equal(ncol(out), 6)
    expect_equal(rownames(out), paste0("mk", seq(10, 6)))
    expect_equal(colnames(out), as.character(seq_len(6)))
    expect_equal(SummarizedExperiment::assayNames(out), "counts")
    expect_equal(out$cluster, colnames(out))
    for (m in rownames(out)) {
        for (cl in colnames(out)) {
            expect_equal(
                SummarizedExperiment::assay(out, "counts")[m, cl],
                median(SummarizedExperiment::assay(
                    d_reg, "counts")[m, d_reg$cluster == cl]))
        }
    }

    ## ... specify assay ID rather than name
    out <- medianByClusterMarker(
        SE = d_reg, assay = 1, marker_in_column = FALSE,
        column_cluster = "cluster",
        use_marker = which(SummarizedExperiment::rowData(d_reg)$type_marker)
    )
    expect_s4_class(out, "SummarizedExperiment")
    expect_equal(nrow(out), 5)
    expect_equal(ncol(out), 6)
    expect_equal(rownames(out), paste0("mk", seq(6, 10)))
    expect_equal(colnames(out), as.character(seq_len(6)))
    expect_equal(SummarizedExperiment::assayNames(out), "counts")
    expect_equal(out$cluster, colnames(out))
    for (m in rownames(out)) {
        for (cl in colnames(out)) {
            expect_equal(
                SummarizedExperiment::assay(out, "counts")[m, cl],
                median(SummarizedExperiment::assay(
                    d_reg, "counts")[m, d_reg$cluster == cl]))
        }
    }

    ## ... specify assay ID rather than name, no assay name in input
    tmp <- d_reg
    SummarizedExperiment::assays(tmp) <- list(SummarizedExperiment::assay(tmp))
    out <- medianByClusterMarker(
        SE = tmp, assay = 1, marker_in_column = FALSE,
        column_cluster = "cluster",
        use_marker = which(SummarizedExperiment::rowData(tmp)$type_marker)
    )
    expect_s4_class(out, "SummarizedExperiment")
    expect_equal(nrow(out), 5)
    expect_equal(ncol(out), 6)
    expect_equal(rownames(out), paste0("mk", seq(6, 10)))
    expect_equal(colnames(out), as.character(seq_len(6)))
    expect_null(SummarizedExperiment::assayNames(out))
    expect_equal(out$cluster, colnames(out))
    for (m in rownames(out)) {
        for (cl in colnames(out)) {
            expect_equal(
                SummarizedExperiment::assay(out)[m, cl],
                median(SummarizedExperiment::assay(
                    d_reg)[m, d_reg$cluster == cl]))
        }
    }

    ## Markers in columns - check that function works with valid input
    ## -------------------------------------------------------------------------
    ## ... 'type' markers
    out <- medianByClusterMarker(
        SE = d_tr, assay = "counts", marker_in_column = TRUE,
        column_cluster = "cluster",
        use_marker = SummarizedExperiment::colData(d_tr)$type_marker
    )
    expect_s4_class(out, "SummarizedExperiment")
    expect_equal(nrow(out), 6)
    expect_equal(ncol(out), 5)
    expect_equal(colnames(out), paste0("mk", seq(6, 10)))
    expect_equal(rownames(out), as.character(seq_len(6)))
    expect_equal(SummarizedExperiment::assayNames(out), "counts")
    expect_equal(SummarizedExperiment::rowData(out)$cluster, rownames(out))
    for (m in colnames(out)) {
        for (cl in rownames(out)) {
            expect_equal(
                SummarizedExperiment::assay(out, "counts")[cl, m],
                median(SummarizedExperiment::assay(
                    d_tr, "counts")[
                        SummarizedExperiment::rowData(d_tr)$cluster == cl, m]))
        }
    }

    ## ... all markers
    out <- medianByClusterMarker(
        SE = d_tr, assay = "counts", marker_in_column = TRUE,
        column_cluster = "cluster", use_marker = NULL
    )
    expect_s4_class(out, "SummarizedExperiment")
    expect_equal(nrow(out), 6)
    expect_equal(ncol(out), 10)
    expect_equal(colnames(out), paste0("mk", seq(1, 10)))
    expect_equal(rownames(out), as.character(seq_len(6)))
    expect_equal(SummarizedExperiment::assayNames(out), "counts")
    expect_equal(SummarizedExperiment::rowData(out)$cluster, rownames(out))
    for (m in colnames(out)) {
        for (cl in rownames(out)) {
            expect_equal(
                SummarizedExperiment::assay(out, "counts")[cl, m],
                median(SummarizedExperiment::assay(
                    d_tr, "counts")[
                        SummarizedExperiment::rowData(d_tr)$cluster == cl, m]))
        }
    }

    ## ... markers specified by numeric indices
    out <- medianByClusterMarker(
        SE = d_tr, assay = "counts", marker_in_column = TRUE,
        column_cluster = "cluster",
        use_marker = which(SummarizedExperiment::colData(d_tr)$type_marker)
    )
    expect_s4_class(out, "SummarizedExperiment")
    expect_equal(nrow(out), 6)
    expect_equal(ncol(out), 5)
    expect_equal(colnames(out), paste0("mk", seq(6, 10)))
    expect_equal(rownames(out), as.character(seq_len(6)))
    expect_equal(SummarizedExperiment::assayNames(out), "counts")
    expect_equal(SummarizedExperiment::rowData(out)$cluster, rownames(out))
    for (m in colnames(out)) {
        for (cl in rownames(out)) {
            expect_equal(
                SummarizedExperiment::assay(out, "counts")[cl, m],
                median(SummarizedExperiment::assay(
                    d_tr, "counts")[
                        SummarizedExperiment::rowData(d_tr)$cluster == cl, m]))
        }
    }

    ## ... markers specified by name
    out <- medianByClusterMarker(
        SE = d_tr, assay = "counts", marker_in_column = TRUE,
        column_cluster = "cluster",
        use_marker = paste0("mk", seq(10, 6))
    )
    expect_s4_class(out, "SummarizedExperiment")
    expect_equal(nrow(out), 6)
    expect_equal(ncol(out), 5)
    expect_equal(colnames(out), paste0("mk", seq(10, 6)))
    expect_equal(rownames(out), as.character(seq_len(6)))
    expect_equal(SummarizedExperiment::assayNames(out), "counts")
    expect_equal(SummarizedExperiment::rowData(out)$cluster, rownames(out))
    for (m in colnames(out)) {
        for (cl in rownames(out)) {
            expect_equal(
                SummarizedExperiment::assay(out, "counts")[cl, m],
                median(SummarizedExperiment::assay(
                    d_tr, "counts")[
                        SummarizedExperiment::rowData(d_tr)$cluster == cl, m]))
        }
    }

    ## ... specify assay ID rather than name
    out <- medianByClusterMarker(
        SE = d_tr, assay = 1, marker_in_column = TRUE,
        column_cluster = "cluster",
        use_marker = which(SummarizedExperiment::colData(d_tr)$type_marker)
    )
    expect_s4_class(out, "SummarizedExperiment")
    expect_equal(nrow(out), 6)
    expect_equal(ncol(out), 5)
    expect_equal(colnames(out), paste0("mk", seq(6, 10)))
    expect_equal(rownames(out), as.character(seq_len(6)))
    expect_equal(SummarizedExperiment::assayNames(out), "counts")
    expect_equal(SummarizedExperiment::rowData(out)$cluster, rownames(out))
    for (m in colnames(out)) {
        for (cl in rownames(out)) {
            expect_equal(
                SummarizedExperiment::assay(out, "counts")[cl, m],
                median(SummarizedExperiment::assay(
                    d_tr, "counts")[
                        SummarizedExperiment::rowData(d_tr)$cluster == cl, m]))
        }
    }

    ## ... specify assay ID rather than name, no assay name in input
    tmp <- d_tr
    SummarizedExperiment::assays(tmp) <- list(SummarizedExperiment::assay(tmp))
    out <- medianByClusterMarker(
        SE = tmp, assay = 1, marker_in_column = TRUE,
        column_cluster = "cluster",
        use_marker = which(SummarizedExperiment::colData(tmp)$type_marker)
    )
    expect_s4_class(out, "SummarizedExperiment")
    expect_equal(nrow(out), 6)
    expect_equal(ncol(out), 5)
    expect_equal(colnames(out), paste0("mk", seq(6, 10)))
    expect_equal(rownames(out), as.character(seq_len(6)))
    expect_null(SummarizedExperiment::assayNames(out))
    expect_equal(SummarizedExperiment::rowData(out)$cluster, rownames(out))
    for (m in colnames(out)) {
        for (cl in rownames(out)) {
            expect_equal(
                SummarizedExperiment::assay(out)[cl, m],
                median(SummarizedExperiment::assay(
                    tmp)[SummarizedExperiment::rowData(d_tr)$cluster == cl, m]))
        }
    }


})
