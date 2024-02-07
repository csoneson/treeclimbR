test_that("selNode works", {
    ## Generate example data
    library(TreeSummarizedExperiment)
    set.seed(1)
    data(tinyTree)
    toyTable <- matrix(rnbinom(40, size = 1, mu = 10), nrow = 10)
    colnames(toyTable) <- paste(rep(LETTERS[seq_len(2)], each = 2),
                                rep(seq_len(2), 2), sep = "_")
    rownames(toyTable) <- tinyTree$tip.label
    lse <- TreeSummarizedExperiment::TreeSummarizedExperiment(
        rowTree = tinyTree, assays = list(counts = toyTable))

    ## Check that function returns error with invalid input
    ## -------------------------------------------------------------------------
    expect_error(selNode(pr = NULL, obj = NULL, assay = 1, data = NULL,
                         tree = NULL, minTip = 0, maxTip = Inf, minPr = 0,
                         maxPr = 1, skip = NULL, all = FALSE),
                 "No valid input was found")
    expect_error(selNode(pr = TRUE, obj = NULL, assay = 1, data = NULL,
                         tree = NULL, minTip = 0, maxTip = Inf, minPr = 0,
                         maxPr = 1, skip = NULL, all = FALSE),
                 "'pr' must be of class 'numeric'")
    expect_error(selNode(pr = c(1, 2), obj = NULL, assay = 1, data = NULL,
                         tree = NULL, minTip = 0, maxTip = Inf, minPr = 0,
                         maxPr = 1, skip = NULL, all = FALSE),
                 "'pr' must be within [0,1]", fixed = TRUE)
    expect_error(selNode(pr = c(0.5, 0.5), obj = NULL, assay = 1, data = NULL,
                         tree = NULL, minTip = 0, maxTip = Inf, minPr = 0,
                         maxPr = 1, skip = NULL, all = FALSE),
                 "'tree' must not be NULL")
    expect_error(selNode(pr = NULL, obj = 3, assay = 1, data = NULL,
                         tree = NULL, minTip = 0, maxTip = Inf, minPr = 0,
                         maxPr = 1, skip = NULL, all = FALSE),
                 "No valid input was found")
    expect_error(selNode(pr = NULL, obj = lse, assay = 10, data = NULL,
                         tree = NULL, minTip = 0, maxTip = Inf, minPr = 0,
                         maxPr = 1, skip = NULL, all = FALSE),
                 "subscript is out of bounds")
    expect_error(selNode(pr = NULL, obj = lse, assay = "missing", data = NULL,
                         tree = NULL, minTip = 0, maxTip = Inf, minPr = 0,
                         maxPr = 1, skip = NULL, all = FALSE),
                 "assay %in% SummarizedExperiment::assayNames(obj)",
                 fixed = TRUE)
    expect_error(selNode(pr = NULL, obj = NULL, assay = 1, data = "x",
                         tree = NULL, minTip = 0, maxTip = Inf, minPr = 0,
                         maxPr = 1, skip = NULL, all = FALSE),
                 "'tree' must not be NULL")
    expect_error(selNode(pr = NULL, obj = NULL, assay = 1, data = "x",
                         tree = tinyTree, minTip = 0, maxTip = Inf, minPr = 0,
                         maxPr = 1, skip = NULL, all = FALSE),
                 'methods::is(obj, "matrix")', fixed = TRUE)
    expect_error(selNode(pr = NULL, obj = NULL, assay = 1, data = toyTable,
                         tree = 1, minTip = 0, maxTip = Inf, minPr = 0,
                         maxPr = 1, skip = NULL, all = FALSE),
                 "'tree' must be of class 'phylo'")
    expect_error(selNode(pr = c(0.5, 0.5), obj = NULL, assay = 1, data = NULL,
                         tree = tinyTree, minTip = 0, maxTip = Inf, minPr = 0,
                         maxPr = 1, skip = NULL, all = FALSE),
                 "!is.null(names(pars)) is not TRUE", fixed = TRUE)
    expect_error(selNode(pr = c(a = 0.5, b = 0.5), obj = NULL,
                         assay = 1, data = NULL,
                         tree = tinyTree, minTip = 0, maxTip = Inf, minPr = 0,
                         maxPr = 1, skip = NULL, all = FALSE),
                 "length(pars) == length(tree$tip.label) is not TRUE",
                 fixed = TRUE)
    expect_error(selNode(pr = structure(rep(0.1, 10),
                                        names = letters[seq_len(10)]),
                         obj = NULL, assay = 1, data = NULL,
                         tree = tinyTree, minTip = 0, maxTip = Inf, minPr = 0,
                         maxPr = 1, skip = NULL, all = FALSE),
                 "all(names(pars) %in% tree$tip.label) is not TRUE",
                 fixed = TRUE)
    expect_error(selNode(pr = structure(rep(0.1, 10),
                                        names = paste0("t", seq_len(10))),
                         obj = NULL, assay = 1, data = NULL,
                         tree = tinyTree, minTip = "x", maxTip = Inf, minPr = 0,
                         maxPr = 1, skip = NULL, all = FALSE),
                 "'minTip' must be of class 'numeric'")
    expect_error(selNode(pr = structure(rep(0.1, 10),
                                        names = paste0("t", seq_len(10))),
                         obj = NULL, assay = 1, data = NULL,
                         tree = tinyTree, minTip = c(1, 2),
                         maxTip = Inf, minPr = 0,
                         maxPr = 1, skip = NULL, all = FALSE),
                 "'minTip' must have length 1")
    expect_error(selNode(pr = structure(rep(0.1, 10),
                                        names = paste0("t", seq_len(10))),
                         obj = NULL, assay = 1, data = NULL,
                         tree = tinyTree, minTip = 0, maxTip = "x", minPr = 0,
                         maxPr = 1, skip = NULL, all = FALSE),
                 "'maxTip' must be of class 'numeric'")
    expect_error(selNode(pr = structure(rep(0.1, 10),
                                        names = paste0("t", seq_len(10))),
                         obj = NULL, assay = 1, data = NULL,
                         tree = tinyTree, minTip = 0, maxTip = c(1, 2),
                         minPr = 0,
                         maxPr = 1, skip = NULL, all = FALSE),
                 "'maxTip' must have length 1")
    expect_error(selNode(pr = structure(rep(0.1, 10),
                                        names = paste0("t", seq_len(10))),
                         obj = NULL, assay = 1, data = NULL,
                         tree = tinyTree, minTip = 0, maxTip = Inf, minPr = "x",
                         maxPr = 1, skip = NULL, all = FALSE),
                 "'minPr' must be of class 'numeric'")
    expect_error(selNode(pr = structure(rep(0.1, 10),
                                        names = paste0("t", seq_len(10))),
                         obj = NULL, assay = 1, data = NULL,
                         tree = tinyTree, minTip = 0, maxTip = Inf,
                         minPr = c(0.2, 0.5),
                         maxPr = 1, skip = NULL, all = FALSE),
                 "'minPr' must have length 1")
    expect_error(selNode(pr = structure(rep(0.1, 10),
                                        names = paste0("t", seq_len(10))),
                         obj = NULL, assay = 1, data = NULL,
                         tree = tinyTree, minTip = 0, maxTip = Inf, minPr = 2,
                         maxPr = 1, skip = NULL, all = FALSE),
                 "'minPr' must be within [0,1]", fixed = TRUE)
    expect_error(selNode(pr = structure(rep(0.1, 10),
                                        names = paste0("t", seq_len(10))),
                         obj = NULL, assay = 1, data = NULL,
                         tree = tinyTree, minTip = 0, maxTip = Inf, minPr = 0,
                         maxPr = "x", skip = NULL, all = FALSE),
                 "'maxPr' must be of class 'numeric'")
    expect_error(selNode(pr = structure(rep(0.1, 10),
                                        names = paste0("t", seq_len(10))),
                         obj = NULL, assay = 1, data = NULL,
                         tree = tinyTree, minTip = 0, maxTip = Inf, minPr = 0,
                         maxPr = c(0.2, 0.5), skip = NULL, all = FALSE),
                 "'maxPr' must have length 1")
    expect_error(selNode(pr = structure(rep(0.1, 10),
                                        names = paste0("t", seq_len(10))),
                         obj = NULL, assay = 1, data = NULL,
                         tree = tinyTree, minTip = 0, maxTip = Inf, minPr = 0,
                         maxPr = 2, skip = NULL, all = FALSE),
                 "'maxPr' must be within [0,1]", fixed = TRUE)
    expect_error(selNode(pr = structure(rep(0.1, 10),
                                        names = paste0("t", seq_len(10))),
                         obj = NULL, assay = 1, data = NULL,
                         tree = tinyTree, minTip = 0, maxTip = Inf, minPr = 0,
                         maxPr = 1, skip = NULL, all = 1),
                 "'all' must be of class 'logical'")
    expect_error(selNode(pr = structure(rep(0.1, 10),
                                        names = paste0("t", seq_len(10))),
                         obj = NULL, assay = 1, data = NULL,
                         tree = tinyTree, minTip = 0, maxTip = Inf, minPr = 0,
                         maxPr = 1, skip = NULL, all = c(TRUE, FALSE)),
                 "'all' must have length 1")


    ## Check that function works as expected with valid input - pr
    ## -------------------------------------------------------------------------
    dat <- parEstimate(obj = toyTable)
    pr <- dat$pi

    ## No restrictions
    out <- selNode(pr = pr, obj = NULL, assay = 1, data = NULL,
                   tree = tinyTree, minTip = 0, maxTip = Inf, minPr = 0,
                   maxPr = 1, skip = NULL, all = TRUE)
    expect_s3_class(out, "data.frame")
    expect_equal(nrow(out), 9)
    expect_named(out, c("nodeNum", "nodeLab", "proportion", "numTip"))
    expect_equal(out$nodeNum, seq(11, 19))
    expect_equal(out$proportion,
                 c(sum(pr), sum(pr[names(pr) != "t3"]),
                   sum(pr[c("t2", "t6", "t7")]), sum(pr[c("t6", "t7")]),
                   sum(pr[c("t1", "t4", "t5", "t8", "t9", "t10")]),
                   sum(pr[c("t1", "t4", "t8", "t9", "t10")]),
                   sum(pr[c("t4", "t8", "t9")]), sum(pr[c("t4", "t9")]),
                   sum(pr[c("t1", "t10")])))
    expect_equal(out$numTip,
                 c(10, 9, 3, 2, 6, 5, 3, 2, 2))

    ## Duplicated node labels
    tree2 <- tinyTree
    tree2$node.label[2] <- "Node_11"
    out <- selNode(pr = pr, obj = NULL, assay = 1, data = NULL,
                   tree = tree2, minTip = 0, maxTip = Inf, minPr = 0,
                   maxPr = 1, skip = NULL, all = TRUE)
    expect_s3_class(out, "data.frame")
    expect_equal(nrow(out), 9)
    expect_named(out, c("nodeNum", "nodeLab", "nodeLab_alias",
                        "proportion", "numTip"))
    expect_equal(out$nodeNum, seq(11, 19))
    expect_equal(out$proportion,
                 c(sum(pr), sum(pr[names(pr) != "t3"]),
                   sum(pr[c("t2", "t6", "t7")]), sum(pr[c("t6", "t7")]),
                   sum(pr[c("t1", "t4", "t5", "t8", "t9", "t10")]),
                   sum(pr[c("t1", "t4", "t8", "t9", "t10")]),
                   sum(pr[c("t4", "t8", "t9")]), sum(pr[c("t4", "t9")]),
                   sum(pr[c("t1", "t10")])))
    expect_equal(out$numTip,
                 c(10, 9, 3, 2, 6, 5, 3, 2, 2))

    ## Restrict proportion
    expect_message({
        out <- selNode(pr = pr, obj = lse, assay = 1, data = NULL,
                       tree = tinyTree, minTip = 0, maxTip = Inf, minPr = 0.2,
                       maxPr = 0.5, skip = NULL, all = TRUE)
    }, "Ignoring obj when input is a proportions vector")
    expect_s3_class(out, "data.frame")
    expect_equal(nrow(out), 2)
    expect_named(out, c("nodeNum", "nodeLab", "proportion", "numTip"))
    expect_equal(out$nodeNum, c(13, 17))
    expect_equal(out$proportion,
                 c(sum(pr[c("t2", "t6", "t7")]),
                   sum(pr[c("t4", "t8", "t9")])))
    expect_equal(out$numTip,
                 c(3, 3))

    ## Restrict proportion too much
    expect_error(selNode(pr = pr, obj = NULL, assay = 1, data = NULL,
                         tree = tinyTree, minTip = 0, maxTip = Inf, minPr = 0,
                         maxPr = 0.1, skip = NULL, all = TRUE),
                 "All nodes have proportions exceeding maxPr")

    expect_error(selNode(pr = pr, obj = NULL, assay = 1, data = NULL,
                         tree = tinyTree, minTip = 0, maxTip = Inf,
                         minPr = 0.6, maxPr = 0.4, skip = NULL, all = TRUE),
                 "No nodes fulfill the requirements")

    expect_error(selNode(pr = pr, obj = NULL, assay = 1, data = NULL,
                         tree = tinyTree, minTip = 0, maxTip = Inf,
                         minPr = 0, maxPr = 1, skip = seq_len(10), all = TRUE),
                 "No nodes fulfill the requirements")

    ## Return only single node
    expect_message({
        out <- selNode(pr = pr, obj = NULL, assay = 1, data = toyTable,
                       tree = tinyTree, minTip = 0, maxTip = Inf, minPr = 0.2,
                       maxPr = 0.5, skip = NULL, all = FALSE)
    }, "Ignoring data when input is a proportions vector")
    expect_s3_class(out, "data.frame")
    expect_equal(nrow(out), 1)
    expect_named(out, c("nodeNum", "nodeLab", "proportion", "numTip"))
    expect_equal(out$nodeNum, c(13))
    expect_equal(out$proportion,
                 c(sum(pr[c("t2", "t6", "t7")])))
    expect_equal(out$numTip,
                 c(3))

    ## Skip some nodes
    out <- selNode(pr = pr, obj = NULL, assay = 1, data = NULL,
                   tree = tinyTree, minTip = 0, maxTip = Inf, minPr = 0,
                   maxPr = 1, skip = c(16), all = TRUE)
    expect_s3_class(out, "data.frame")
    expect_equal(nrow(out), 2)
    expect_named(out, c("nodeNum", "nodeLab", "proportion", "numTip"))
    expect_equal(out$nodeNum, seq(13, 14))
    expect_equal(out$proportion,
                 c(sum(pr[c("t2", "t6", "t7")]), sum(pr[c("t6", "t7")])))
    expect_equal(out$numTip,
                 c(3, 2))

    ## Skip, but with node label
    out <- selNode(pr = pr, obj = NULL, assay = 1, data = NULL,
                   tree = tinyTree, minTip = 0, maxTip = Inf, minPr = 0,
                   maxPr = 1, skip = c("Node_16"), all = TRUE)
    expect_s3_class(out, "data.frame")
    expect_equal(nrow(out), 2)
    expect_named(out, c("nodeNum", "nodeLab", "proportion", "numTip"))
    expect_equal(out$nodeNum, seq(13, 14))
    expect_equal(out$proportion,
                 c(sum(pr[c("t2", "t6", "t7")]), sum(pr[c("t6", "t7")])))
    expect_equal(out$numTip,
                 c(3, 2))

    ## Skip tips
    out <- selNode(pr = pr, obj = NULL, assay = 1, data = NULL,
                   tree = tinyTree, minTip = 0, maxTip = Inf, minPr = 0,
                   maxPr = 1, skip = c("t10", "t6"), all = TRUE)
    expect_s3_class(out, "data.frame")
    expect_equal(nrow(out), 2)
    expect_named(out, c("nodeNum", "nodeLab", "proportion", "numTip"))
    expect_equal(out$nodeNum, seq(17, 18))
    expect_equal(out$proportion,
                 c(sum(pr[c("t4", "t8", "t9")]), sum(pr[c("t4", "t9")])))
    expect_equal(out$numTip,
                 c(3, 2))

    ## Check that function works as expected with valid input - TSE
    ## -------------------------------------------------------------------------
    dat <- parEstimate(obj = toyTable)
    pr <- dat$pi

    ## No restrictions
    expect_message({
        out <- selNode(pr = NULL, obj = lse, assay = 1, data = NULL,
                       tree = tinyTree, minTip = 0, maxTip = Inf, minPr = 0,
                       maxPr = 1, skip = NULL, all = TRUE)
    }, "Ignoring tree when input is a TSE")
    expect_s3_class(out, "data.frame")
    expect_equal(nrow(out), 9)
    expect_named(out, c("nodeNum", "nodeLab", "proportion", "numTip"))
    expect_equal(out$nodeNum, seq(11, 19))
    expect_equal(out$proportion,
                 c(sum(pr), sum(pr[names(pr) != "t3"]),
                   sum(pr[c("t2", "t6", "t7")]), sum(pr[c("t6", "t7")]),
                   sum(pr[c("t1", "t4", "t5", "t8", "t9", "t10")]),
                   sum(pr[c("t1", "t4", "t8", "t9", "t10")]),
                   sum(pr[c("t4", "t8", "t9")]), sum(pr[c("t4", "t9")]),
                   sum(pr[c("t1", "t10")])))
    expect_equal(out$numTip,
                 c(10, 9, 3, 2, 6, 5, 3, 2, 2))

    ## Restrict proportion
    expect_message({
        out <- selNode(pr = NULL, obj = lse, assay = "counts", data = dat,
                       tree = NULL, minTip = 0, maxTip = Inf, minPr = 0.2,
                       maxPr = 0.5, skip = NULL, all = TRUE)
    }, "Ignoring data when input is a TSE")
    expect_s3_class(out, "data.frame")
    expect_equal(nrow(out), 2)
    expect_named(out, c("nodeNum", "nodeLab", "proportion", "numTip"))
    expect_equal(out$nodeNum, c(13, 17))
    expect_equal(out$proportion,
                 c(sum(pr[c("t2", "t6", "t7")]),
                   sum(pr[c("t4", "t8", "t9")])))
    expect_equal(out$numTip,
                 c(3, 3))

    ## Check that function works as expected with valid input - data + tree
    ## -------------------------------------------------------------------------
    dat <- parEstimate(obj = toyTable)
    pr <- dat$pi

    ## No restrictions
    out <- selNode(pr = NULL, obj = NULL, assay = 1, data = dat,
                       tree = tinyTree, minTip = 0, maxTip = Inf, minPr = 0,
                       maxPr = 1, skip = NULL, all = TRUE)
    expect_s3_class(out, "data.frame")
    expect_equal(nrow(out), 9)
    expect_named(out, c("nodeNum", "nodeLab", "proportion", "numTip"))
    expect_equal(out$nodeNum, seq(11, 19))
    expect_equal(out$proportion,
                 c(sum(pr), sum(pr[names(pr) != "t3"]),
                   sum(pr[c("t2", "t6", "t7")]), sum(pr[c("t6", "t7")]),
                   sum(pr[c("t1", "t4", "t5", "t8", "t9", "t10")]),
                   sum(pr[c("t1", "t4", "t8", "t9", "t10")]),
                   sum(pr[c("t4", "t8", "t9")]), sum(pr[c("t4", "t9")]),
                   sum(pr[c("t1", "t10")])))
    expect_equal(out$numTip,
                 c(10, 9, 3, 2, 6, 5, 3, 2, 2))

    ## Restrict proportion
    out <- selNode(pr = NULL, obj = NULL, assay = 1, data = dat,
                   tree = tinyTree, minTip = 0, maxTip = Inf, minPr = 0.2,
                   maxPr = 0.5, skip = NULL, all = TRUE)
    expect_s3_class(out, "data.frame")
    expect_equal(nrow(out), 2)
    expect_named(out, c("nodeNum", "nodeLab", "proportion", "numTip"))
    expect_equal(out$nodeNum, c(13, 17))
    expect_equal(out$proportion,
                 c(sum(pr[c("t2", "t6", "t7")]),
                   sum(pr[c("t4", "t8", "t9")])))
    expect_equal(out$numTip,
                 c(3, 3))
})
