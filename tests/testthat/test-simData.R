test_that("simData works", {
    ## Generate data to use as the starting point
    set.seed(1L)
    y <- matrix(rnbinom(120, size = 1, mu = 10), nrow = 10)
    colnames(y) <- paste("S", seq_len(12), sep = "")
    rownames(y) <- tinyTree$tip.label
    toy_lse <- TreeSummarizedExperiment(rowTree = tinyTree,
                                        assays = list(counts = y))

    ## Check that function returns error if provided with invalid input
    ## -------------------------------------------------------------------------
    .args <- list(tree = NULL, data = NULL, obj = NULL, assay = NULL,
                  scenario = "BS", from.A = NULL, from.B = NULL,
                  minTip.A = 0, maxTip.A = Inf, minTip.B = 0, maxTip.B = Inf,
                  minPr.A = 0, maxPr.A = 1, ratio = 4, adjB = NULL,
                  pct = 0.6, nSam = c(50, 50), mu = 10000, size = NULL,
                  n = 1, FUN = sum, message = FALSE)
    expect_error(do.call(simData, .args),
                 "(!is.null(tree) && !is.null(data)) || !is.null(obj)",
                 fixed = TRUE)

    args <- .args
    args$tree <- tinyTree
    args$data <- y[seq_len(5), ]
    expect_error(do.call(simData, args),
                 "The rownames of data do not match the leaf labels.")
    args$data <- "x"
    expect_error(do.call(simData, args),
                 "data should be a matrix")

    args <- .args
    args$obj <- 1
    expect_error(do.call(simData, args),
                 "'obj' must be of class 'TreeSummarizedExperiment'")

    args <- .args
    args$obj <- toy_lse
    args$assay <- "missing"
    expect_error(do.call(simData, args),
                 "assay %in% SummarizedExperiment::assayNames(obj) is not TRUE",
                 fixed = TRUE)

    args0 <- .args
    args0$obj <- toy_lse
    args0$assay <- "counts"

    args <- args0
    args$scenario <- 1
    expect_error(do.call(simData, args),
                 "'scenario' must be of class 'character'")
    args$scenario <- c("BS", "US")
    expect_error(do.call(simData, args),
                 "'scenario' must have length 1")
    args$scenario <- "missing"
    expect_error(do.call(simData, args),
                 "'scenario' must be one of")

    args <- args0
    args$from.A <- "missing"
    expect_error(do.call(simData, args),
                 "The provided from.A is not a node in the tree")
    args <- args0
    args$from.A <- 100
    expect_error(do.call(simData, args),
                 "The provided from.A is not a node in the tree")
    args$from.A <- TRUE
    expect_error(do.call(simData, args),
                 "from.A must be a character or numeric value")

    args <- args0
    args$from.A <- 14
    args$from.B <- "missing"
    expect_error(do.call(simData, args),
                 "The provided from.B is not a node in the tree")
    args$from.B <- 100
    expect_error(do.call(simData, args),
                 "The provided from.B is not a node in the tree")
    args$from.B <- TRUE
    expect_error(do.call(simData, args),
                 "from.B must be a character or numeric value")

    args <- args0
    args$minTip.A <- "x"
    expect_error(do.call(simData, args),
                 "'minTip.A' must be of class 'numeric'")
    args$minTip.A <- c(1, 2)
    expect_error(do.call(simData, args),
                 "'minTip.A' must have length 1")

    args <- args0
    args$minTip.B <- "x"
    expect_error(do.call(simData, args),
                 "'minTip.B' must be of class 'numeric'")
    args$minTip.B <- c(1, 2)
    expect_error(do.call(simData, args),
                 "'minTip.B' must have length 1")

    args <- args0
    args$minPr.A <- "x"
    expect_error(do.call(simData, args),
                 "'minPr.A' must be of class 'numeric'")
    args$minPr.A <- c(0.1, 0.2)
    expect_error(do.call(simData, args),
                 "'minPr.A' must have length 1")
    args$minPr.A <- 3
    expect_error(do.call(simData, args),
                 "'minPr.A' must be within [0,1]", fixed = TRUE)

    args <- args0
    args$maxPr.A <- "x"
    expect_error(do.call(simData, args),
                 "'maxPr.A' must be of class 'numeric'")
    args$maxPr.A <- c(0.1, 0.2)
    expect_error(do.call(simData, args),
                 "'maxPr.A' must have length 1")
    args$maxPr.A <- 3
    expect_error(do.call(simData, args),
                 "'maxPr.A' must be within [0,1]", fixed = TRUE)
    args$maxPr.A <- 0
    expect_error(do.call(simData, args),
                 "maxPr.A is lower than the minimum", fixed = TRUE)

    args <- args0
    args$ratio <- "x"
    expect_error(do.call(simData, args),
                 "'ratio' must be of class 'numeric'")
    args$ratio <- c(1, 2)
    expect_error(do.call(simData, args),
                 "'ratio' must have length 1")
    args$ratio <- 100
    expect_error(do.call(simData, args),
                 "Could not find two branches which fulfill the requirement")
    args$ratio <- 0.0001
    expect_error(do.call(simData, args),
                 "Could not find two branches which fulfill the requirement")
    args$ratio <- 4
    args$minPr.A <- 0.8
    expect_error(do.call(simData, args),
                 "minPr.A*ratio is above the maximum value of", fixed = TRUE)

    args <- args0
    args$adjB <- "x"
    expect_error(do.call(simData, args),
                 "'adjB' must be of class 'numeric'")
    args$adjB <- c(0.1, 0.2)
    expect_error(do.call(simData, args),
                 "'adjB' must have length 1")
    args$adjB <- 3
    expect_error(do.call(simData, args),
                 "'adjB' must be within [0,1]", fixed = TRUE)

    args <- args0
    args$pct <- "x"
    expect_error(do.call(simData, args),
                 "'pct' must be of class 'numeric'")
    args$pct <- c(0.1, 0.2)
    expect_error(do.call(simData, args),
                 "'pct' must have length 1")
    args$pct <- 3
    expect_error(do.call(simData, args),
                 "'pct' must be within [0,1]", fixed = TRUE)

    args <- args0
    args$nSam <- "x"
    expect_error(do.call(simData, args),
                 "'nSam' must be of class 'numeric'")

    args <- args0
    args$mu <- "x"
    expect_error(do.call(simData, args),
                 "'mu' must be of class 'numeric'")

    args <- args0
    args$size <- "x"
    expect_error(do.call(simData, args),
                 "'size' must be of class 'numeric'")

    args <- args0
    args$n <- "x"
    expect_error(do.call(simData, args),
                 "'n' must be of class 'numeric'")
    args$n <- c(1, 2)
    expect_error(do.call(simData, args),
                 "'n' must have length 1")

    args <- args0
    args$FUN <- "x"
    expect_error(do.call(simData, args),
                 "'FUN' must be of class 'function'")

    args <- args0
    args$message <- "x"
    expect_error(do.call(simData, args),
                 "'message' must be of class 'logical'")
    args$message <- c(TRUE, FALSE)
    expect_error(do.call(simData, args),
                 "'message' must have length 1")

    ## Check that function works as expected with valid input
    ## -------------------------------------------------------------------------
    ## Let function choose branches
    set.seed(1)
    out <- simData(tree = NULL, data = NULL, obj = toy_lse, assay = "counts",
                   scenario = "BS", from.A = NULL, from.B = NULL,
                   minTip.A = 0, maxTip.A = Inf, minTip.B = 0, maxTip.B = Inf,
                   minPr.A = 0, maxPr.A = 1, ratio = 4, adjB = NULL,
                   pct = 0.6, nSam = c(50, 50), mu = 10000, size = NULL,
                   n = 1, FUN = sum, message = FALSE)
    expect_s4_class(out, "TreeSummarizedExperiment")
    expect_equal(nrow(out), nrow(toy_lse))
    expect_equal(ncol(out), 100)
    expect_equal(SummarizedExperiment::assayNames(out), "counts")
    expect_equal(TreeSummarizedExperiment::rowTree(out),
                 TreeSummarizedExperiment::rowTree(toy_lse))
    expect_equal(rownames(out), rownames(toy_lse))
    expect_equal(out$group, rep(c("C1", "C2"), each = 50))
    expect_type(S4Vectors::metadata(out), "list")
    expect_equal(prod(unique(S4Vectors::metadata(out)$FC)), 1)
    expect_s3_class(S4Vectors::metadata(out)$branch, "data.frame")
    expect_equal(S4Vectors::metadata(out)$branch$A, 14)
    expect_equal(S4Vectors::metadata(out)$branch$B, 15)

    ## Let function choose branches - with restrictions
    set.seed(1)
    out <- simData(tree = NULL, data = NULL, obj = toy_lse, assay = "counts",
                   scenario = "BS", from.A = NULL, from.B = NULL,
                   minTip.A = 1, maxTip.A = 5, minTip.B = 2, maxTip.B = 6,
                   minPr.A = 0.1, maxPr.A = 0.9, ratio = 2, adjB = NULL,
                   pct = 0.6, nSam = c(50, 50), mu = 10000, size = NULL,
                   n = 1, FUN = sum, message = FALSE)
    expect_s4_class(out, "TreeSummarizedExperiment")
    expect_equal(nrow(out), nrow(toy_lse))
    expect_equal(ncol(out), 100)
    expect_equal(SummarizedExperiment::assayNames(out), "counts")
    expect_equal(TreeSummarizedExperiment::rowTree(out),
                 TreeSummarizedExperiment::rowTree(toy_lse))
    expect_equal(rownames(out), rownames(toy_lse))
    expect_equal(out$group, rep(c("C1", "C2"), each = 50))
    expect_type(S4Vectors::metadata(out), "list")
    expect_equal(prod(unique(S4Vectors::metadata(out)$FC)), 1)
    expect_s3_class(S4Vectors::metadata(out)$branch, "data.frame")
    expect_equal(S4Vectors::metadata(out)$branch$A, 13)
    expect_equal(S4Vectors::metadata(out)$branch$B, 15)

    ## Fix branches - numbers
    set.seed(1)
    out <- simData(tree = NULL, data = NULL, obj = toy_lse, assay = "counts",
                   scenario = "BS", from.A = 18, from.B = 19,
                   minTip.A = 0, maxTip.A = Inf, minTip.B = 0, maxTip.B = Inf,
                   minPr.A = 0, maxPr.A = 1, ratio = 4, adjB = NULL,
                   pct = 0.6, nSam = c(50, 50), mu = 10000, size = NULL,
                   n = 1, FUN = sum, message = FALSE)
    expect_s4_class(out, "TreeSummarizedExperiment")
    expect_equal(nrow(out), nrow(toy_lse))
    expect_equal(ncol(out), 100)
    expect_equal(SummarizedExperiment::assayNames(out), "counts")
    expect_equal(TreeSummarizedExperiment::rowTree(out),
                 TreeSummarizedExperiment::rowTree(toy_lse))
    expect_equal(rownames(out), rownames(toy_lse))
    expect_equal(out$group, rep(c("C1", "C2"), each = 50))
    expect_type(S4Vectors::metadata(out), "list")
    expect_equal(prod(unique(S4Vectors::metadata(out)$FC)), 1)
    expect_s3_class(S4Vectors::metadata(out)$branch, "data.frame")
    expect_equal(S4Vectors::metadata(out)$branch$A, 18)
    expect_equal(S4Vectors::metadata(out)$branch$B, 19)
    expect_true(round(S4Vectors::metadata(out)$branch$ratio, 5) %in%
                    round(S4Vectors::metadata(out)$FC, 5))

    ## Fix branch A - numbers
    set.seed(1)
    out <- simData(tree = NULL, data = NULL, obj = toy_lse, assay = "counts",
                   scenario = "BS", from.A = 18, from.B = NULL,
                   minTip.A = 0, maxTip.A = Inf, minTip.B = 0, maxTip.B = Inf,
                   minPr.A = 0, maxPr.A = 1, ratio = 4, adjB = NULL,
                   pct = 0.6, nSam = c(50, 50), mu = 10000, size = NULL,
                   n = 1, FUN = sum, message = FALSE)
    expect_s4_class(out, "TreeSummarizedExperiment")
    expect_equal(nrow(out), nrow(toy_lse))
    expect_equal(ncol(out), 100)
    expect_equal(SummarizedExperiment::assayNames(out), "counts")
    expect_equal(TreeSummarizedExperiment::rowTree(out),
                 TreeSummarizedExperiment::rowTree(toy_lse))
    expect_equal(rownames(out), rownames(toy_lse))
    expect_equal(out$group, rep(c("C1", "C2"), each = 50))
    expect_type(S4Vectors::metadata(out), "list")
    expect_equal(prod(unique(S4Vectors::metadata(out)$FC)), 1)
    expect_s3_class(S4Vectors::metadata(out)$branch, "data.frame")
    expect_equal(S4Vectors::metadata(out)$branch$A, 18)
    expect_equal(S4Vectors::metadata(out)$branch$B, 13)
    expect_true(round(S4Vectors::metadata(out)$branch$ratio, 5) %in%
                    round(S4Vectors::metadata(out)$FC, 5))

    ## Fix branches - node labels
    set.seed(1)
    out <- simData(tree = NULL, data = NULL, obj = toy_lse, assay = "counts",
                   scenario = "BS", from.A = "Node_18", from.B = "Node_19",
                   minTip.A = 0, maxTip.A = Inf, minTip.B = 0, maxTip.B = Inf,
                   minPr.A = 0, maxPr.A = 1, ratio = 4, adjB = NULL,
                   pct = 0.6, nSam = c(50, 50), mu = 10000, size = NULL,
                   n = 1, FUN = sum, message = FALSE)
    expect_s4_class(out, "TreeSummarizedExperiment")
    expect_equal(nrow(out), nrow(toy_lse))
    expect_equal(ncol(out), 100)
    expect_equal(SummarizedExperiment::assayNames(out), "counts")
    expect_equal(TreeSummarizedExperiment::rowTree(out),
                 TreeSummarizedExperiment::rowTree(toy_lse))
    expect_equal(rownames(out), rownames(toy_lse))
    expect_equal(out$group, rep(c("C1", "C2"), each = 50))
    expect_type(S4Vectors::metadata(out), "list")
    expect_equal(prod(unique(S4Vectors::metadata(out)$FC)), 1)
    expect_s3_class(S4Vectors::metadata(out)$branch, "data.frame")
    expect_equal(S4Vectors::metadata(out)$branch$A, "alias_18")
    expect_equal(S4Vectors::metadata(out)$branch$B, "alias_19")
    expect_true(round(S4Vectors::metadata(out)$branch$ratio, 5) %in%
                    round(S4Vectors::metadata(out)$FC, 5))

    ## Fix branches - aliases
    set.seed(1)
    out <- simData(tree = NULL, data = NULL, obj = toy_lse, assay = "counts",
                   scenario = "BS", from.A = "alias_18", from.B = "alias_19",
                   minTip.A = 0, maxTip.A = Inf, minTip.B = 0, maxTip.B = Inf,
                   minPr.A = 0, maxPr.A = 1, ratio = 4, adjB = NULL,
                   pct = 0.6, nSam = c(50, 50), mu = 10000, size = NULL,
                   n = 1, FUN = sum, message = FALSE)
    expect_s4_class(out, "TreeSummarizedExperiment")
    expect_equal(nrow(out), nrow(toy_lse))
    expect_equal(ncol(out), 100)
    expect_equal(SummarizedExperiment::assayNames(out), "counts")
    expect_equal(TreeSummarizedExperiment::rowTree(out),
                 TreeSummarizedExperiment::rowTree(toy_lse))
    expect_equal(rownames(out), rownames(toy_lse))
    expect_equal(out$group, rep(c("C1", "C2"), each = 50))
    expect_type(S4Vectors::metadata(out), "list")
    expect_equal(prod(unique(S4Vectors::metadata(out)$FC)), 1)
    expect_s3_class(S4Vectors::metadata(out)$branch, "data.frame")
    expect_equal(S4Vectors::metadata(out)$branch$A, "alias_18")
    expect_equal(S4Vectors::metadata(out)$branch$B, "alias_19")
    expect_true(round(S4Vectors::metadata(out)$branch$ratio, 5) %in%
                    round(S4Vectors::metadata(out)$FC, 5))

    ## Multiple matrices
    set.seed(1)
    out <- simData(tree = NULL, data = NULL, obj = toy_lse, assay = "counts",
                   scenario = "BS", from.A = "alias_18", from.B = "alias_19",
                   minTip.A = 0, maxTip.A = Inf, minTip.B = 0, maxTip.B = Inf,
                   minPr.A = 0, maxPr.A = 1, ratio = 4, adjB = NULL,
                   pct = 0.6, nSam = c(50, 50), mu = 10000, size = NULL,
                   n = 2, FUN = sum, message = FALSE)
    expect_s4_class(out, "TreeSummarizedExperiment")
    expect_equal(nrow(out), nrow(toy_lse))
    expect_equal(ncol(out), 100)
    expect_null(SummarizedExperiment::assayNames(out))
    expect_length(SummarizedExperiment::assays(out), 2)
    expect_equal(TreeSummarizedExperiment::rowTree(out),
                 TreeSummarizedExperiment::rowTree(toy_lse))
    expect_equal(rownames(out), rownames(toy_lse))
    expect_equal(out$group, rep(c("C1", "C2"), each = 50))
    expect_type(S4Vectors::metadata(out), "list")
    expect_equal(prod(unique(S4Vectors::metadata(out)$FC)), 1)
    expect_s3_class(S4Vectors::metadata(out)$branch, "data.frame")
    expect_equal(S4Vectors::metadata(out)$branch$A, "alias_18")
    expect_equal(S4Vectors::metadata(out)$branch$B, "alias_19")
    expect_true(round(S4Vectors::metadata(out)$branch$ratio, 5) %in%
                    round(S4Vectors::metadata(out)$FC, 5))

    ## Scenario US; fix branches - node labels
    set.seed(1)
    out <- simData(tree = NULL, data = NULL, obj = toy_lse, assay = "counts",
                   scenario = "US", from.A = "Node_18", from.B = "Node_19",
                   minTip.A = 0, maxTip.A = Inf, minTip.B = 0, maxTip.B = Inf,
                   minPr.A = 0, maxPr.A = 1, ratio = 4, adjB = NULL,
                   pct = 0.6, nSam = c(50, 50), mu = 10000, size = NULL,
                   n = 1, FUN = sum, message = FALSE)
    expect_s4_class(out, "TreeSummarizedExperiment")
    expect_equal(nrow(out), nrow(toy_lse))
    expect_equal(ncol(out), 100)
    expect_equal(SummarizedExperiment::assayNames(out), "counts")
    expect_equal(TreeSummarizedExperiment::rowTree(out),
                 TreeSummarizedExperiment::rowTree(toy_lse))
    expect_equal(rownames(out), rownames(toy_lse))
    expect_equal(out$group, rep(c("C1", "C2"), each = 50))
    expect_type(S4Vectors::metadata(out), "list")
    expect_s3_class(S4Vectors::metadata(out)$branch, "data.frame")
    expect_equal(S4Vectors::metadata(out)$branch$A, "alias_18")
    expect_equal(S4Vectors::metadata(out)$branch$B, "alias_19")

    ## Scenario US; switch branches
    set.seed(1)
    out <- simData(tree = NULL, data = NULL, obj = toy_lse, assay = "counts",
                   scenario = "US", from.A = "Node_19", from.B = "Node_18",
                   minTip.A = 0, maxTip.A = Inf, minTip.B = 0, maxTip.B = Inf,
                   minPr.A = 0, maxPr.A = 1, ratio = 4, adjB = NULL,
                   pct = 0.6, nSam = c(50, 50), mu = 10000, size = NULL,
                   n = 1, FUN = sum, message = FALSE)
    expect_s4_class(out, "TreeSummarizedExperiment")
    expect_equal(nrow(out), nrow(toy_lse))
    expect_equal(ncol(out), 100)
    expect_equal(SummarizedExperiment::assayNames(out), "counts")
    expect_equal(TreeSummarizedExperiment::rowTree(out),
                 TreeSummarizedExperiment::rowTree(toy_lse))
    expect_equal(rownames(out), rownames(toy_lse))
    expect_equal(out$group, rep(c("C1", "C2"), each = 50))
    expect_type(S4Vectors::metadata(out), "list")
    expect_s3_class(S4Vectors::metadata(out)$branch, "data.frame")
    expect_equal(S4Vectors::metadata(out)$branch$A, "alias_19")
    expect_equal(S4Vectors::metadata(out)$branch$B, "alias_18")

    ## Scenario SS; switch branches
    set.seed(1)
    out <- simData(tree = NULL, data = NULL, obj = toy_lse, assay = "counts",
                   scenario = "SS", from.A = "Node_19", from.B = "Node_18",
                   minTip.A = 0, maxTip.A = Inf, minTip.B = 0, maxTip.B = Inf,
                   minPr.A = 0, maxPr.A = 1, ratio = 4, adjB = NULL,
                   pct = 0.6, nSam = c(50, 50), mu = 10000, size = NULL,
                   n = 1, FUN = sum, message = FALSE)
    expect_s4_class(out, "TreeSummarizedExperiment")
    expect_equal(nrow(out), nrow(toy_lse))
    expect_equal(ncol(out), 100)
    expect_equal(SummarizedExperiment::assayNames(out), "counts")
    expect_equal(TreeSummarizedExperiment::rowTree(out),
                 TreeSummarizedExperiment::rowTree(toy_lse))
    expect_equal(rownames(out), rownames(toy_lse))
    expect_equal(out$group, rep(c("C1", "C2"), each = 50))
    expect_type(S4Vectors::metadata(out), "list")
    expect_s3_class(S4Vectors::metadata(out)$branch, "data.frame")
    expect_equal(S4Vectors::metadata(out)$branch$A, "alias_19")
    expect_equal(S4Vectors::metadata(out)$branch$B, "alias_18")

    ## Scenario SS; switch branches - specify size
    set.seed(1)
    out <- simData(tree = NULL, data = NULL, obj = toy_lse, assay = "counts",
                   scenario = "SS", from.A = "Node_19", from.B = "Node_18",
                   minTip.A = 0, maxTip.A = Inf, minTip.B = 0, maxTip.B = Inf,
                   minPr.A = 0, maxPr.A = 1, ratio = 4, adjB = NULL,
                   pct = 0.6, nSam = c(50, 50), mu = 10000, size = 0.5,
                   n = 1, FUN = sum, message = FALSE)
    expect_s4_class(out, "TreeSummarizedExperiment")
    expect_equal(nrow(out), nrow(toy_lse))
    expect_equal(ncol(out), 100)
    expect_equal(SummarizedExperiment::assayNames(out), "counts")
    expect_equal(TreeSummarizedExperiment::rowTree(out),
                 TreeSummarizedExperiment::rowTree(toy_lse))
    expect_equal(rownames(out), rownames(toy_lse))
    expect_equal(out$group, rep(c("C1", "C2"), each = 50))
    expect_type(S4Vectors::metadata(out), "list")
    expect_s3_class(S4Vectors::metadata(out)$branch, "data.frame")
    expect_equal(S4Vectors::metadata(out)$branch$A, "alias_19")
    expect_equal(S4Vectors::metadata(out)$branch$B, "alias_18")

    ## Scenario SS; set adjB
    set.seed(1)
    out <- simData(tree = NULL, data = NULL, obj = toy_lse, assay = "counts",
                   scenario = "SS", from.A = NULL, from.B = NULL,
                   minTip.A = 0, maxTip.A = Inf, minTip.B = 0, maxTip.B = Inf,
                   minPr.A = 0, maxPr.A = 1, ratio = 4, adjB = 0.8,
                   pct = 0.6, nSam = c(50, 50), mu = 10000, size = 0.5,
                   n = 1, FUN = sum, message = FALSE)
    expect_s4_class(out, "TreeSummarizedExperiment")
    expect_equal(nrow(out), nrow(toy_lse))
    expect_equal(ncol(out), 100)
    expect_equal(SummarizedExperiment::assayNames(out), "counts")
    expect_equal(TreeSummarizedExperiment::rowTree(out),
                 TreeSummarizedExperiment::rowTree(toy_lse))
    expect_equal(rownames(out), rownames(toy_lse))
    expect_equal(out$group, rep(c("C1", "C2"), each = 50))
    expect_type(S4Vectors::metadata(out), "list")
    expect_s3_class(S4Vectors::metadata(out)$branch, "data.frame")
    expect_equal(S4Vectors::metadata(out)$branch$A, 14)
    expect_equal(S4Vectors::metadata(out)$branch$B, 15)



})
