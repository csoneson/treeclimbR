test_that("treeScore works", {
    ## Generate some data
    library(TreeSummarizedExperiment)
    library(ggtree)
    library(dplyr)
    data(tinyTree)
    exScore <- data.frame(nodeNum = seq_len(19), score = (seq_len(19))/10)

    ## Check that function returns error for invalid input
    ## -------------------------------------------------------------------------
    expect_error(treeScore(tree = 1, score_data = exScore,
                           node_column = "nodeNum", score_column = "score",
                           new_score = "wScore"),
                 "'tree' must be of class 'phylo'")
    expect_error(treeScore(tree = tinyTree, score_data = 1,
                           node_column = "nodeNum", score_column = "score",
                           new_score = "wScore"),
                 "'score_data' must be of class 'data.frame'")
    expect_error(treeScore(tree = tinyTree, score_data = exScore,
                           node_column = 1, score_column = "score",
                           new_score = "wScore"),
                 "'node_column' must be of class 'character'")
    expect_error(treeScore(tree = tinyTree, score_data = exScore,
                           node_column = c("score", "nodeNum"),
                           score_column = "score", new_score = "wScore"),
                 "'node_column' must have length 1")
    expect_error(treeScore(tree = tinyTree, score_data = exScore,
                           node_column = "missing", score_column = "score",
                           new_score = "wScore"),
                 "All values in 'node_column' must be one of")
    expect_error(treeScore(tree = tinyTree, score_data = exScore,
                           node_column = "nodeNum", score_column = 1,
                           new_score = "wScore"),
                 "'score_column' must be of class 'character'")
    expect_error(treeScore(tree = tinyTree, score_data = exScore,
                           node_column = "nodeNum",
                           score_column = c("score", "nodeNum"),
                           new_score = "wScore"),
                 "'score_column' must have length 1")
    expect_error(treeScore(tree = tinyTree, score_data = exScore,
                           node_column = "nodeNum", score_column = "missing",
                           new_score = "wScore"),
                 "All values in 'score_column' must be one of")
    expect_error(treeScore(tree = tinyTree, score_data = exScore,
                           node_column = "nodeNum", score_column = "score",
                           new_score = 1),
                 "'new_score' must be of class 'character'")
    expect_error(treeScore(tree = tinyTree, score_data = exScore,
                           node_column = "nodeNum", score_column = "score",
                           new_score = c("wScore", "wScore2")),
                 "'new_score' must have length 1")

    ## Check that function works as expected with valid input
    ## -------------------------------------------------------------------------
    newScore <- treeScore(tree = tinyTree, score_data = exScore,
                          node_column = "nodeNum", score_column = "score",
                          new_score = "wScore")
    expect_s3_class(newScore, "data.frame")
    expect_equal(nrow(newScore), nrow(exScore))
    expect_equal(ncol(newScore), ncol(exScore) + 1)
    expect_equal(colnames(newScore), c(colnames(exScore), "wScore"))
    expect_equal(newScore$nodeNum, seq_len(19))
    expect_equal(newScore$score, seq_len(19) / 10)
    expect_equal(round(newScore$wScore, 6),
                 c(seq_len(10) / 10, 1.155781, 1.235069, 0.941667, 0.825,
                   1.434375, 1.4625, 1.325, 1.125, 1.325))

})
