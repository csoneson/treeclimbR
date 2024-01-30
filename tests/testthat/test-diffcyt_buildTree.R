test_that("buildTree works", {
     ## Function to create random data (one sample)
    library(diffcyt)
     d_random <- function(n = 20000, mean = 0, sd = 1, ncol = 20, cofactor = 5) {
         d <- sinh(matrix(rnorm(n, mean, sd), ncol = ncol)) * cofactor
         colnames(d) <- paste0("marker", sprintf("%02d", 1:ncol))
         d
     }
     ## Create random data (without differential signal)
     set.seed(1L)
     d_input <- list(
         sample1 = d_random(), sample2 = d_random(),
         sample3 = d_random(), sample4 = d_random()
     )
     experiment_info <- data.frame(
         sample_id = factor(paste0("sample", seq_len(4))),
         group_id = factor(c("group1", "group1", "group2", "group2"))
     )
     marker_info <- data.frame(
         channel_name = paste0("channel", sprintf("%03d", seq_len(20))),
         marker_name = paste0("marker", sprintf("%02d", seq_len(20))),
         marker_class = factor(c(rep("type", 10), rep("state", 10)),
                               levels = c("type", "state", "none"))
     )
     d_se <- prepareData(d_input, experiment_info, marker_info)
     d_se <- transformData(d_se)
     d_se <- generateClusters(d_se)

     ## Test that function returns error for misspecified input
     ## -------------------------------------------------------------------------
     expect_error(buildTree(d_se = 1, dist_method = "euclidean",
                            hclust_method = "average"),
                  "'d_se' must be of class 'SummarizedExperiment'")
     expect_error(buildTree(d_se = d_se, dist_method = 1,
                            hclust_method = "average"),
                  "'dist_method' must be of class 'character'")
     expect_error(buildTree(d_se = d_se, dist_method = c("euclidean", "cosine"),
                            hclust_method = "average"),
                  "'dist_method' must have length 1")
     expect_error(buildTree(d_se = d_se, dist_method = "missing",
                            hclust_method = "average"),
                  "invalid distance method")
     expect_error(buildTree(d_se = d_se, dist_method = "euclidean",
                            hclust_method = 1),
                  "'hclust_method' must be of class 'character'")
     expect_error(buildTree(d_se = d_se, dist_method = "euclidean",
                            hclust_method = c("average", "complete")),
                  "'hclust_method' must have length 1")
     expect_error(buildTree(d_se = d_se, dist_method = "euclidean",
                            hclust_method = "missing"),
                  "invalid clustering method missing")

     ## Test that the function works as expected
     ## ------------------------------------------------------------------------
     tr <- buildTree(d_se, dist_method = "euclidean", hclust_method = "average")
     expect_s3_class(tr, "phylo")
     expect_equal(tr$tip.label, as.character(seq_len(100)))
     expect_equal(tr$node.label, paste0("Node_", seq(from = 101, to = 199,
                                                     by = 1)))
     expect_equal(tr$Nnode, 99)
})
