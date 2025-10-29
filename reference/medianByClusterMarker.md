# Calculate median values of markers for each cluster

Calculate median value of each marker in each cluster.

## Usage

``` r
medianByClusterMarker(
  SE,
  assay = 1,
  marker_in_column = TRUE,
  column_cluster = "cluster_id",
  use_marker = NULL
)
```

## Arguments

- SE:

  A
  [`SummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  object.

- assay:

  A numeric index or assay name indicating with assay of `SE` to use to
  calculate medians.

- marker_in_column:

  A logical scalar, indicating whether markers (genes, features) are in
  the columns of `SE` or not.

- column_cluster:

  The name of the column of `colData(SE)` that contains the cluster
  assignment of each sample.

- use_marker:

  A logical or numeric vector such that `SE[use_marker, ]` (if
  `marker_in_column = FALSE`) or `SE[, use_marker]` (if
  `marker_in_column = TRUE`) subsets `SE` to the markers that should be
  retained. If `NULL` (default), all markers are used.

## Value

A
[`SummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
object containing the median value of each marker in each cluster.

## Author

Ruizhu Huang, Charlotte Soneson

## Examples

``` r
suppressPackageStartupMessages({
    library(SummarizedExperiment)
})

## Simulate data with 100 cells and 10 markers (5 type, 5 state markers)
set.seed(1)
count <- matrix(rpois(n = 1000, lambda = 10), nrow = 100)
colnames(count) <- paste0("mk", 1:10)
rowD <- data.frame("cluster" = sample(seq_len(6), 100, replace = TRUE))
colD <- data.frame(type_marker = rep(c(FALSE, TRUE), each = 5))

## SE with markers in columns
d_se <- SummarizedExperiment(assays = list(counts = count),
                             rowData = rowD,
                             colData = colD)
medianByClusterMarker(SE = d_se, marker_in_column = TRUE,
                      column_cluster = "cluster",
                      use_marker = colData(d_se)$type_marker)
#> class: SummarizedExperiment 
#> dim: 6 5 
#> metadata(0):
#> assays(1): counts
#> rownames(6): 1 2 ... 5 6
#> rowData names(1): cluster
#> colnames(5): mk6 mk7 mk8 mk9 mk10
#> colData names(0):

## SE with markers in rows
d_se <- SummarizedExperiment(assays = list(counts = t(count)),
                             rowData = colD,
                             colData = rowD)
medianByClusterMarker(SE = d_se, marker_in_column = FALSE,
                      column_cluster = "cluster",
                      use_marker = rowData(d_se)$type_marker)
#> class: SummarizedExperiment 
#> dim: 5 6 
#> metadata(0):
#> assays(1): counts
#> rownames(5): mk6 mk7 mk8 mk9 mk10
#> rowData names(0):
#> colnames(6): 1 2 ... 5 6
#> colData names(1): cluster
```
