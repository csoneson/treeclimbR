# Select branches meeting certain criteria

Select branches in a tree meeting the specified criteria in terms of
number of leaves and the count proportion. Note that only internal
branch nodes are considered - no individual leaves will be returned.

## Usage

``` r
selNode(
  pr = NULL,
  obj = NULL,
  assay = 1,
  data = NULL,
  tree = NULL,
  minTip = 0,
  maxTip = Inf,
  minPr = 0,
  maxPr = 1,
  skip = NULL,
  all = FALSE
)
```

## Arguments

- pr:

  A named numeric vector to provide proportions of entities. If this is
  provided, `obj` and `data` will be ignored.

- obj:

  A `TreeSummarizedExperiment` object. Only used if `pr` is `NULL`.

- assay:

  The index or name of the assay of `obj` to use for estimating node
  count proportions. Only used if `obj` is not `NULL`.

- data:

  Either a count table with entities in rows and samples in columns, or
  a list with `pi` and `theta` estimates (the output of
  [`parEstimate`](https://csoneson.github.io/treeclimbR/reference/parEstimate.md)).
  Only used if `pr` and `obj` are `NULL`.

- tree:

  A `phylo` object. If `obj` is used as input, the tree will be
  extracted from the `rowTree` of `obj`.

- minTip:

  the minimum number of leaves in the selected branch.

- maxTip:

  The maximum number of leaves in the selected branch.

- minPr:

  The minimum count proportion of the selected branch in a sample. A
  value between 0 and 1.

- maxPr:

  The maximum count proportion of the selected branch in a sample. A
  value between 0 and 1.

- skip:

  A character vector of node labels. These nodes can not be descendants
  or the ancestors of the selected branch.

- all:

  A logical scalar. If `FALSE` (default), the branch node of a single
  branch, which meets the requirements and has the minimum count
  proportion of branches meeting the requirements, is returned;
  otherwise branch nodes of all branches meeting the requirements are
  returned.

## Value

A `data.frame` with node information for the selected internal node(s).

## Author

Ruizhu Huang, Charlotte Soneson

## Examples

``` r
suppressPackageStartupMessages({
    library(TreeSummarizedExperiment)
})

## Generate example data
set.seed(1)
data(tinyTree)
toyTable <- matrix(rnbinom(40, size = 1, mu = 10), nrow = 10)
colnames(toyTable) <- paste(rep(LETTERS[seq_len(2)], each = 2),
                            rep(seq_len(2), 2), sep = "_")
rownames(toyTable) <- tinyTree$tip.label

## Estimate entity proportions from count matrix under a Dirichlet
## Multinomial framework, and use this as the input for selNode
dat <- parEstimate(obj = toyTable)
#> Iteration 1: Log-likelihood value: -633.210632859034
#> Iteration 2: Log-likelihood value: -631.489192200762
#> Iteration 3: Log-likelihood value: -631.155358434505
#> Iteration 4: Log-likelihood value: -631.132407245137
#> Iteration 5: Log-likelihood value: -631.132252461193
#> Iteration 6: Log-likelihood value: -631.132252452366
selNode(tree = tinyTree, data = dat, all = TRUE)
#>          nodeNum nodeLab proportion numTip
#> alias_11      11 Node_11  1.0000000     10
#> alias_12      12 Node_12  0.8644560      9
#> alias_13      13 Node_13  0.2148379      3
#> alias_14      14 Node_14  0.1488320      2
#> alias_15      15 Node_15  0.6496181      6
#> alias_16      16 Node_16  0.5177530      5
#> alias_17      17 Node_17  0.3484034      3
#> alias_18      18 Node_18  0.1716573      2
#> alias_19      19 Node_19  0.1693496      2
selNode(tree = tinyTree, data = dat,
        minTip = 4, maxTip = 9, minPr = 0, maxPr = 0.8, all = TRUE)
#>          nodeNum nodeLab proportion numTip
#> alias_15      15 Node_15  0.6496181      6
#> alias_16      16 Node_16  0.5177530      5

## Alternatively, directly provide the proportions vector
selNode(tree = tinyTree, pr = dat$pi, all = TRUE)
#>          nodeNum nodeLab proportion numTip
#> alias_11      11 Node_11  1.0000000     10
#> alias_12      12 Node_12  0.8644560      9
#> alias_13      13 Node_13  0.2148379      3
#> alias_14      14 Node_14  0.1488320      2
#> alias_15      15 Node_15  0.6496181      6
#> alias_16      16 Node_16  0.5177530      5
#> alias_17      17 Node_17  0.3484034      3
#> alias_18      18 Node_18  0.1716573      2
#> alias_19      19 Node_19  0.1693496      2

## Return only branch with lowest proportion among valid ones
selNode(tree = tinyTree, pr = dat$pi, all = FALSE)
#>          nodeNum nodeLab proportion numTip
#> alias_14      14 Node_14   0.148832      2

## Start instead from a TreeSummarizedExperiment object
lse <- TreeSummarizedExperiment(rowTree = tinyTree,
                                assays = list(counts = toyTable))
selNode(obj = lse, assay = "counts", all = TRUE)
#> Iteration 1: Log-likelihood value: -633.210632859034
#> Iteration 2: Log-likelihood value: -631.489192200762
#> Iteration 3: Log-likelihood value: -631.155358434505
#> Iteration 4: Log-likelihood value: -631.132407245137
#> Iteration 5: Log-likelihood value: -631.132252461193
#> Iteration 6: Log-likelihood value: -631.132252452366
#>          nodeNum nodeLab proportion numTip
#> alias_11      11 Node_11  1.0000000     10
#> alias_12      12 Node_12  0.8644560      9
#> alias_13      13 Node_13  0.2148379      3
#> alias_14      14 Node_14  0.1488320      2
#> alias_15      15 Node_15  0.6496181      6
#> alias_16      16 Node_16  0.5177530      5
#> alias_17      17 Node_17  0.3484034      3
#> alias_18      18 Node_18  0.1716573      2
#> alias_19      19 Node_19  0.1693496      2

## Don't allow node 1 to be included
selNode(obj = lse, assay = "counts", skip = 1, all = TRUE)
#> Iteration 1: Log-likelihood value: -633.210632859034
#> Iteration 2: Log-likelihood value: -631.489192200762
#> Iteration 3: Log-likelihood value: -631.155358434505
#> Iteration 4: Log-likelihood value: -631.132407245137
#> Iteration 5: Log-likelihood value: -631.132252461193
#> Iteration 6: Log-likelihood value: -631.132252452366
#>          nodeNum nodeLab proportion numTip
#> alias_14      14 Node_14  0.1488320      2
#> alias_15      15 Node_15  0.6496181      6
#> alias_16      16 Node_16  0.5177530      5
#> alias_17      17 Node_17  0.3484034      3
#> alias_18      18 Node_18  0.1716573      2
#> alias_19      19 Node_19  0.1693496      2
```
