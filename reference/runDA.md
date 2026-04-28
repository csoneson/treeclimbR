# Test for differential abundance using edgeR

Test for differential abundance of entities using functions from the
[`edgeR`](https://rdrr.io/pkg/edgeR/man/edgeR-package.html) package.
This adapts
[`edgerWrp`](https://csoneson.github.io/treeclimbR/reference/edgerWrp.md)
to accept input as a
[`TreeSummarizedExperiment`](https://rdrr.io/pkg/TreeSummarizedExperiment/man/TreeSummarizedExperiment-constructor.html)
(TSE) object instead of a `matrix`. Features could be represented in
either rows or columns. By default, features are in the rows. Then,
samples are in columns and the sample information is in `colData`. The
tree that stores the hierarchical information about features is in
`rowTree`. Each row of the `assays` can be mapped to a node of the tree.
Data on rows that are mapped to internal nodes is generated from data on
leaf nodes. Normalization for samples is automatically performed by
`edgeR` and the library size is calculated using features that are
mapped to leaf nodes.

## Usage

``` r
runDA(
  TSE,
  feature_on_row = TRUE,
  assay = NULL,
  option = c("glm", "glmQL"),
  design = NULL,
  contrast = NULL,
  filter_min_count = 10,
  filter_min_total_count = 15,
  filter_large_n = 10,
  filter_min_prop = 0.7,
  normalize = TRUE,
  normalize_method = "TMM",
  group_column = "group",
  design_terms = "group",
  ...
)
```

## Arguments

- TSE:

  A `TreeSummarizedExperiment` object.

- feature_on_row:

  A logical scalar. If `TRUE` (default), features or entities (e.g.
  genes, OTUs) are in rows of the `assays` tables, and samples are in
  columns; otherwise, it's the other way around.

- assay:

  A numeric index or assay name to specify which assay from `assays` is
  used for analysis.

- option:

  Either `"glm"` or `"glmQL"`. If `"glm"`,
  [`glmFit`](https://rdrr.io/pkg/edgeR/man/glmfit.html) and
  [`glmLRT`](https://rdrr.io/pkg/edgeR/man/glmLRT.html) are used;
  otherwise, [`glmQLFit`](https://rdrr.io/pkg/edgeR/man/glmQLFit.html)
  and [`glmQLFTest`](https://rdrr.io/pkg/edgeR/man/glmQLFTest.html) are
  used. Details about the difference between two options are in the help
  page of [`glmQLFit`](https://rdrr.io/pkg/edgeR/man/glmQLFit.html).

- design:

  A numeric design matrix. If `NULL`, all columns of the sample
  annotation will be used to create the design matrix.

- contrast:

  A numeric vector specifying one contrast of the linear model
  coefficients to be tested equal to zero. Its length must equal to the
  number of columns of design. If `NULL`, the last coefficient will be
  tested equal to zero.

- filter_min_count:

  A numeric value, passed to **min.count** of
  [`filterByExpr`](https://rdrr.io/pkg/edgeR/man/filterByExpr.html).

- filter_min_total_count:

  A numeric value, passed to **min.total.count** of
  [`filterByExpr`](https://rdrr.io/pkg/edgeR/man/filterByExpr.html).

- filter_large_n:

  A numeric value, passed to **large.n** of
  [`filterByExpr`](https://rdrr.io/pkg/edgeR/man/filterByExpr.html).

- filter_min_prop:

  A numeric value, passed to **min.prop** of
  [`filterByExpr`](https://rdrr.io/pkg/edgeR/man/filterByExpr.html).

- normalize:

  A logical scalar indicating whether to estimate normalization factors
  (using
  [`calcNormFactors`](https://rdrr.io/pkg/edgeR/man/calcNormFactors.html)).

- normalize_method:

  Normalization method to be used. See
  [`calcNormFactors`](https://rdrr.io/pkg/edgeR/man/calcNormFactors.html)
  for more details.

- group_column:

  The name of the column in the sample annotation providing group labels
  for samples (currently not used).

- design_terms:

  The names of columns from the sample annotation that will be used to
  generate the design matrix. This is ignored if **design** is provided.

- ...:

  More arguments to pass to
  [`glmFit`](https://rdrr.io/pkg/edgeR/man/glmfit.html)
  (`option = "glm"` or
  [`glmQLFit`](https://rdrr.io/pkg/edgeR/man/glmQLFit.html)
  (`option = "glmQL"`).

## Value

A list with entries **edgeR_results**, **tree**, and **nodes_drop**.

- edgeR_results:

  The output of
  [`glmQLFTest`](https://rdrr.io/pkg/edgeR/man/glmQLFTest.html) or
  [`glmLRT`](https://rdrr.io/pkg/edgeR/man/glmLRT.html) depending on the
  specified `option`.

- tree:

  The hierarchical structure of entities that was stored in the input
  `TSE`.

- nodes_drop:

  A vector storing the alias node labels of entities that are filtered
  before analysis due to low counts.

## Details

The experimental design is specified by a design matrix and provided via
the argument `design`. More details about the calculation of
normalization factor could be found from
[`calcNormFactors`](https://rdrr.io/pkg/edgeR/man/calcNormFactors.html).

## Author

Ruizhu Huang

## Examples

``` r
suppressPackageStartupMessages({
    library(TreeSummarizedExperiment)
})

## Load example data set
lse <- readRDS(system.file("extdata", "da_sim_100_30_18de.rds",
                           package = "treeclimbR"))

## Aggregate counts on internal nodes
nodes <- showNode(tree = tinyTree, only.leaf = FALSE)
tse <- aggTSE(x = lse, rowLevel = nodes)

dd <- model.matrix(~ group, data = colData(tse))
out <- runDA(TSE = tse, feature_on_row = TRUE,
             assay = 1, option = "glmQL",
             design = dd, contrast = NULL,
             normalize = TRUE, filter_min_count = 2)
#> calcNormFactors has been renamed to normLibSizes
names(out)
#> [1] "edgeR_results" "nodes_drop"    "tree"         
out$nodes_drop
#> [1] "alias_1"  "alias_8"  "alias_12" "alias_18"
edgeR::topTags(out$edgeR_results, sort.by = "PValue")
#> Coefficient:  groupB 
#>               logFC   logCPM          F       PValue          FDR
#> alias_16  1.0465856 16.09379 25.4242776 1.552057e-06 2.328085e-05
#> alias_5   0.7599008 15.55706  8.4414625 4.325592e-03 3.244194e-02
#> alias_6   0.6596795 15.67131  7.3818991 7.507354e-03 3.753677e-02
#> alias_9  -0.1659047 16.97145  0.8703796 3.526151e-01 6.079871e-01
#> alias_10 -0.1693171 16.75636  0.8173532 3.676633e-01 6.079871e-01
#> alias_15 -0.1540461 16.64917  0.7279468 3.951532e-01 6.079871e-01
#> alias_17 -0.1601577 16.60601  0.6906133 4.075134e-01 6.079871e-01
#> alias_7  -0.1489379 16.69156  0.6334974 4.275578e-01 6.079871e-01
#> alias_4  -0.1488372 16.41854  0.5717246 4.509738e-01 6.079871e-01
#> alias_19 -0.1549683 16.43653  0.5652093 4.535590e-01 6.079871e-01
```
