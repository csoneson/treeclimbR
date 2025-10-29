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
names(out)
#> [1] "edgeR_results" "nodes_drop"    "tree"         
out$nodes_drop
#> [1] "alias_1"  "alias_8"  "alias_12" "alias_18"
edgeR::topTags(out$edgeR_results, sort.by = "PValue")
#> Coefficient:  groupB 
#>               logFC   logCPM          F       PValue          FDR
#> alias_16  1.0468351 16.09379 25.4239062 1.541511e-06 2.312267e-05
#> alias_5   0.7600070 15.55706  8.4473672 4.308004e-03 3.231003e-02
#> alias_6   0.6596231 15.67131  7.3778576 7.517163e-03 3.758582e-02
#> alias_9  -0.1652670 16.97145  0.8615460 3.550504e-01 6.088225e-01
#> alias_10 -0.1684029 16.75636  0.8066779 3.707893e-01 6.088225e-01
#> alias_15 -0.1539229 16.64917  0.7248203 3.961568e-01 6.088225e-01
#> alias_17 -0.1598991 16.60601  0.6872196 4.086549e-01 6.088225e-01
#> alias_7  -0.1492673 16.69156  0.6354511 4.268383e-01 6.088225e-01
#> alias_4  -0.1490846 16.41854  0.5726784 4.505882e-01 6.088225e-01
#> alias_19 -0.1549905 16.43653  0.5647012 4.537526e-01 6.088225e-01
```
