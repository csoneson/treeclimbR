# Wrapper applying an edgeR differential analysis workflow

`edgerWrp` is a wrapper using functions from the
[`edgeR`](https://rdrr.io/pkg/edgeR/man/edgeR-package.html) package
(Robinson et al. 2010, *Bioinformatics*; McCarthy et al. 2012, *Nucleic
Acids Research*) to fit models and perform a moderated test for each
entity.

## Usage

``` r
edgerWrp(
  count,
  lib_size = NULL,
  option = c("glm", "glmQL"),
  design,
  contrast = NULL,
  normalize = TRUE,
  normalize_method = "TMM",
  ...
)
```

## Arguments

- count:

  A matrix with features (e.g., genes or microbes) in rows and samples
  in columns.

- lib_size:

  A numeric vector with library sizes for each sample. If `NULL`
  (default), the column sums of `count` are used.

- option:

  Either `"glm"` or `"glmQL"`. If `"glm"`,
  [`glmFit`](https://rdrr.io/pkg/edgeR/man/glmfit.html) and
  [`glmLRT`](https://rdrr.io/pkg/edgeR/man/glmLRT.html) are used;
  otherwise, [`glmQLFit`](https://rdrr.io/pkg/edgeR/man/glmQLFit.html)
  and [`glmQLFTest`](https://rdrr.io/pkg/edgeR/man/glmQLFTest.html) are
  used. Details about the difference between the two options can be
  found in the help pages of
  [`glmQLFit`](https://rdrr.io/pkg/edgeR/man/glmQLFit.html).

- design:

  A numeric design matrix, e.g. created by
  [`model.matrix`](https://rdrr.io/r/stats/model.matrix.html). Please
  refer to `design` in
  [`glmQLFit`](https://rdrr.io/pkg/edgeR/man/glmQLFit.html) and
  [`glmFit`](https://rdrr.io/pkg/edgeR/man/glmfit.html) for more
  details.

- contrast:

  A numeric vector specifying one contrast of the linear model
  coefficients to be tested. Its length must equal the number of columns
  of `design`. If `NULL`, the last coefficient will be tested. Please
  refer to `contrast` in
  [`glmQLFTest`](https://rdrr.io/pkg/edgeR/man/glmQLFTest.html) and
  [`glmLRT`](https://rdrr.io/pkg/edgeR/man/glmLRT.html) for more
  details.

- normalize:

  A logical scalar, specifying whether normalization factors should be
  calculated (using
  [`calcNormFactors`](https://rdrr.io/pkg/edgeR/man/calcNormFactors.html)).

- normalize_method:

  Normalization method to be used. Please refer to `method` in
  [`calcNormFactors`](https://rdrr.io/pkg/edgeR/man/calcNormFactors.html)
  for more details.

- ...:

  More arguments to pass to
  [`glmFit`](https://rdrr.io/pkg/edgeR/man/glmfit.html) (if
  `option = "glm"` or
  [`glmQLFit`](https://rdrr.io/pkg/edgeR/man/glmQLFit.html) (if
  `option = "glmQL"`).

## Value

The output of
[`glmQLFTest`](https://rdrr.io/pkg/edgeR/man/glmQLFTest.html) or
[`glmLRT`](https://rdrr.io/pkg/edgeR/man/glmLRT.html) depending on the
specified `option`.

## Details

The function performs the following steps:

- Create a [`DGEList`](https://rdrr.io/pkg/edgeR/man/DGEList.html)
  object. If `lib_size` is given, set the library sizes to these values,
  otherwise use the column sums of the count matrix.

- If `normalize` is `TRUE`, estimate normalization factors using
  [`calcNormFactors`](https://rdrr.io/pkg/edgeR/man/calcNormFactors.html).

- Estimate dispersions with
  [`estimateDisp`](https://rdrr.io/pkg/edgeR/man/estimateDisp.html).

- Depending on the value of `option`, apply either the LRT or QLF edgeR
  workflows (i.e., either
  [`glmFit`](https://rdrr.io/pkg/edgeR/man/glmfit.html) +
  [`glmLRT`](https://rdrr.io/pkg/edgeR/man/glmLRT.html) or
  [`glmQLFit`](https://rdrr.io/pkg/edgeR/man/glmQLFit.html) +
  [`glmQLFTest`](https://rdrr.io/pkg/edgeR/man/glmQLFTest.html)),
  testing for the specified contrast.

## Author

Ruizhu Huang

## Examples

``` r
suppressPackageStartupMessages({
    library(TreeSummarizedExperiment)
})
## Read example data
x <- readRDS(system.file("extdata/da_sim_100_30_18de.rds",
                         package = "treeclimbR"))

## Run differential abundance analysis
out <- edgerWrp(count = assay(x), option = "glm",
                design = model.matrix(~ group, data = colData(x)),
                contrast = c(0, 1))

## The output is an edgeR DGELRT object
class(out)
#> [1] "DGELRT"
#> attr(,"package")
#> [1] "edgeR"
```
