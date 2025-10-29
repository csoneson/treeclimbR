# Parameter estimation for Dirichlet-multinomial distribution

`parEstimate` is a wrapper of the function
[`dirmult`](https://rdrr.io/pkg/dirmult/man/dirmult.html) with default
settings for `init`, `initscalar`, `epsilon`, `trace` and `mode`. It
allows the input `obj` to be either a `matrix` or a
`TreeSummarizedExperiment` and outputs the estimated values of `pi` and
`theta`.

## Usage

``` r
parEstimate(obj, assay = NULL)
```

## Arguments

- obj:

  A matrix or `TreeSummarizedExperiment`, with samples in the columns
  and entities in the rows.

- assay:

  If `obj` is a `TreeSummarizedExperiment`, the name or index of the
  assay to use to estimate Dirichlet multinomial parameters. If `NULL`,
  the first assay will be used.

## Value

A list including the estimates of `pi` (a vector with one element per
row in `obj`) and `theta` (a scalar).

## Author

Ruizhu Huang, Charlotte Soneson

## Examples

``` r
suppressPackageStartupMessages({
    library(TreeSummarizedExperiment)
})

set.seed(1L)
y <- matrix(rnbinom(200, size = 1, mu = 10), nrow = 10)
colnames(y) <- paste("S", seq_len(20), sep = "")
rownames(y) <- tinyTree$tip.label
toy_tse <- TreeSummarizedExperiment(rowTree = tinyTree,
                                    assays = list(y))
res <- parEstimate(obj = toy_tse)
#> Iteration 1: Log-likelihood value: -3707.79061811115
#> Iteration 2: Log-likelihood value: -3705.14216240896
#> Iteration 3: Log-likelihood value: -3704.92754722943
#> Iteration 4: Log-likelihood value: -3704.92506122156
#> Iteration 5: Log-likelihood value: -3704.92506080308
metadata(res)$assays.par
#> $pi
#>         t2         t7         t6         t9         t4         t8        t10 
#> 0.10811579 0.06185462 0.10655326 0.09701295 0.12987201 0.10127906 0.09596747 
#>         t1         t5         t3 
#> 0.08350587 0.10676714 0.10907184 
#> 
#> $theta
#> [1] 0.1016139
#> 
```
