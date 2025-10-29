# Simulate different scenarios of abundance change in entities

Simulate a data set with different abundance patterns for entities under
different conditions. These entities have their corresponding nodes on a
tree.

## Usage

``` r
simData(
  tree = NULL,
  data = NULL,
  obj = NULL,
  assay = NULL,
  scenario = "BS",
  from.A = NULL,
  from.B = NULL,
  minTip.A = 0,
  maxTip.A = Inf,
  minTip.B = 0,
  maxTip.B = Inf,
  minPr.A = 0,
  maxPr.A = 1,
  ratio = 4,
  adjB = NULL,
  pct = 0.6,
  nSam = c(50, 50),
  mu = 10000,
  size = NULL,
  n = 1,
  FUN = sum,
  message = FALSE
)
```

## Arguments

- tree:

  A `phylo` object. Only used when `obj` is `NULL`.

- data:

  A count matrix with entities corresponding to tree leaves in the rows
  and samples in the columns. Only used when `obj` is `NULL`.

- obj:

  A `TreeSummarizedExperiment` object with observed data to use as the
  input for the simulation. If `NULL`, `data` and \\ `tree` must be
  provided instead.

- assay:

  If `obj` is not `NULL`, a numeric index or character scalar indicating
  which assay of the object to use as the basis for simulation. If
  `assay` is `NULL`, the first assay in the object is used.

- scenario:

  The simulation scenario, either “BS”, “US”, or “SS” (see **Details**).

- from.A, from.B:

  The branch node labels of branches A and B for which the signal will
  be swapped. By default, both are `NULL`, in which case they will be
  chosen based on the restrictions provided (`minTip.A`, `maxTip.A`,
  `minTip.B`, `maxTip.B`, `minPr.A`, `maxPr.A`, `ratio`). Note: If
  `from.A` is `NULL`, `from.B` is also set to `NULL`.

- minTip.A:

  The minimum number of leaves allowed in branch A.

- maxTip.A:

  The maximum number of leaves allowed in branch A.

- minTip.B:

  The minimum number of leaves allowed in branch B.

- maxTip.B:

  The maximum number of leaves allowed in branch B.

- minPr.A:

  A numeric value in \[0, 1\]. The minimum abundance proportion of
  leaves in branch A.

- maxPr.A:

  A numeric value in \[0, 1\]. The maximum abundance proportion of
  leaves in branch A.

- ratio:

  A numeric value. The proportion ratio of branch B to branch A. This
  value is used to select branches(see **Details**). If there are no
  branches having exactly this ratio, the pair with the value closest to
  `ratio` will be selected.

- adjB:

  A numeric value in \[0, 1\] (only for `scenario` “SS”), or `NULL`. If
  `NULL`, branch A and the selected part of branch B swap their
  proportions. If a numeric value, e.g. 0.1, then the counts for the
  selected part of branch B decreases to 10 the original value, and this
  decrease is added to branch A. For example, assume there are two
  experimental conditions (C1 & C2), branch A has a count of 10 and
  branch B has a count of 40 in C1. If adjB is set to 0.1, then in C2
  branch B becomes 4 and branch A 46 so that the total count of the two
  branches stays the same.

- pct:

  The percentage of leaves in branch B that have differential abundance
  under different conditions (only for scenario “SS”).

- nSam:

  A numeric vector of length 2, indicating the sample size for each of
  the two simulated conditions.

- mu, size:

  The parameters of the Negative Binomial distribution. (see mu and size
  in [`rnbinom`](https://rdrr.io/r/stats/NegBinomial.html)). These
  parameters are used to generate the library size for each simulated
  sample. If `size` is not specified, `mu` should be a vector of numbers
  from which the library size is sampled with replacement.

- n:

  A numeric value to specify how many count tables would be generated
  with the same settings. The default is 1, i.e., one count table would
  be obtained at the end. If greater than 1, the output is a list of
  matrices.

- FUN:

  A function to calculate the aggregated count at each internal node
  based on its descendant leaves (e.g., `sum`, `mean`). The argument of
  the function should be a numeric vector with the counts of an internal
  node's descendant leaves.

- message:

  A logical scalar, indicating whether progress messages should be
  printed to the console.

## Value

a TreeSummarizedExperiment object.

- **assays** A list of count matrices, with entities in rows and samples
  in columns. Each row can be mapped to a node of the tree.

- **rowData** Annotation data for the rows.

- **colData** Annotation data for the columns.

- **rowTree** The tree structure of entities.

- **rowLinks** The link between rows and nodes on the tree.

- **metadata** More details about the simulation.

  - **FC** the fold change of entities corresponding to the tree leaves.

  - **Branch** the information about two selected branches.

    - **A** The branch node label (or number) of branch A.

    - **B** The branch node label (or number) of branch B.

    - **ratio** The count proportion ratio of branch B to branch A.

    - **A_tips** The number of leaves on branch A.

    - **B_tips** The number of leaves on branch B.

    - **A_prop** The count proportion of branch A.

    - **B_prop** The count proportion of branch B.

## Details

Simulate a count table for entities which are corresponding to the nodes
of a tree. The entities are in rows and the samples from different
groups or conditions are in columns. The library size of each sample is
sampled from a Negative Binomial distribution with mean and size
specified by the arguments `mu` and `size`. The counts of entities, that
are mapped to the leaf nodes, in a sample are assumed to follow a
Dirichlet-Multinomial distribution. The parameters for the
Dirichlet-Multinomial distribution are estimated from a real data set
specified by `data` via the function `dirmult` (see
[`dirmult`](https://rdrr.io/pkg/dirmult/man/dirmult.html)). To generate
different abundance patterns under different conditions, we provide
three different scenarios, “BS”, “US”, and “SS” (specified via
`scenario`).

- BS: two branches are selected to swap their proportions, and leaves on
  the same branch have the same fold change.

- US: two branches are selected to swap their proportions. Leaves in the
  same branch have different fold changes but same direction (either
  increase or decrease).

- SS: two branches are selected. One branch has its proportion swapped
  with the proportion of some leaves from the other branch.

## Author

Ruizhu Huang, Charlotte Soneson

## Examples

``` r
suppressPackageStartupMessages({
    library(TreeSummarizedExperiment)
})
## Generate data to use as the starting point (this would usually be a
## real data set)
set.seed(1L)
y <- matrix(rnbinom(120, size = 1, mu = 10), nrow = 10)
colnames(y) <- paste("S", seq_len(12), sep = "")
rownames(y) <- tinyTree$tip.label

toy_lse <- TreeSummarizedExperiment(rowTree = tinyTree,
                                    assays = list(counts = y))
simData(obj = toy_lse, ratio = 2, scenario = "BS", pct = 0.5)
#> Iteration 1: Log-likelihood value: -2142.56163506266
#> Iteration 2: Log-likelihood value: -2140.22864492689
#> Iteration 3: Log-likelihood value: -2139.9586572514
#> Iteration 4: Log-likelihood value: -2139.95301764506
#> Iteration 5: Log-likelihood value: -2139.95301439773
#> class: TreeSummarizedExperiment 
#> dim: 10 100 
#> metadata(3): FC branch scenario
#> assays(1): counts
#> rownames(10): t2 t7 ... t5 t3
#> rowData names(0):
#> colnames(100): C1_1 C1_2 ... C2_49 C2_50
#> colData names(1): group
#> reducedDimNames(0):
#> mainExpName: NULL
#> altExpNames(0):
#> rowLinks: a LinkDataFrame (10 rows)
#> rowTree: 1 phylo tree(s) (10 leaves)
#> colLinks: NULL
#> colTree: NULL
```
