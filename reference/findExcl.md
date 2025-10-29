# Find branches that are non-overlapping with specified branches in a tree

Find all branches whose leaves do not overlap with those of the
specified branches.

## Usage

``` r
findExcl(tree, node, use.alias = FALSE)
```

## Arguments

- tree:

  A `phylo` object.

- node:

  A numeric or character vector specifying the nodes.

- use.alias:

  A logical scalar. If `TRUE`, the alias name is used to name the output
  vector.

## Value

A vector of node numbers

## Author

Ruizhu Huang

## Examples

``` r
suppressPackageStartupMessages({
    library(ggtree)
    library(TreeSummarizedExperiment)
})

data(tinyTree)
ggtree(tinyTree, branch.length = "none") +
    geom_text2(aes(label = node)) +
    geom_hilight(node = 17, fill = "blue", alpha = 0.3) +
    geom_hilight(node = 13, fill = "orange", alpha = 0.3)


## Find branches whose leaves do not overlap with the two colored branches.
## The returned branches are represented at the highest tree level
## possible without including any of the forbidden branches.
findExcl(tree = tinyTree, node = c(17, 13))
#> Node_19      t5      t3 
#>      19       9      10 
```
