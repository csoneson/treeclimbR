# Find the children of an internal node in a tree

Find the direct children of an internal node in a tree.

## Usage

``` r
findChild(tree, node, use.alias = FALSE)
```

## Arguments

- tree:

  A `phylo` object.

- node:

  Either the node number or node label of an internal node of `tree`.

- use.alias:

  A logical scalar. If `FALSE` (default), the node label is used to name
  the output; otherwise, the alias of the node label is used. The alias
  of the node label is created by adding a prefix `"alias_"` to the node
  number.

## Value

A vector of nodes. The numeric value is the node number, and the vector
name is the corresponding node label. If a node has no label, it would
have NA as name when `use.alias = FALSE`, and have the alias of the node
label as name when `use.alias = TRUE`.

## Author

Ruizhu Huang, Charlotte Soneson

## Examples

``` r
suppressPackageStartupMessages({
    library(ggtree)
    library(TreeSummarizedExperiment)
})

data(tinyTree)
ggtree(tinyTree, branch.length = "none") +
    geom_text2(aes(label = node), color = "darkblue",
               hjust = -0.5, vjust = 0.7) +
    geom_text2(aes(label = label), color = "darkorange",
               hjust = -0.1, vjust = -0.7)


## Specify node numbers
findChild(tree = tinyTree, node = c(17, 12))
#> Node_13      t8 
#>      13       6 

## Name return values using aliases
findChild(tree = tinyTree, node = c(17, 12), use.alias = TRUE)
#> alias_13  alias_6 
#>       13        6 

## Specify node labels
findChild(tree = tinyTree, node = c("Node_17", "Node_12"))
#> Node_13      t8 
#>      13       6 

## Tips have no children
findChild(tree = tinyTree, node = "t4")
#> named integer(0)
```
