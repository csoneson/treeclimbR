# treeclimbR

<!-- badges: start -->
[![R-CMD-check](https://github.com/csoneson/treeclimbR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/csoneson/treeclimbR/actions)
<!-- badges: end -->

`treeclimbR` is an algorithm to pinpoint the optimal data-dependent resolution for interpreting hierarchical hypotheses.

## Installation

``` r
remotes::install_github("csoneson/treeclimbR")
```

## Citation

If you use `treeclimbR`, please cite

Huang R, Soneson C, Germain PL, Schmidt TSB, Mering CV, Robinson MD: [treeclimbR pinpoints the data-dependent resolution of hierarchical hypotheses](https://doi.org/10.1186/s13059-021-02368-1). _Genome Biology_ 22(1):157 (2021). 

## Basic example ([here](https://fionarhuang.github.io/treeclimbR_toy_example/toy_signal.html))

We first generate a random tree with 100 leaves, each leaf representing an entity. 
Then we sample counts of entities in each sample from a multinomial distribution. 
In total, there are 40 samples (20 in group A and the other 20 in group B). 
18 entities are randomly selected to have their counts multiplied by 2 in group A.

To run `treeclimbR`, we perform a Wilcoxon rank sum test on all nodes of the 
tree to obtain P-values and directions of change. Details are [here](https://fionarhuang.github.io/treeclimbR_toy_example/toy_signal.html). 
Candidates proposed at different `t` values are shown in the animation below.

<p align="center"> 
<img src="https://raw.githubusercontent.com/fionarhuang/treeclimbR_toy_example/master/output/signal_cands.gif">
</p>

The heatmap shows counts of entities (rows) in samples (columns) split by groups. 
Branches that include entities with signals are colored in orange.

### treeclimbR vs BH

Nodes identified by `treeclimbR` are compared to those identified by `BH` under FDR 0.05.

<p align="center"> 
<img src="https://raw.githubusercontent.com/fionarhuang/treeclimbR_toy_example/master/output/signal_result.png">
</p>


## More scenarios ([here](https://htmlpreview.github.io/?https://github.com/fionarhuang/treeclimbR_animation/blob/master/docs/index.html))

1. BS

<p align="center">
<img src="https://raw.githubusercontent.com/fionarhuang/treeclimbR_animation/master/output/pk_BS.gif">
</p>

2. SS

<p align="center">
<img src="https://raw.githubusercontent.com/fionarhuang/treeclimbR_animation/master/output/pk_SS.gif">
</p>

## Learn more
More examples of using `treeclimbR` can be found [here](https://github.com/fionarhuang/treeclimbR_article).
