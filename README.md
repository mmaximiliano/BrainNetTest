# BrainNetTest

<!-- badges: start -->
<!-- badges: end -->

**BrainNetTest** provides non-parametric hypothesis testing for populations of
brain networks represented as graphs, following the L1-distance ANOVA
framework of Fraiman and Fraiman (2018,
[doi:10.1038/s41598-018-21688-0](https://doi.org/10.1038/s41598-018-21688-0)).

The package includes:

* `compute_central_graph()` and `compute_distance()` for building central
  (mean) graphs and measuring Manhattan (L1) distance between adjacency
  matrices.
* `compute_test_statistic()` for the group test statistic T.
* `identify_critical_links()` for a fast permutation-based identification of
  the edges driving between-group differences, using a prefix-sum
  decomposition that reduces the complexity from O(K * B * |E| * m) to
  O(B * |E| * m).
* `get_critical_nodes()` to summarise the critical edges at the node level.
* `generate_category_graphs()` and `generate_community_graph()` to simulate
  populations of community-structured graphs.
* `plot_critical_edges()` for a multi-panel visualisation of the
  per-population central graphs and the critical edges identified by
  `identify_critical_links()`.

## Installation

```r
# From CRAN (once released)
install.packages("BrainNetTest")

# Development version
# install.packages("remotes")
remotes::install_github("mmaximiliano/BrainNetTest")
```

## Quick example

```r
library(BrainNetTest)

set.seed(1)
control <- generate_category_graphs(
  n_graphs = 20, n_nodes = 10, n_communities = 2,
  base_intra_prob = 0.8, base_inter_prob = 0.2, seed = 1)
patient <- generate_category_graphs(
  n_graphs = 20, n_nodes = 10, n_communities = 2,
  base_intra_prob = 0.6, base_inter_prob = 0.4, seed = 2)

populations <- list(Control = control, Patient = patient)

# Global test and critical-edge identification (permutation test)
result <- identify_critical_links(
  populations, alpha = 0.05, method = "fisher",
  n_permutations = 1000, seed = 42)

head(result$critical_edges)
get_critical_nodes(result)
```

## References

Fraiman, D. and Fraiman, R. (2018) An ANOVA approach for statistical
comparisons of brain networks. *Scientific Reports*, 8, 4746.
<https://doi.org/10.1038/s41598-018-21688-0>

## License

MIT (c) 2026 Maximiliano Martino, Daniel Fraiman.
