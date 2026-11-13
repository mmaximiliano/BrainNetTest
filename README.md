# BrainNetTest

<!-- badges: start -->
<!-- badges: end -->

**BrainNetTest** provides non-parametric hypothesis testing for populations of
brain networks represented as graphs, following the L1-distance ANOVA
framework of Fraiman and Fraiman (2018,
[doi:10.1038/s41598-018-23152-5](https://doi.org/10.1038/s41598-018-23152-5)).

The package includes:

* `compute_central_graph()` and `compute_distance()` for building central
  (mean) graphs and measuring Manhattan (L1) distance between adjacency
  matrices.
* `compute_test_statistic()` for the group test statistic T, and
  `global_test()` for its permutation p-value: the test of whether the
  populations differ at all.
* `identify_critical_links()` for exploring edges associated with group
  differences. It reuses per-edge contributions and one set of permutations
  throughout the removal steps. The selected edges are exploratory.
  The result is a `critical_links` object with `print()`, `summary()` and
  `plot()` methods.
* `get_critical_nodes()` to summarise the critical edges at the node level.
* `generate_category_graphs()` and `generate_community_graph()` to simulate
  populations of community-structured graphs.
* `plot_critical_edges()` for a multi-panel visualisation of the
  per-population central graphs and the critical edges identified by
  `identify_critical_links()`.

## Installation

```r
# From CRAN
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

# Do the populations differ at all? Statistic T and its permutation p-value
global_test(populations, n_permutations = 1000, seed = 42)

# Explore which edges are associated with the observed difference
result <- identify_critical_links(
  populations, alpha = 0.05, method = "fisher",
  n_permutations = 1000, seed = 42)

result                  # compact report
summary(result)         # settings and node-level ranking
head(result$critical_edges)
get_critical_nodes(result)
```

## Interpreting results

The initial permutation test requires exchangeable group labels under the
null: shuffling whole graphs must preserve the joint distribution of the
observations. Independent graphs from the same graph distribution meet this
condition. Equal edge means alone do not. Paired data, repeated measures,
blocks and covariates need a different permutation design.

The test uses a permutation distribution, not a normal approximation.
`criterion = "combined"` (the default) combines T with the smallest marginal
edge p-value. Use `criterion = "original"` to test with T alone, and report
the criterion, number of permutations, seed and package version.

Report the initial global p-value separately from the selected edges.
Edge removal is an exploratory procedure: it does not control the false
discovery rate or family-wise error rate of the selected set. Adjusting the
marginal p-values with `adjust_method` does not change that, because the
adjustment preserves their order. A non-significant
stopping result does not prove that the remaining networks are equal.

Distances sum both matrix triangles. For symmetric graphs with zero diagonal,
the reported T is twice its value under an upper-triangle-only convention;
this common scale factor does not change the permutation p-value.

For n graphs, m groups, E possible edges and B permutations, the dense
permuted-count calculation costs O(m * n * B * E). It is done once, instead
of repeating the distance calculation at every removal step. Ranking and
marginal tests add work: the balanced two-group Fisher lookup can cost
O(n^3). See `?identify_critical_links` for time and memory details.

## References

Fraiman, D. and Fraiman, R. (2018) An ANOVA approach for statistical
comparisons of brain networks. *Scientific Reports*, 8, 4746.
<https://doi.org/10.1038/s41598-018-23152-5>

## License

MIT (c) 2026 Maximiliano Martino, Daniel Fraiman.
