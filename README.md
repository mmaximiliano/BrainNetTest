# BrainNetTest

**BrainNetTest** provides finite-sample inference for independent populations
of aligned binary, undirected brain networks. Its main use is the comparison
of subject-level functional or structural connectomes across clinical,
behavioral, or demographic groups.

The package separates two questions:

1. A whole-graph label-randomization test asks whether subject-level brain
   networks are exchangeable across groups. Its L1 ANOVA score is targeted to
   marginal edge-frequency shifts; it is not an omnibus test for every
   topological change.
2. Edge-wise Fisher exact tests identify marginal frequency differences across
   the complete edge family. Holm adjustment is the default, so selected edges
   have family-wise error control under arbitrary edge dependence.

Validated S3 objects provide a standard R workflow:

`construct -> test -> print/summary -> extract -> plot`

## Installation

```r
install.packages("BrainNetTest")

# Development version
# install.packages("remotes")
remotes::install_github("mmaximiliano/BrainNetTest")
```

## Quick example

```r
library("BrainNetTest")

group_a <- generate_category_graphs(
  n_graphs = 20,
  n_nodes = 10,
  n_communities = 2,
  base_intra_prob = 0.8,
  base_inter_prob = 0.2,
  seed = 1
)
group_b <- generate_category_graphs(
  n_graphs = 20,
  n_nodes = 10,
  n_communities = 2,
  base_intra_prob = 0.5,
  base_inter_prob = 0.5,
  seed = 2
)

networks <- brainnet_data(list(GroupA = group_a, GroupB = group_b))
networks

result <- brainnet_test(
  networks,
  n_permutations = 999,
  adjust = "holm",
  seed = 42
)

result
summary(result)
selected_edges(result)
selected_nodes(result)
plot(result)
```

`as.data.frame(result)` returns all tested edges, including group proportions,
effect sizes, raw p-values, adjusted p-values, and selection decisions.

## Exploratory ablation

```r
result_with_ablation <- brainnet_test(
  networks,
  n_permutations = 999,
  ablation = TRUE,
  seed = 42
)
ablation_path(result_with_ablation)
```

The ablation path is descriptive. Its `descriptive_tail_fraction` is not a
post-selection p-value, and failure to reject after removing an adaptively
ordered prefix is not evidence of equivalence.

## Scope

Version 1.0 supports subject-level binary, undirected brain networks without
self-loops, observed on the same labeled atlas. Atlas selection, registration,
connectivity estimation, nuisance correction, thresholding, and quality
control happen before this package is used. Weighted or directed graphs,
missing regions, paired or clustered observations, and covariate-adjusted
designs are not currently supported.

The statistical contract also applies to other fields where each independent
unit yields a network over the same ordered nodes. Selected edges are group
associations under this contract; they are not, by themselves, causal
connections or validated clinical biomarkers.

## Reference

Fraiman, D. and Fraiman, R. (2018). An ANOVA approach for statistical
comparisons of brain networks. *Scientific Reports*, 8, 4746.
<https://doi.org/10.1038/s41598-018-23152-5>

## License

MIT (c) 2026 Maximiliano Martino and Daniel Fraiman.
