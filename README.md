# BrainNetTest <img src="man/figures/logo.png" align="right" height="139" alt="" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/mmaximiliano/BrainNetTest/workflows/R-CMD-check/badge.svg)](https://github.com/mmaximiliano/BrainNetTest/actions)
[![CRAN status](https://www.r-pkg.org/badges/version/BrainNetTest)](https://CRAN.R-project.org/package=BrainNetTest)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

**BrainNetTest** is a comprehensive R package for statistical analysis and hypothesis testing of brain network populations. The package treats brain networks as undirected graphs where nodes represent brain regions and edges represent functional or structural connections between regions.

## Overview

BrainNetTest addresses the critical need for robust statistical methods in computational neuroscience by providing tools to:

- **Detect population-level differences** between brain network groups (e.g., healthy controls vs. patients)
- **Identify critical connections** that drive observed differences between populations  
- **Perform edge-level statistical inference** with multiple testing corrections
- **Generate realistic synthetic brain networks** with community structures for method validation

The package implements both global hypothesis testing for overall population differences and localized inference to pinpoint specific connections responsible for group distinctions.

## Key Features

### 🧠 **Population-Level Analysis**
- Compute representative central graphs for each population
- Calculate Manhattan norm distances between networks
- Statistical hypothesis testing using the T-statistic for global population differences

### 🔍 **Edge-Level Inference** 
- Fisher's exact test, chi-squared test, and proportion test for individual connections
- Multiple testing correction methods (Bonferroni, FDR, etc.)
- Ranking and visualization of statistically significant edges

### 🎯 **Critical Link Identification**
- Novel iterative algorithm to identify edges that explain population differences
- Bootstrap-based hypothesis testing with adaptive edge removal
- Comprehensive output including modified networks after critical edge removal

### 📊 **Network Generation**
- Generate random brain networks with specified properties
- Create networks with realistic community structures
- Support for multiple populations with controlled inter-group differences

## Installation

### Development Version

Install the latest development version from GitHub:

```r
# Install remotes if not already available
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

# Install BrainNetTest
remotes::install_github("mmaximiliano/BrainNetTest", build_vignettes = TRUE)
```

### CRAN Version

```r
# Will be available once published to CRAN
install.packages("BrainNetTest")
```

## Quick Start

```r
library(BrainNetTest)

# Generate synthetic populations with different network properties
control_nets <- generate_category_graphs(
  n_graphs = 50, 
  n_nodes = 100, 
  n_communities = 4,
  base_intra_prob = 0.7, 
  base_inter_prob = 0.1
)

patient_nets <- generate_category_graphs(
  n_graphs = 50, 
  n_nodes = 100, 
  n_communities = 4,
  base_intra_prob = 0.5, 
  base_inter_prob = 0.2
)

# Combine into populations list
populations <- list(Control = control_nets, Patient = patient_nets)

# Test for global population differences
T_statistic <- compute_test_statistic(populations, a = 1)
cat("T-statistic:", T_statistic, "\n")

# Identify critical connections
critical_results <- identify_critical_links(
  populations = populations,
  alpha = 0.05,
  method = "fisher",
  n_bootstrap = 1000
)

# View critical edges
print(critical_results$critical_edges)
```

## Documentation

- **Get started** with the [Introduction vignette](https://mmaximiliano.github.io/BrainNetTest/articles/BrainNetTest_Introduction.html)
- **Learn about community structures** in the [Community Structures vignette](https://mmaximiliano.github.io/BrainNetTest/articles/BrainNetTest_CommunityStructures.html)
- **Browse function documentation** at [https://mmaximiliano.github.io/BrainNetTest/](https://mmaximiliano.github.io/BrainNetTest/)

## Use Cases

### Neuroscientific Applications
- **Clinical research**: Compare brain connectivity between healthy controls and patients with neurological/psychiatric disorders
- **Developmental studies**: Analyze changes in brain network organization across age groups
- **Treatment efficacy**: Assess network-level changes following therapeutic interventions

### Methodological Applications  
- **Method validation**: Generate ground-truth data with known differences for testing new analysis approaches
- **Power analysis**: Determine sample sizes needed to detect specific effect sizes
- **Null model generation**: Create appropriate null distributions for statistical testing

## Algorithm Details

The package implements several key algorithms:

1. **Global Population Testing**: Uses a T-statistic based on Manhattan distances from population centroids
2. **Edge-Level Testing**: Multiple statistical tests (Fisher's exact, chi-squared, proportion) with correction for multiple comparisons  
3. **Critical Link Identification**: Iterative bootstrap-based algorithm that removes edges until population differences are no longer significant

## Dependencies

BrainNetTest has minimal dependencies and integrates well with the R ecosystem:

- **Base R**: Core functionality uses only base R
- **stats**: Statistical functions (automatically available)
- **Suggested packages**: `igraph` (network visualization), `knitr` and `rmarkdown` (vignettes), `testthat` (testing)

## Contributing

We welcome contributions! Please see our [Contributing Guidelines](CONTRIBUTING.md) for details on:

- Reporting bugs and feature requests
- Submitting pull requests  
- Code style and testing requirements

## Citation

If you use BrainNetTest in your research, please cite:

```
Martino, M. (2024). BrainNetTest: Statistical Testing and Edge Identification 
in Brain Network Populations. R package version 0.0.0.9000.
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Support

- **Issues**: Report bugs or request features on [GitHub Issues](https://github.com/mmaximiliano/BrainNetTest/issues)
- **Discussions**: Ask questions on [GitHub Discussions](https://github.com/mmaximiliano/BrainNetTest/discussions)
- **Email**: Contact the maintainer at maxii.martino@gmail.com

---

*BrainNetTest is developed and maintained by Maximiliano Martino.*
