# BrainNetTest

**BrainNetTest** is an R package designed to perform hypothesis testing on brain networks represented as graphs. It includes functionalities to compute central graphs, calculate distances using the Manhattan norm, and compute the test statistic \( T \) to determine if different populations of brain networks originate from the same distribution.

## Features

- **Compute Central Graphs:** Create representative central graphs for each population by averaging adjacency matrices.
- **Calculate Distances:** Measure the distance between individual graphs and central graphs using the Manhattan norm.
- **Compute Test Statistic \( T \):** Assess whether the distributions of distances are significantly different across populations.
- **Generate Synthetic Data:** (Optional) Generate random symmetric adjacency matrices for testing purposes.

## Installation

You can install the development version of **BrainNetTest** from GitHub using the `devtools` package:

```r
# Install devtools if not already installed
install.packages("devtools")

# Install BrainNetTest from GitHub
devtools::install_github("mmaximiliano/BrainNetTest")
