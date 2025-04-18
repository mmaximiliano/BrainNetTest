% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rank_edges.R
\name{rank_edges}
\alias{rank_edges}
\title{Rank Edges Based on P-Values}
\usage{
rank_edges(edge_pvalues)
}
\arguments{
\item{edge_pvalues}{A square matrix of p-values for each edge, typically obtained from \code{compute_edge_pvalues}.}
}
\value{
A data frame with edges and their corresponding p-values, ordered from most significant.
}
\description{
This function ranks edges based on their p-values obtained from statistical tests,
ordering them from lowest to highest p-value (most to least significant).
}
\examples{
# Generate synthetic populations
control_graphs <- generate_category_graphs(
  n_graphs = 5,
  n_nodes = 10,
  n_communities = 2,
  base_intra_prob = 0.8,
  base_inter_prob = 0.2,
  seed = 1
)
disease_graphs <- generate_category_graphs(
  n_graphs = 5,
  n_nodes = 10,
  n_communities = 2,
  base_intra_prob = 0.6,
  base_inter_prob = 0.4,
  seed = 2
)
populations <- list(Control = control_graphs, Disease = disease_graphs)

# Compute edge frequencies
frequencies <- compute_edge_frequencies(populations)
edge_counts <- frequencies$edge_counts
N <- sapply(populations, length)

# Compute p-values for edge differences
edge_pvalues <- compute_edge_pvalues(edge_counts, N)

# Rank edges based on p-values
edge_df <- rank_edges(edge_pvalues)
# View the top ranked edges
head(edge_df)
}
