# tests/testthat/test-generate_category_graphs.R

test_that("generate_category_graphs works correctly", {
  set.seed(123)
  n_graphs <- 5
  n_nodes <- 50
  n_communities <- 3
  graphs <- generate_category_graphs(n_graphs = n_graphs, n_nodes = n_nodes, n_communities = n_communities,
                                     base_intra_prob = 0.8, base_inter_prob = 0.2,
                                     intra_prob_variation = 0.05, inter_prob_variation = 0.05)
  
  expect_equal(length(graphs), n_graphs)
  
  for (G in graphs) {
    expect_equal(dim(G), c(n_nodes, n_nodes))
    expect_true(isSymmetric(G))
    expect_equal(diag(G), rep(0, n_nodes))
  }
  
  # Check that community assignments are consistent across graphs
  community_sizes <- NULL
  if (is.null(community_sizes)) {
    base_size <- n_nodes %/% n_communities
    remainder <- n_nodes %% n_communities
    community_sizes <- rep(base_size, n_communities)
    if (remainder > 0) {
      community_sizes[1:remainder] <- community_sizes[1:remainder] + 1
    }
  }
  community_assignments <- rep(1:n_communities, times = community_sizes)
  
  # Check that intra-community edges are generally more frequent than inter-community edges
  intra_edge_counts <- numeric(n_graphs)
  inter_edge_counts <- numeric(n_graphs)
  
  for (i in seq_along(graphs)) {
    G <- graphs[[i]]
    intra_edges <- 0
    inter_edges <- 0
    for (c in 1:n_communities) {
      nodes_in_c <- which(community_assignments == c)
      intra_edges <- intra_edges + sum(G[nodes_in_c, nodes_in_c])
    }
    inter_edges <- sum(G) - intra_edges
    intra_edge_counts[i] <- intra_edges
    inter_edge_counts[i] <- inter_edges
    expect_true(intra_edges > inter_edges)
  }
  
  # Expect some variation in edge counts due to probability variations
  expect_true(sd(intra_edge_counts) > 0)
  expect_true(sd(inter_edge_counts) > 0)
})
