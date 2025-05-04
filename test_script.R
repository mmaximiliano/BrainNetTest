library(testthat)
library(BrainNetTest)

# Helper function to suppress specific warnings during testing
suppress_warning <- function(expr, pattern) {
  withCallingHandlers(
    expr,
    warning = function(w) {
      if (grepl(pattern, w$message, fixed = TRUE)) {
        invokeRestart("muffleWarning")
      }
    }
  )
}

# Test 1: Basic Functionality: No Crash with Minimal Inputs
test_that("identify_critical_links does not crash with minimal inputs", {
  # Create small adjacency matrices
  set.seed(123)
  GA1 <- generate_random_graph(n_nodes = 3, edge_prob = 0.5)
  GA2 <- generate_random_graph(n_nodes = 3, edge_prob = 0.5)
  GB1 <- generate_random_graph(n_nodes = 3, edge_prob = 0.5)
  GB2 <- generate_random_graph(n_nodes = 3, edge_prob = 0.5)

  # Create populations
  populations <- list(
    A = list(GA1, GA2),
    B = list(GB1, GB2)
  )

  # Run function without error
  result <- expect_no_error(suppress_warning(
    identify_critical_links(populations, alpha = 0.05, method = "fisher"),
    "Initial test is not significant"
    ))

  # Check structure
  expect_type(result, "list")
  expect_true(all(c("critical_edges", "edges_removed", "modified_populations") %in% names(result)))
})

# Test 2: No Edge Removal Expected: Identical Graphs Across Populations
test_that("identify_critical_links issues warning with identical populations", {
  # Generate one graph and replicate it
  set.seed(456)
  G <- generate_random_graph(n_nodes = 5, edge_prob = 0.3)
  A <- list(G, G, G)  # 3 copies
  B <- list(G, G, G)
  populations <- list(A = A, B = B)

  # Run with warning expected
  expect_warning(
    result <- identify_critical_links(populations, alpha = 0.05, method = "fisher"),
    "Initial test is not significant"
  )

  # Check results
  expect_null(result$critical_edges)
  expect_equal(length(result$edges_removed), 0)
  expect_equal(result$modified_populations, populations)
})

# Test 3: Edge Removal Expected: Two Distinct Populations with Clear Differences
test_that("identify_critical_links removes edges with distinct populations", {
  # Generate distinctly different graph populations
  A <- generate_category_graphs(
    n_graphs = 5,
    n_nodes = 5,
    n_communities = 1,
    base_intra_prob = 0.0,
    base_inter_prob = 0.9,
    intra_prob_variation = 0.0,
    inter_prob_variation = 0.0,
    seed = 1
  )

  B <- generate_category_graphs(
    n_graphs = 5,
    n_nodes = 5,
    n_communities = 1,
    base_intra_prob = 0.0,
    base_inter_prob = 0.1,
    intra_prob_variation = 0.0,
    inter_prob_variation = 0.0,
    seed = 2
  )

  populations <- list(A = A, B = B)

  # Use Fisher's exact test which handles small sample issues better
  result <- identify_critical_links(populations, alpha = 0.05, method = "fisher")

  # Check results
  expect_true(!is.null(result$critical_edges))
  expect_true(length(result$edges_removed) > 0)

  # Check edges were actually removed
  if (length(result$edges_removed) > 0) {
    first_edge_removed <- result$edges_removed[[1]]
    node1 <- first_edge_removed[1]
    node2 <- first_edge_removed[2]

    for (pop in result$modified_populations) {
      for (graph in pop) {
        expect_equal(graph[node1, node2], 0)
        expect_equal(graph[node2, node1], 0)
      }
    }
  }
})

# Test 4: Multiple Populations (3 Groups) with Noticeable Differences
test_that("identify_critical_links handles 3 populations", {
  # ----------- generate data -----------
  Control <- generate_category_graphs(
    n_graphs = 20, n_nodes = 12, n_communities = 3,
    base_intra_prob = 0.90, base_inter_prob = 0.05,
    intra_prob_variation = 0.02, inter_prob_variation = 0.02,
    seed = 1
  )

  Disease1 <- generate_category_graphs(
    n_graphs = 20, n_nodes = 12, n_communities = 3,
    base_intra_prob = 0.60, base_inter_prob = 0.30,
    intra_prob_variation = 0.02, inter_prob_variation = 0.02,
    seed = 2
  )

  Disease2 <- generate_category_graphs(
    n_graphs = 20, n_nodes = 12, n_communities = 3,
    base_intra_prob = 0.35, base_inter_prob = 0.55,
    intra_prob_variation = 0.02, inter_prob_variation = 0.02,
    seed = 3
  )

  populations <- list(Control = Control, Disease1 = Disease1, Disease2 = Disease2)

  # ----------- run -----------
  result <- identify_critical_links(populations,
                                    alpha  = 0.05,
                                    method = "fisher",
                                    n_bootstrap = 1000)   # default but explicit

  # ----------- checks -----------
  expect_true(is.data.frame(result$critical_edges) ||
                is.null(result$critical_edges))   # tolerate either branch
  expect_equal(length(result$modified_populations), 3)

  # Optional sanity check – make sure *something* was found
  expect_true(!is.null(result$critical_edges) && nrow(result$critical_edges) > 0)
})


# Test 5 – Minimal edge case: 2‑node graphs, but with enough samples to be significant
test_that("identify_critical_links works with minimal 2‑node graphs", {
  ## ---- Build two extreme populations -----------------------------------
  edge_on  <- matrix(c(0, 1,   # A‑type graphs have the edge
                       1, 0), 2, 2, byrow = TRUE)
  edge_off <- matrix(0, 2, 2)   # B‑type graphs do not

  # Replicate each tiny graph to raise n_i (=> more power for Fisher)
  n_graphs <- 20                # <- key lever
  populations <- list(
    A = replicate(n_graphs, edge_on ,  simplify = FALSE),
    B = replicate(n_graphs, edge_off, simplify = FALSE)
  )

  ## ---- Run -------------------------------------------------------------
  result <- suppressWarnings(
    identify_critical_links(populations,
                            alpha  = 0.05,   # keep default
                            method = "fisher")
  )

  ## ---- Checks ----------------------------------------------------------
  # Fisher should now flag the single edge as significant, so a data‑frame
  # of critical links must be returned.
  expect_true(!is.null(result$critical_edges) &&
                is.data.frame(result$critical_edges) &&
                nrow(result$critical_edges) > 0)

  # That edge must have been removed from *all* graphs
  for (pop in result$modified_populations) {
    for (g in pop) {
      expect_equal(g[1, 2], 0)
      expect_equal(g[2, 1], 0)
    }
  }
})

# Test 6: Larger graphs – stress / scalability with a truly significant difference
test_that("identify_critical_links handles larger graphs and finds significant edges", {
  skip_on_cran()                 # avoid long tests on CRAN

  ## -------- Generate two clearly different 20‑node populations ----------
  n_graphs <- 12                 # more samples -> more power
  n_nodes  <- 20
  n_comms  <- 2

  # Population A : strong community structure (dense intra, sparse inter)
  LargePopA <- generate_category_graphs(
    n_graphs   = n_graphs,
    n_nodes    = n_nodes,
    n_communities = n_comms,
    base_intra_prob = 0.85,
    base_inter_prob = 0.05,
    intra_prob_variation = 0.02,
    inter_prob_variation = 0.02,
    seed = 123
  )

  # Population B : community structure degraded (sparser intra, denser inter)
  LargePopB <- generate_category_graphs(
    n_graphs   = n_graphs,
    n_nodes    = n_nodes,
    n_communities = n_comms,
    base_intra_prob = 0.55,
    base_inter_prob = 0.35,
    intra_prob_variation = 0.02,
    inter_prob_variation = 0.02,
    seed = 456
  )

  populations <- list(A = LargePopA, B = LargePopB)

  ## -------- Run ---------------------------------------------------------
  result <- suppressWarnings(
    identify_critical_links(populations,
                            alpha  = 0.05,
                            method = "fisher",
                            n_bootstrap = 1000)   # default, explicit
  )

  ## -------- Checks ------------------------------------------------------
  expect_type(result, "list")
  expect_true(all(c("critical_edges", "edges_removed", "modified_populations")
                  %in% names(result)))

  # The global test must have been significant, so critical_edges cannot be NULL
  expect_true(is.data.frame(result$critical_edges) &&
                nrow(result$critical_edges) > 0)
})


# Test 7: Testing Different p-Value Adjustment Methods
test_that("identify_critical_links handles different p-value adjustments", {
  # Generate populations
  A <- generate_category_graphs(
    n_graphs = 5,
    n_nodes = 5,
    n_communities = 1,
    base_intra_prob = 0.8,
    base_inter_prob = 0.1,
    seed = 1
  )

  B <- generate_category_graphs(
    n_graphs = 5,
    n_nodes = 5,
    n_communities = 1,
    base_intra_prob = 0.4,
    base_inter_prob = 0.3,
    seed = 2
  )

  populations <- list(A = A, B = B)

  # Run with different adjustment methods
  result_none <- identify_critical_links(
    populations, alpha = 0.05, method = "fisher", adjust_method = "none"
  )

  result_bonferroni <- identify_critical_links(
    populations, alpha = 0.05, method = "fisher", adjust_method = "bonferroni"
  )

  result_bh <- identify_critical_links(
    populations, alpha = 0.05, method = "fisher", adjust_method = "BH"
  )

  # Check all methods completed without error
  expect_true(is.list(result_none))
  expect_true(is.list(result_bonferroni))
  expect_true(is.list(result_bh))

  # Check for expected behavior
  expect_true(is.data.frame(result_none$critical_edges) || is.null(result_none$critical_edges))
  expect_true(is.data.frame(result_bonferroni$critical_edges) || is.null(result_bonferroni$critical_edges))
  expect_true(is.data.frame(result_bh$critical_edges) || is.null(result_bh$critical_edges))
})

# Test 8: DenseVsSparse_3x3
test_that("identify_critical_links works with dense vs sparse 3x3 graphs", {
  # Create fully connected matrices for population A
  A1 <- matrix(c(
    0, 1, 1,
    1, 0, 1,
    1, 1, 0
  ), nrow = 3, byrow = TRUE)

  A2 <- matrix(c(
    0, 1, 1,
    1, 0, 1,
    1, 1, 0
  ), nrow = 3, byrow = TRUE)

  # Create empty matrices for population B
  B1 <- matrix(0, nrow = 3, ncol = 3)
  B2 <- matrix(0, nrow = 3, ncol = 3)

  # Create populations
  populations <- list(A = list(A1, A2), B = list(B1, B2))

  # Run function
  result <- identify_critical_links(populations, alpha = 0.05, method = "fisher")

  # Check results
  expect_true(!is.null(result$critical_edges))
  expect_true(nrow(result$critical_edges) >= 2) # At least 2 nodes should be removed

  # Check that all edges in modified populations are now 0
  for (pop in result$modified_populations) {
    for (graph in pop) {
      expect_equal(sum(graph), 0)
    }
  }
})

# Test 9: SingleEdgeVsNone_2x2
test_that("identify_critical_links works with single edge vs none 2x2", {
  # Create matrices with single edge for population A
  A1 <- matrix(c(
    0, 1,
    1, 0
  ), nrow = 2)
  A2 <- matrix(c(
    0, 1,
    1, 0
  ), nrow = 2)

  # Create empty matrices for population B
  B1 <- matrix(0, nrow = 2, ncol = 2)
  B2 <- matrix(0, nrow = 2, ncol = 2)

  # Create populations
  populations <- list(A = list(A1, A2), B = list(B1, B2))

  # Run function
  result <- identify_critical_links(populations, alpha = 0.05, method = "fisher")

  # Check results
  expect_true(!is.null(result$critical_edges))
  expect_true(nrow(result$critical_edges) == 1) # Single edge should be identified

  # Check that edge (1,2) was removed in population A
  for (pop in result$modified_populations) {
    for (graph in pop) {
      expect_equal(graph[1, 2], 0)
      expect_equal(graph[2, 1], 0)
    }
  }
})

# Test 10: StrongEdgeVsMultipleWeak_5x5
test_that("identify_critical_links identifies strong consistent edge differences", {
  set.seed(123)

  # Create function to generate random adjacency matrices with fixed edge
  generate_with_fixed_edge <- function(n_nodes, edge_prob, i, j, value) {
    G <- generate_random_graph(n_nodes, edge_prob)
    G[i, j] <- G[j, i] <- value
    return(G)
  }

  # Generate population A - always has edge (1,2)
  A <- list()
  for (i in 1:3) {
    A[[i]] <- generate_with_fixed_edge(5, 0.3, 1, 2, 1)
  }

  # Generate population B - never has edge (1,2)
  B <- list()
  for (i in 1:3) {
    B[[i]] <- generate_with_fixed_edge(5, 0.3, 1, 2, 0)
  }

  # Create populations
  populations <- list(A = A, B = B)

  # Run function
  result <- identify_critical_links(populations, alpha = 0.05, method = "fisher")

  # Check results
  expect_true(!is.null(result$critical_edges))

  # Check that edge (1,2) was removed (should be the first or among first removed)
  edge_found <- FALSE
  for (i in seq_along(result$edges_removed)) {
    if ((result$edges_removed[[i]][1] == 1 && result$edges_removed[[i]][2] == 2) ||
        (result$edges_removed[[i]][1] == 2 && result$edges_removed[[i]][2] == 1)) {
      edge_found <- TRUE
      break
    }
  }
  expect_true(edge_found)
})

# Test 11: FullVsSparse_4x4
test_that("identify_critical_links works with full vs Sparse-filled 4x4", {
  # Create fully connected matrices for population A
  A1 <- matrix(1, nrow = 4, ncol = 4)
  diag(A1) <- 0
  A2 <- A1 # same as A1

  # Create half-filled matrices for population B
  set.seed(42)
  B1 <- generate_random_graph(n_nodes = 4, edge_prob = 0.3)
  B2 <- generate_random_graph(n_nodes = 4, edge_prob = 0.3)

  # Create populations
  populations <- list(A = list(A1, A2), B = list(B1, B2))

  # Run function
  result <- identify_critical_links(populations, alpha = 0.05, method = "fisher")

  # Check results
  expect_true(!is.null(result$critical_edges))
  expect_true(length(result$edges_removed) > 0)

  # Check that sufficient edges were removed to make populations similar
  # (we expect the modified A graphs to lose most edges that B doesn't have)
  for (pop_name in names(result$modified_populations)) {
    for (g in seq_along(result$modified_populations[[pop_name]])) {
      graph <- result$modified_populations[[pop_name]][[g]]
      expect_true(sum(graph) < sum(A1)) # Fewer edges than original A
      expect_true(sum(graph) > 1) # At least one edge remains
    }
  }
})

# Test 12: AlmostIdentical_5x5
test_that("identify_critical_links identifies single different edge", {
  set.seed(123)

  # Create base matrix
  baseA <- generate_random_graph(n_nodes = 5, edge_prob = 0.8)

  # Create baseB as copy of baseA, but flip edge (2,3)
  baseB <- baseA
  baseB[2, 3] <- baseB[3, 2] <- 0  # Ensure this edge is different

  # Each population has two copies
  A <- list(baseA, baseA)
  B <- list(baseB, baseB)
  populations <- list(A = A, B = B)

  # Run function
  result <- identify_critical_links(populations, alpha = 0.05, method = "fisher")

  # Check results
  expect_true(!is.null(result$critical_edges))

  # Check that edge (2,3) was removed
  edge_found <- FALSE
  for (i in seq_along(result$edges_removed)) {
    if ((result$edges_removed[[i]][1] == 2 && result$edges_removed[[i]][2] == 3) ||
        (result$edges_removed[[i]][1] == 3 && result$edges_removed[[i]][2] == 2)) {
      edge_found <- TRUE
      break
    }
  }
  expect_true(edge_found)
})

# Test 13: HighProbVsLowProb_4x4
test_that("identify_critical_links works with high vs low probability", {
  set.seed(123)

  # Create high probability matrices for population A
  A <- replicate(3, generate_random_graph(n_nodes = 4, edge_prob = 0.9), simplify = FALSE)

  # Create low probability matrices for population B
  B <- replicate(3, generate_random_graph(n_nodes = 4, edge_prob = 0.2), simplify = FALSE)

  # Create populations
  populations <- list(A = A, B = B)

  # Run function
  result <- identify_critical_links(populations, alpha = 0.05, method = "fisher")

  # Check results
  expect_true(!is.null(result$critical_edges))
  expect_true(length(result$edges_removed) > 0)

  # Check that edges were actually removed
  if (length(result$edges_removed) > 0) {
    first_edge_removed <- result$edges_removed[[1]]
    node1 <- first_edge_removed[1]
    node2 <- first_edge_removed[2]

    for (pop in result$modified_populations) {
      for (graph in pop) {
        expect_equal(graph[node1, node2], 0)
        expect_equal(graph[node2, node1], 0)
      }
    }
  }
})

# Test 14: MissingRowColumn_6x6
test_that("identify_critical_links handles structural differences", {
  # Create function to generate with specific pattern
  create_node1_pattern <- function(n_nodes, node1_connected) {
    G <- matrix(0, nrow = n_nodes, ncol = n_nodes)
    if (node1_connected) {
      G[1, 2:n_nodes] <- 1
      G[2:n_nodes, 1] <- 1
    }
    # Add some random connections between other nodes
    for (i in 2:(n_nodes-1)) {
      for (j in (i+1):n_nodes) {
        G[i, j] <- G[j, i] <- sample(c(0, 1), 1, prob = c(0.7, 0.3))
      }
    }
    return(G)
  }

  # Generate population A with node 1 connected to all others
  A1 <- create_node1_pattern(6, TRUE)
  A2 <- create_node1_pattern(6, TRUE)

  # Generate population B with node 1 disconnected from all others
  B1 <- create_node1_pattern(6, FALSE)
  B2 <- create_node1_pattern(6, FALSE)

  # Create populations
  populations <- list(A = list(A1, A2), B = list(B1, B2))

  # Node 1's connections
  node_1_connections <- max(sum(A1[1, ]), sum(A2[1, ]))

  # Run function
  result <- identify_critical_links(populations, alpha = 0.05, method = "fisher")

  # Check results
  expect_true(!is.null(result$critical_edges))

  # Check that node 1's connections were removed
  for (pop in result$modified_populations) {
    for (graph in pop) {
      # At least some Node 1's connections should removed
      expect_true(sum(graph[1, ]) < node_1_connections)
      expect_true(sum(graph[, 1]) < node_1_connections)
    }
  }
})

# Test 15: DrasticPattern_7x7
test_that("identify_critical_links handles different graph structures", {
  # Create ring structure function
  create_ring_graph <- function(n_nodes) {
    G <- matrix(0, nrow = n_nodes, ncol = n_nodes)
    for (i in 1:n_nodes) {
      # Connect to next node (circular)
      next_node <- if (i == n_nodes) 1 else i + 1
      G[i, next_node] <- G[next_node, i] <- 1

      # Connect to node after next (circular)
      next_next_node <- if (i >= n_nodes - 1) (i + 2) %% n_nodes else i + 2
      if (next_next_node == 0) next_next_node <- n_nodes
      G[i, next_next_node] <- G[next_next_node, i] <- 1
    }
    return(G)
  }

  # Create star structure function
  create_star_graph <- function(n_nodes) {
    G <- matrix(0, nrow = n_nodes, ncol = n_nodes)
    G[1, 2:n_nodes] <- 1
    G[2:n_nodes, 1] <- 1
    return(G)
  }

  # Create ring and star graphs
  A1 <- create_ring_graph(7)
  A2 <- create_ring_graph(7)
  B1 <- create_star_graph(7)
  B2 <- create_star_graph(7)

  # Create populations
  populations <- list(A = list(A1, A2), B = list(B1, B2))

  # Run function
  result <- identify_critical_links(populations, alpha = 0.05, method = "fisher")

  # Check results
  expect_true(!is.null(result$critical_edges))
  expect_true(length(result$edges_removed) > 0)

  # Original structures should be disrupted
  for (pop_name in names(result$modified_populations)) {
    for (g in seq_along(result$modified_populations[[pop_name]])) {
      graph <- result$modified_populations[[pop_name]][[g]]
      # The modified graph shouldn't match either original structure
      if (pop_name == "A") {
        expect_false(identical(graph, A1))
      } else {
        expect_false(identical(graph, B1))
      }
    }
  }
})

# Test 16: MostlySimilarExceptTwo_8x8
test_that("identify_critical_links identifies specific different edges", {
  set.seed(456)

  # Create common base with random edges
  common_base <- generate_random_graph(n_nodes = 8, edge_prob = 0.3)

  # Copy to A_mat and set specific edges to 1
  A_mat <- common_base
  A_mat[2, 3] <- A_mat[3, 2] <- 1
  A_mat[5, 6] <- A_mat[6, 5] <- 1

  # Copy to B_mat and set same edges to 0
  B_mat <- common_base
  B_mat[2, 3] <- B_mat[3, 2] <- 0
  B_mat[5, 6] <- B_mat[6, 5] <- 0

  # Create populations
  populations <- list(A = list(A_mat, A_mat), B = list(B_mat, B_mat))

  # Run function
  result <- identify_critical_links(populations, alpha = 0.05, method = "fisher")

  # Check results
  expect_true(!is.null(result$critical_edges))

  # Check if both specific edges were removed
  edges_found <- c(FALSE, FALSE)
  for (i in seq_along(result$edges_removed)) {
    if ((result$edges_removed[[i]][1] == 2 && result$edges_removed[[i]][2] == 3) ||
        (result$edges_removed[[i]][1] == 3 && result$edges_removed[[i]][2] == 2)) {
      edges_found[1] <- TRUE
    }
    if ((result$edges_removed[[i]][1] == 5 && result$edges_removed[[i]][2] == 6) ||
        (result$edges_removed[[i]][1] == 6 && result$edges_removed[[i]][2] == 5)) {
      edges_found[2] <- TRUE
    }
  }
  expect_true(all(edges_found))
})

# Test 17: TenGraphsOpposite_5x5
test_that("identify_critical_links works with larger populations of different densities", {
  set.seed(789)

  # Create 10 dense graphs for population A
  A <- replicate(10, generate_random_graph(n_nodes = 5, edge_prob = 0.9), simplify = FALSE)

  # Create 10 sparse graphs for population B
  B <- replicate(10, generate_random_graph(n_nodes = 5, edge_prob = 0.1), simplify = FALSE)

  # Create populations
  populations <- list(A = A, B = B)

  # Run function
  result <- identify_critical_links(populations, alpha = 0.05, method = "fisher")

  # Check results
  expect_true(!is.null(result$critical_edges))
  expect_true(length(result$edges_removed) > 0)

  # Check that significant density difference was addressed
  orig_density_diff <- mean(sapply(A, sum)) - mean(sapply(B, sum))
  final_A <- result$modified_populations$A
  final_B <- result$modified_populations$B
  final_density_diff <- mean(sapply(final_A, sum)) - mean(sapply(final_B, sum))

  # Final density difference should be smaller than original
  expect_true(abs(final_density_diff) < abs(orig_density_diff))
})

# Test 18: SingleEdgeDiff_3x3
test_that("identify_critical_links finds single edge difference in 3x3 matrices", {
  # Base adjacency matrix
  base_mat <- matrix(c(
    0, 1, 1,
    1, 0, 0,
    1, 0, 0
  ), nrow = 3, ncol = 3, byrow = TRUE)

  # Population A: Force edge (2,3) = 1
  A1 <- base_mat
  A1[2, 3] <- 1
  A1[3, 2] <- 1  # Keep symmetry
  A2 <- A1  # Identical second graph

  # Population B: Keep that edge at 0
  B1 <- base_mat
  B2 <- base_mat

  # Create populations
  populations <- list(A = list(A1, A2), B = list(B1, B2))

  # Run function
  result <- identify_critical_links(populations, alpha = 0.05, method = "fisher")

  # Check results
  expect_true(!is.null(result$critical_edges))

  # Only one edge should be removed
  expect_equal(nrow(result$critical_edges), 1)
  expect_equal(length(result$edges_removed), 1)

  # Check that edge (2,3) was identified and removed
  edge_removed <- result$edges_removed[[1]]
  edge_found <- ((edge_removed[1] == 2 && edge_removed[2] == 3) ||
                 (edge_removed[1] == 3 && edge_removed[2] == 2))
  expect_true(edge_found)

  # After removal, both populations should be identical
  for (g in seq_along(result$modified_populations$A)) {
    expect_equal(result$modified_populations$A[[g]], result$modified_populations$B[[g]])
  }
})

# Test 19: OneEdgeOff_5x5
test_that("identify_critical_links finds single edge difference in 5x5 matrices", {
  # Create base adjacency matrix with moderate density
  set.seed(999)
  base_mat <- matrix(rbinom(25, 1, 0.4), 5, 5)
  diag(base_mat) <- 0
  # Force symmetry
  base_mat[lower.tri(base_mat)] <- base_mat[upper.tri(base_mat)]

  # Population A: Force edge (1,4) = 1
  A_mat <- base_mat
  A_mat[1, 4] <- 1
  A_mat[4, 1] <- 1  # Keep symmetry
  A <- list(A_mat, A_mat)

  # Population B: Force edge (1,4) = 0
  B_mat <- base_mat
  B_mat[1, 4] <- 0
  B_mat[4, 1] <- 0  # Keep symmetry
  B <- list(B_mat, B_mat)

  # Create populations
  populations <- list(A = A, B = B)

  # Run function
  result <- identify_critical_links(populations, alpha = 0.05, method = "fisher")

  # Check results
  expect_true(!is.null(result$critical_edges))

  # Only one edge should be removed
  expect_equal(nrow(result$critical_edges), 1)
  expect_equal(length(result$edges_removed), 1)

  # Check that edge (1,4) was identified and removed
  edge_removed <- result$edges_removed[[1]]
  edge_found <- ((edge_removed[1] == 1 && edge_removed[2] == 4) ||
                 (edge_removed[1] == 4 && edge_removed[2] == 1))
  expect_true(edge_found)

  # After removal, both populations should be identical
  for (g in seq_along(result$modified_populations$A)) {
    expect_equal(result$modified_populations$A[[g]], result$modified_populations$B[[g]])
  }

  # The edge (1,4) should now be 0 in population A
  expect_equal(result$modified_populations$A[[1]][1, 4], 0)
  expect_equal(result$modified_populations$A[[1]][4, 1], 0)
})
