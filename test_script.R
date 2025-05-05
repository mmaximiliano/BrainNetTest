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


# Test 7 – p‑value adjustment methods (all must yield significant edges)
test_that("identify_critical_links handles different p‑value adjustments", {

  ## ------------ Strongly separated 10‑node populations -----------------
  n_graphs <- 30              # many graphs => high Fisher power
  n_nodes  <- 10
  n_comms  <- 2

  A <- generate_category_graphs(
    n_graphs = n_graphs,
    n_nodes  = n_nodes,
    n_communities = n_comms,
    base_intra_prob = 0.90,    # almost‑clique within communities
    base_inter_prob = 0.05,    # sparse between communities
    intra_prob_variation  = 0.02,
    inter_prob_variation  = 0.02,
    seed = 1
  )

  B <- generate_category_graphs(
    n_graphs = n_graphs,
    n_nodes  = n_nodes,
    n_communities = n_comms,
    base_intra_prob = 0.30,    # much weaker community structure
    base_inter_prob = 0.60,
    intra_prob_variation  = 0.02,
    inter_prob_variation  = 0.02,
    seed = 2
  )

  populations <- list(A = A, B = B)

  ## ------------ Run with three adjustment methods ----------------------
  result_none <- identify_critical_links(
    populations, alpha = 0.05, method = "fisher", adjust_method = "none"
  )

  result_bonferroni <- identify_critical_links(
    populations, alpha = 0.05, method = "fisher", adjust_method = "bonferroni"
  )

  result_bh <- identify_critical_links(
    populations, alpha = 0.05, method = "fisher", adjust_method = "BH"
  )

  ## ------------ Common assertions --------------------------------------
  for (res in list(result_none, result_bonferroni, result_bh)) {
    expect_type(res, "list")
    expect_true(all(c("critical_edges", "edges_removed", "modified_populations")
                    %in% names(res)))

    # Global test must have been significant → critical_edges is a non‑empty df
    expect_true(is.data.frame(res$critical_edges) &&
                  nrow(res$critical_edges) > 0)
  }
})


# Test 8 – Dense vs sparse 3×3 graphs
test_that("identify_critical_links works with dense vs sparse 3×3 graphs", {

  ## ---------- Build the base 3×3 patterns ------------------------------
  dense <- matrix(c(
    0, 1, 1,
    1, 0, 1,
    1, 1, 0
  ), nrow = 3, byrow = TRUE)

  sparse <- matrix(0, nrow = 3, ncol = 3)

  ## ---------- Replicate to raise n_i (=> much more power) --------------
  n_rep <- 25                    # ← key lever
  populations <- list(
    A = replicate(n_rep, dense ,  simplify = FALSE),
    B = replicate(n_rep, sparse, simplify = FALSE)
  )

  ## ---------- Run -------------------------------------------------------
  result <- suppressWarnings(
    identify_critical_links(populations,
                            alpha  = 0.05,
                            method = "fisher")
  )

  ## ---------- Checks ----------------------------------------------------
  # 1. A significant global test must have produced a non‑empty data‑frame
  expect_true(!is.null(result$critical_edges) &&
                is.data.frame(result$critical_edges) &&
                nrow(result$critical_edges) >= 2)

  # 2. All graphs in *both* populations should now be edge‑free
  for (pop in result$modified_populations) {
    for (graph in pop) {
      expect_equal(sum(graph), 0)
    }
  }
})


# Test 9 – Single edge vs none (2 × 2)
test_that("identify_critical_links works with single edge vs none 2×2", {

  ## -------- Base patterns ---------------------------------------------
  edge_on  <- matrix(c(0, 1,
                       1, 0), nrow = 2, byrow = TRUE)
  edge_off <- matrix(0, 2, 2)

  ## -------- Replicate to raise n_i ------------------------------------
  n_rep <- 30                    # key lever: boosts Fisher power
  populations <- list(
    A = replicate(n_rep, edge_on ,  simplify = FALSE),
    B = replicate(n_rep, edge_off, simplify = FALSE)
  )

  ## -------- Run --------------------------------------------------------
  result <- suppressWarnings(
    identify_critical_links(populations,
                            alpha  = 0.05,
                            method = "fisher")
  )

  ## -------- Checks -----------------------------------------------------
  # Global test must have been significant → non‑NULL, exactly 1 edge
  expect_true(!is.null(result$critical_edges) &&
                is.data.frame(result$critical_edges))
  expect_equal(nrow(result$critical_edges), 1)   # (1, 2) only

  # Edge (1‑2) must be 0 in every modified graph
  for (pop in result$modified_populations) {
    for (graph in pop) {
      expect_equal(graph[1, 2], 0)
      expect_equal(graph[2, 1], 0)
    }
  }
})


# Test 10 – Strong edge vs multiple weak edges (5 × 5)
test_that("identify_critical_links identifies strong consistent edge differences", {
  set.seed(123)

  ## ----- Helper: random graph with (i,j) forced to `value` -------------
  generate_with_fixed_edge <- function(n_nodes, edge_prob, i, j, value) {
    G <- generate_random_graph(n_nodes, edge_prob)
    G[i, j] <- G[j, i] <- value
    G
  }

  n_graphs <- 30           # key lever: Fisher tables 30 + 30 counts
  n_nodes  <- 5
  edge_p   <- 0.25         # modest density so (1,2) stands out

  ## ----- Build populations ---------------------------------------------
  A <- lapply(seq_len(n_graphs),
              \(x) generate_with_fixed_edge(n_nodes, edge_p, 1, 2, 1))  # always edge
  B <- lapply(seq_len(n_graphs),
              \(x) generate_with_fixed_edge(n_nodes, edge_p, 1, 2, 0))  # never edge

  populations <- list(A = A, B = B)

  ## ----- Run -----------------------------------------------------------
  result <- suppressWarnings(
    identify_critical_links(populations,
                            alpha  = 0.05,
                            method = "fisher")
  )

  ## ----- Checks --------------------------------------------------------
  # Non‑empty list with critical edges detected
  expect_true(!is.null(result$critical_edges) &&
                is.data.frame(result$critical_edges) &&
                nrow(result$critical_edges) > 0)

  # Edge (1,2) must be among the identified/removed edges
  edge_found <- any(
    (result$critical_edges[, 1] == 1 & result$critical_edges[, 2] == 2) |
      (result$critical_edges[, 1] == 2 & result$critical_edges[, 2] == 1)
  )
  expect_true(edge_found)

  # Edge (1,2) should be 0 in every modified graph
  for (pop in result$modified_populations) {
    for (graph in pop) {
      expect_equal(graph[1, 2], 0)
      expect_equal(graph[2, 1], 0)
    }
  }
})


# Test 11 – Full vs sparse‑filled 4×4 (revised so sum(graph) > 1)
test_that("identify_critical_links works with full vs sparse‑filled 4×4", {
  ## ---------- Build deterministic populations --------------------------
  full4 <- matrix(1, 4, 4); diag(full4) <- 0      # 12 ones
  sparse4 <- full4
  sparse4[1, 2] <- sparse4[2, 1] <- 0             # remove three edges
  sparse4[1, 3] <- sparse4[3, 1] <- 0
  sparse4[1, 4] <- sparse4[4, 1] <- 0             # 9 ones left

  n_rep <- 10                                     # enough for Fisher power

  populations <- list(
    A = replicate(n_rep, full4 ,  simplify = FALSE),
    B = replicate(n_rep, sparse4, simplify = FALSE)
  )

  ## ---------- Run -------------------------------------------------------
  result <- suppressWarnings(
    identify_critical_links(populations,
                            alpha  = 0.05,
                            method = "fisher")
  )

  ## ---------- Checks ----------------------------------------------------
  expect_true(!is.null(result$critical_edges) &&
                is.data.frame(result$critical_edges) &&
                nrow(result$critical_edges) > 0)       # some edges removed
  expect_true(length(result$edges_removed) > 0)

  # Each modified graph in both populations:
  for (graph_list in result$modified_populations) {
    for (graph in graph_list) {
      expect_true(sum(graph) < sum(full4))   # fewer than 12 edges
      expect_true(sum(graph) > 1)            # at least two edges remain
    }
  }
})



# Test 12 – Almost‑identical 5×5 graphs: one edge differs, guaranteed significant
test_that("identify_critical_links identifies single different edge", {
  set.seed(123)

  ## ------------ Build deterministic base graphs ------------------------
  baseA <- generate_random_graph(n_nodes = 5, edge_prob = 0.8)
  baseA[2, 3] <- baseA[3, 2] <- 1        # ensure the edge is present

  baseB <- baseA
  baseB[2, 3] <- baseB[3, 2] <- 0        # edge absent in B

  ## ------------ Replicate to raise n_i ---------------------------------
  n_rep <- 40                            # key lever: high Fisher power
  populations <- list(
    A = replicate(n_rep, baseA, simplify = FALSE),
    B = replicate(n_rep, baseB, simplify = FALSE)
  )

  ## ------------ Run -----------------------------------------------------
  result <- suppressWarnings(
    identify_critical_links(populations,
                            alpha  = 0.05,
                            method = "fisher")
  )

  ## ------------ Checks --------------------------------------------------
  expect_true(!is.null(result$critical_edges) &&
                is.data.frame(result$critical_edges) &&
                nrow(result$critical_edges) > 0)

  # Edge (2,3) must be among the removed/identified edges
  edge_found <- any(
    (result$critical_edges[, 1] == 2 & result$critical_edges[, 2] == 3) |
      (result$critical_edges[, 1] == 3 & result$critical_edges[, 2] == 2)
  )
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

# Test 14 – Structural difference: node 1 fully connected vs isolated (6 × 6)
test_that("identify_critical_links handles structural differences", {
  set.seed(123)                             # reproducibility

  n_nodes <- 6
  n_rep   <- 40                            # key lever → high Fisher power

  ## ----- Helper: create deterministic pattern --------------------------
  create_node1_pattern <- function(node1_connected) {
    G <- matrix(0, n_nodes, n_nodes)
    if (node1_connected) {
      G[1, 2:n_nodes] <- 1
      G[2:n_nodes, 1] <- 1                # node 1 connected to all others
    }
    # Random edges among nodes 2…6 (same probability in both groups)
    for (i in 2:(n_nodes - 1)) {
      for (j in (i + 1):n_nodes) {
        val <- rbinom(1, 1, 0.3)          # P(edge) = 0.3
        G[i, j] <- G[j, i] <- val
      }
    }
    G
  }

  ## ----- Build populations ---------------------------------------------
  A_graphs <- replicate(n_rep, create_node1_pattern(TRUE ), simplify = FALSE)
  B_graphs <- replicate(n_rep, create_node1_pattern(FALSE), simplify = FALSE)

  populations <- list(A = A_graphs, B = B_graphs)
  node_1_connections <- n_nodes - 1        # = 5 (original degree in A)

  ## ----- Run ------------------------------------------------------------
  result <- suppressWarnings(
    identify_critical_links(populations,
                            alpha  = 0.05,
                            method = "fisher")
  )

  ## ----- Checks ---------------------------------------------------------
  expect_true(!is.null(result$critical_edges) &&
                is.data.frame(result$critical_edges) &&
                nrow(result$critical_edges) > 0)

  # Node‑1 degree must drop strictly below 5 in *every* modified graph of A,
  # and remain ≤ 5 in B (always true).
  for (pop in result$modified_populations) {
    for (graph in pop) {
      expect_true(sum(graph[1, ]) < node_1_connections)
      expect_true(sum(graph[, 1]) < node_1_connections)
    }
  }
})


# Test 15 – Ring vs star structures (7 × 7) – guaranteed significant
test_that("identify_critical_links handles different graph structures", {
  ## ----- Helper: graph generators --------------------------------------
  create_ring_graph <- function(n_nodes) {
    G <- matrix(0, nrow = n_nodes, ncol = n_nodes)
    for (i in seq_len(n_nodes)) {
      nxt  <- if (i == n_nodes) 1 else i + 1
      nxt2 <- if (i >= n_nodes - 1) (i + 2) %% n_nodes else i + 2
      if (nxt2 == 0) nxt2 <- n_nodes
      G[i, nxt ] <- G[nxt , i] <- 1
      G[i, nxt2] <- G[nxt2, i] <- 1
    }
    G
  }

  create_star_graph <- function(n_nodes) {
    G <- matrix(0, nrow = n_nodes, ncol = n_nodes)
    G[1, 2:n_nodes] <- 1
    G[2:n_nodes, 1] <- 1
    G
  }

  ## ----- Build canonical ring/star graphs ------------------------------
  A_proto <- create_ring_graph(7)
  B_proto <- create_star_graph(7)

  ## ----- Replicate to raise sample size (=> Fisher power) --------------
  n_rep <- 30                                    # key lever
  populations <- list(
    A = replicate(n_rep, A_proto, simplify = FALSE),
    B = replicate(n_rep, B_proto, simplify = FALSE)
  )

  ## ----- Run ------------------------------------------------------------
  result <- suppressWarnings(
    identify_critical_links(populations,
                            alpha  = 0.05,
                            method = "fisher")
  )

  ## ----- Checks ---------------------------------------------------------
  expect_true(!is.null(result$critical_edges) &&
                is.data.frame(result$critical_edges) &&
                nrow(result$critical_edges) > 0)

  expect_true(length(result$edges_removed) > 0)

  # Modified graphs must differ from their originals
  for (pop_name in names(result$modified_populations)) {
    for (graph in result$modified_populations[[pop_name]]) {
      if (pop_name == "A") {
        expect_false(identical(graph, A_proto))
      } else {
        expect_false(identical(graph, B_proto))
      }
    }
  }
})


# Test 16 – Mostly similar 8×8 graphs: two specific edges differ
test_that("identify_critical_links identifies specific different edges", {
  set.seed(456)

  ## ---------- Common random base ---------------------------------------
  common_base <- generate_random_graph(n_nodes = 8, edge_prob = 0.3)

  # Force edges (2,3) and (5,6) to differ
  A_mat <- common_base
  A_mat[2, 3] <- A_mat[3, 2] <- 1
  A_mat[5, 6] <- A_mat[6, 5] <- 1

  B_mat <- common_base
  B_mat[2, 3] <- B_mat[3, 2] <- 0
  B_mat[5, 6] <- B_mat[6, 5] <- 0

  ## ---------- Replicate to raise sample size ---------------------------
  n_rep <- 50                              # key lever → high Fisher power
  populations <- list(
    A = replicate(n_rep, A_mat, simplify = FALSE),
    B = replicate(n_rep, B_mat, simplify = FALSE)
  )

  ## ---------- Run -------------------------------------------------------
  result <- suppressWarnings(
    identify_critical_links(populations,
                            alpha  = 0.05,
                            method = "fisher")
  )

  ## ---------- Checks ----------------------------------------------------
  expect_true(!is.null(result$critical_edges) &&
                is.data.frame(result$critical_edges) &&
                nrow(result$critical_edges) > 0)

  # Both specific edges must be among those removed
  edges_found <- c(FALSE, FALSE)
  for (edge in seq_len(nrow(result$critical_edges))) {
    i <- result$critical_edges[edge, 1]
    j <- result$critical_edges[edge, 2]
    if ((i == 2 && j == 3) || (i == 3 && j == 2)) edges_found[1] <- TRUE
    if ((i == 5 && j == 6) || (i == 6 && j == 5)) edges_found[2] <- TRUE
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

# Test 18 – Single edge difference in 3×3 matrices (guaranteed significant)
test_that("identify_critical_links finds single edge difference in 3×3 matrices", {
  ## -------- Base pattern ----------------------------------------------
  base_mat <- matrix(c(
    0, 1, 1,
    1, 0, 0,
    1, 0, 0
  ), nrow = 3, byrow = TRUE)

  # Population A: force edge (2,3) = 1
  A_mat <- base_mat
  A_mat[2, 3] <- A_mat[3, 2] <- 1

  # Population B: edge (2,3) remains 0
  B_mat <- base_mat

  ## -------- Replicate to raise sample size -----------------------------
  n_rep <- 40                                    # key lever → high Fisher power
  populations <- list(
    A = replicate(n_rep, A_mat, simplify = FALSE),
    B = replicate(n_rep, B_mat, simplify = FALSE)
  )

  ## -------- Run --------------------------------------------------------
  result <- suppressWarnings(
    identify_critical_links(populations,
                            alpha  = 0.05,
                            method = "fisher")   # default batch_size = 1
  )

  ## -------- Checks -----------------------------------------------------
  expect_true(!is.null(result$critical_edges) &&
                is.data.frame(result$critical_edges))

  # Exactly one critical edge should be reported and removed
  expect_equal(nrow(result$critical_edges), 1)
  expect_equal(length(result$edges_removed), 1)

  # That edge must be (2,3)
  edge_removed <- result$edges_removed[[1]]
  edge_found <- (edge_removed[1] == 2 && edge_removed[2] == 3) ||
    (edge_removed[1] == 3 && edge_removed[2] == 2)
  expect_true(edge_found)

  # After removal, the two populations must be identical
  for (g in seq_along(result$modified_populations$A)) {
    expect_equal(result$modified_populations$A[[g]],
                 result$modified_populations$B[[g]])
  }
})


# Test 19 – Single edge difference in 5×5 matrices (guaranteed significant)
test_that("identify_critical_links finds single edge difference in 5×5 matrices", {
  ## -------- Base adjacency matrix -------------------------------------
  set.seed(999)
  base_mat <- matrix(rbinom(25, 1, 0.4), 5, 5)
  diag(base_mat) <- 0
  base_mat[lower.tri(base_mat)] <- base_mat[upper.tri(base_mat)]  # symmetry

  # Population A: edge (1,4) = 1
  A_mat <- base_mat
  A_mat[1, 4] <- A_mat[4, 1] <- 1

  # Population B: edge (1,4) = 0
  B_mat <- base_mat
  B_mat[1, 4] <- B_mat[4, 1] <- 0

  ## -------- Replicate to raise sample size ----------------------------
  n_rep <- 50                        # key lever → high Fisher power
  populations <- list(
    A = replicate(n_rep, A_mat, simplify = FALSE),
    B = replicate(n_rep, B_mat, simplify = FALSE)
  )

  ## -------- Run -------------------------------------------------------
  result <- suppressWarnings(
    identify_critical_links(populations,
                            alpha  = 0.05,
                            method = "fisher")   # batch_size = 1 (default)
  )

  ## -------- Checks ----------------------------------------------------
  expect_true(!is.null(result$critical_edges) &&
                is.data.frame(result$critical_edges))

  # Exactly one critical edge should be reported and removed
  expect_equal(nrow(result$critical_edges), 1)
  expect_equal(length(result$edges_removed), 1)

  # That edge must be (1,4)
  edge_removed <- result$edges_removed[[1]]
  edge_found <- ((edge_removed[1] == 1 && edge_removed[2] == 4) ||
                   (edge_removed[1] == 4 && edge_removed[2] == 1))
  expect_true(edge_found)

  # After removal, the two populations must be identical
  for (g in seq_along(result$modified_populations$A)) {
    expect_equal(result$modified_populations$A[[g]],
                 result$modified_populations$B[[g]])
  }

  # And edge (1,4) must now be 0 in the modified A graphs
  expect_equal(result$modified_populations$A[[1]][1, 4], 0)
  expect_equal(result$modified_populations$A[[1]][4, 1], 0)
})

