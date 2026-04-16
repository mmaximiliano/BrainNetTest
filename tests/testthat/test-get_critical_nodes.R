test_that("get_critical_nodes returns correct structure without labels", {
  skip_on_cran()

  set.seed(42)
  ctrl <- generate_category_graphs(n_graphs = 15, n_nodes = 10,
    n_communities = 2, base_intra_prob = 0.8, base_inter_prob = 0.2, seed = 1)
  dis  <- generate_category_graphs(n_graphs = 15, n_nodes = 10,
    n_communities = 2, base_intra_prob = 0.5, base_inter_prob = 0.5, seed = 2)
  pops <- list(Control = ctrl, Disease = dis)

  result <- identify_critical_links(pops, alpha = 0.05, method = "fisher",
    n_permutations = 200, seed = 99)

  nodes_df <- get_critical_nodes(result)

  expect_true(is.data.frame(nodes_df))
  expect_true("node" %in% names(nodes_df))
  expect_true("critical_degree" %in% names(nodes_df))
  expect_false("label" %in% names(nodes_df))

  # Every node must appear in at least one critical edge
  expect_true(all(nodes_df$critical_degree >= 1))

  # Nodes must be a subset of nodes in critical_edges
  ce <- result$critical_edges
  all_edge_nodes <- unique(c(ce$node1, ce$node2))
  expect_true(all(nodes_df$node %in% all_edge_nodes))
  expect_equal(sort(nodes_df$node), sort(all_edge_nodes))

  # Sorted by critical_degree descending
  expect_true(all(diff(nodes_df$critical_degree) <= 0))
})


test_that("get_critical_nodes includes labels when provided", {
  skip_on_cran()

  set.seed(42)
  ctrl <- generate_category_graphs(n_graphs = 15, n_nodes = 10,
    n_communities = 2, base_intra_prob = 0.8, base_inter_prob = 0.2, seed = 1)
  dis  <- generate_category_graphs(n_graphs = 15, n_nodes = 10,
    n_communities = 2, base_intra_prob = 0.5, base_inter_prob = 0.5, seed = 2)
  pops <- list(Control = ctrl, Disease = dis)

  result <- identify_critical_links(pops, alpha = 0.05, method = "fisher",
    n_permutations = 200, seed = 99)

  labels <- paste0("R", 1:10)
  nodes_df <- get_critical_nodes(result, node_labels = labels)

  expect_true("label" %in% names(nodes_df))
  # Column order: node, label, critical_degree
  expect_equal(names(nodes_df), c("node", "label", "critical_degree"))
  # Labels must match
  for (i in seq_len(nrow(nodes_df))) {
    expect_equal(nodes_df$label[i], labels[nodes_df$node[i]])
  }
})


test_that("get_critical_nodes handles NULL critical_edges", {
  # Simulate a result where no differences were found
  result <- list(
    critical_edges       = NULL,
    edges_removed        = list(),
    modified_populations = list()
  )

  nodes_df <- get_critical_nodes(result)
  expect_true(is.data.frame(nodes_df))
  expect_equal(nrow(nodes_df), 0)
  expect_true("node" %in% names(nodes_df))
  expect_true("critical_degree" %in% names(nodes_df))

  # Also with labels
  nodes_df2 <- get_critical_nodes(result, node_labels = c("A", "B", "C"))
  expect_equal(nrow(nodes_df2), 0)
  expect_true("label" %in% names(nodes_df2))
})


test_that("get_critical_nodes validates inputs", {
  # Not a list
  expect_error(get_critical_nodes("not_a_list"), "identify_critical_links")

  # Missing critical_edges component
  expect_error(get_critical_nodes(list(a = 1)), "identify_critical_links")

  # node_labels not character
  result <- list(critical_edges = data.frame(node1 = 1, node2 = 2, p_value = 0.01))
  expect_error(get_critical_nodes(result, node_labels = 1:5), "character")

  # node_labels too short
  expect_error(get_critical_nodes(result, node_labels = "A"), "length 1")
})


test_that("get_critical_nodes critical_degree sums correctly", {
  # Manual example: edges (1,2), (1,3), (2,3)
  ce <- data.frame(node1 = c(1L, 1L, 2L), node2 = c(2L, 3L, 3L),
                   p_value = c(0.001, 0.01, 0.02))
  result <- list(critical_edges = ce, edges_removed = list(), modified_populations = list())

  nodes_df <- get_critical_nodes(result)

  expect_equal(nrow(nodes_df), 3)
  # Node 1 appears in edges (1,2) and (1,3) → degree 2
  expect_equal(nodes_df$critical_degree[nodes_df$node == 1], 2L)
  # Node 2 appears in edges (1,2) and (2,3) → degree 2
  expect_equal(nodes_df$critical_degree[nodes_df$node == 2], 2L)
  # Node 3 appears in edges (1,3) and (2,3) → degree 2
  expect_equal(nodes_df$critical_degree[nodes_df$node == 3], 2L)
})
