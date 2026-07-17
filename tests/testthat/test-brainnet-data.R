test_that("brainnet_data constructs a validated grouped-network object", {
  graph <- graph_from_edges(3L, rbind(c(1L, 2L), c(2L, 3L)))
  populations <- list(
    GroupA = constant_group(graph, 2L),
    GroupB = constant_group(graph, 3L)
  )

  data <- brainnet_data(
    populations,
    node_labels = c("A", "B", "C"),
    communities = c("front", "front", "back")
  )

  expect_s3_class(data, "brainnet_data")
  expect_identical(data$group_sizes, c(GroupA = 2L, GroupB = 3L))
  expect_identical(data$n_nodes, 3L)
  expect_identical(data$node_labels, c("A", "B", "C"))
  expect_identical(data$communities, c("front", "front", "back"))
})

test_that("brainnet_data requires named groups and at least two observations", {
  graph <- graph_from_edges(3L)

  expect_error(
    brainnet_data(list(constant_group(graph), constant_group(graph))),
    "named"
  )
  expect_error(
    brainnet_data(list(A = list(graph), B = constant_group(graph))),
    "at least 2"
  )
  expect_error(
    brainnet_data(list(A = constant_group(graph))),
    "at least 2 groups"
  )
})

test_that("brainnet_data enforces the binary undirected graph contract", {
  graph <- graph_from_edges(3L)
  valid <- list(A = constant_group(graph), B = constant_group(graph))

  non_square <- matrix(0L, nrow = 2L, ncol = 3L)
  expect_error(
    brainnet_data(list(A = list(non_square, non_square), B = valid$B)),
    "square"
  )

  mismatched <- graph_from_edges(4L)
  expect_error(
    brainnet_data(list(A = valid$A, B = constant_group(mismatched))),
    "same dimensions"
  )

  asymmetric <- graph
  asymmetric[1L, 2L] <- 1L
  expect_error(
    brainnet_data(list(A = list(graph, asymmetric), B = valid$B)),
    "symmetric"
  )

  looped <- graph
  looped[1L, 1L] <- 1L
  expect_error(
    brainnet_data(list(A = list(graph, looped), B = valid$B)),
    "zero diagonal"
  )

  weighted <- graph
  weighted[1L, 2L] <- weighted[2L, 1L] <- 0.5
  expect_error(
    brainnet_data(list(A = list(graph, weighted), B = valid$B)),
    "binary"
  )

  missing <- graph
  missing[1L, 2L] <- missing[2L, 1L] <- NA_real_
  expect_error(
    brainnet_data(list(A = list(graph, missing), B = valid$B)),
    "finite"
  )
})

test_that("brainnet_data validates and preserves node ordering", {
  graph <- graph_from_edges(3L)
  dimnames(graph) <- list(c("A", "B", "C"), c("A", "B", "C"))
  reordered <- graph[c("B", "A", "C"), c("B", "A", "C")]

  expect_error(
    brainnet_data(list(
      GroupA = list(graph, reordered),
      GroupB = constant_group(graph)
    )),
    "node ordering"
  )

  data <- brainnet_data(list(
    GroupA = constant_group(graph),
    GroupB = constant_group(graph)
  ))
  expect_identical(data$node_labels, c("A", "B", "C"))

  expect_error(
    brainnet_data(
      list(GroupA = constant_group(graph), GroupB = constant_group(graph)),
      node_labels = c("A", "A", "C")
    ),
    "unique"
  )
})

test_that("brainnet_data methods reject corrupted derived metadata", {
  data <- tiny_brainnet_data()

  corrupted_nodes <- data
  corrupted_nodes$n_nodes <- 4L
  expect_error(summary(corrupted_nodes), "internally inconsistent")

  corrupted_sizes <- data
  corrupted_sizes$group_sizes[[1L]] <- 99L
  expect_error(print(corrupted_sizes), "internally inconsistent")
})
