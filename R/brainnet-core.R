.brainnet_encode <- function(data) {
  .validate_brainnet_data(data)
  n_nodes <- data$n_nodes
  edge_indices <- which(
    upper.tri(matrix(FALSE, n_nodes, n_nodes)),
    arr.ind = TRUE
  )
  graphs <- unlist(data$populations, recursive = FALSE)
  n_edges <- nrow(edge_indices)

  encoded_values <- vapply(
    graphs,
    function(graph) as.numeric(graph[edge_indices]),
    numeric(n_edges)
  )
  edge_matrix <- t(matrix(
    encoded_values,
    nrow = n_edges,
    ncol = length(graphs)
  ))
  storage.mode(edge_matrix) <- "double"
  colnames(edge_matrix) <- paste(
    edge_indices[, 1L],
    edge_indices[, 2L],
    sep = "-"
  )

  group <- rep.int(seq_along(data$group_sizes), data$group_sizes)

  list(
    edge_matrix = edge_matrix,
    group = group,
    edge_indices = edge_indices
  )
}

.brainnet_group_counts <- function(edge_matrix, group,
                                   n_groups = max(group)) {
  counts <- vapply(
    seq_len(n_groups),
    function(group_index) {
      colSums(edge_matrix[group == group_index, , drop = FALSE])
    },
    numeric(ncol(edge_matrix))
  )
  if (!is.matrix(counts)) {
    counts <- matrix(counts, ncol = n_groups)
  }
  counts
}

.brainnet_edge_contributions <- function(counts, group_sizes) {
  group_sizes <- as.numeric(group_sizes)
  n_groups <- length(group_sizes)
  n_total <- sum(group_sizes)

  if (!is.matrix(counts) || ncol(counts) != n_groups) {
    stop("`counts` must have one column per group.", call. = FALSE)
  }

  group_proportions <- sweep(counts, 2L, group_sizes, "/")
  pooled_proportions <- rowSums(counts) / n_total
  pooled_matrix <- matrix(
    pooled_proportions,
    nrow = nrow(counts),
    ncol = n_groups
  )

  within <- 2 * group_proportions * (1 - group_proportions)
  pooled <- group_proportions + pooled_matrix -
    2 * group_proportions * pooled_matrix

  within_weights <- sqrt(group_sizes) *
    group_sizes / (group_sizes - 1)
  pooled_weights <- sqrt(group_sizes) *
    n_total / (n_total - 1)

  as.vector(
    sqrt(n_groups) *
      (within %*% within_weights - pooled %*% pooled_weights)
  )
}

.brainnet_score <- function(counts, group_sizes) {
  sum(.brainnet_edge_contributions(counts, group_sizes))
}

.number_group_assignments <- function(group_sizes) {
  log_count <- lgamma(sum(group_sizes) + 1) -
    sum(lgamma(group_sizes + 1))
  if (log_count > log(2^53)) {
    return(Inf)
  }
  round(exp(log_count))
}

.enumerate_group_assignments <- function(group_sizes) {
  group_sizes <- as.integer(group_sizes)
  n_groups <- length(group_sizes)
  n_total <- sum(group_sizes)
  n_assignments <- .number_group_assignments(group_sizes)
  if (!is.finite(n_assignments) || n_assignments > 1000000L) {
    stop(
      "Exhaustive assignment enumeration is limited to 1,000,000 rows.",
      call. = FALSE
    )
  }
  assignments <- vector("list", n_assignments)
  output_index <- 0L

  recurse <- function(remaining, group_index, labels) {
    if (group_index == n_groups) {
      labels[remaining] <- group_index
      output_index <<- output_index + 1L
      assignments[[output_index]] <<- labels
      return(invisible(NULL))
    }

    choices <- utils::combn(
      remaining,
      group_sizes[[group_index]],
      simplify = FALSE
    )
    for (chosen in choices) {
      next_labels <- labels
      next_labels[chosen] <- group_index
      recurse(
        remaining[!remaining %in% chosen],
        group_index + 1L,
        next_labels
      )
    }
    invisible(NULL)
  }

  recurse(seq_len(n_total), 1L, integer(n_total))
  matrix(
    unlist(assignments, use.names = FALSE),
    nrow = length(assignments),
    byrow = TRUE
  )
}

.sample_group_assignments <- function(group_sizes, n_permutations) {
  base_labels <- rep.int(seq_along(group_sizes), group_sizes)
  assignments <- matrix(
    0L,
    nrow = n_permutations,
    ncol = length(base_labels)
  )

  for (index in seq_len(n_permutations)) {
    assignments[index, ] <- sample(base_labels, replace = FALSE)
  }
  assignments
}

.brainnet_scores_for_assignments <- function(edge_matrix, assignments,
                                             group_sizes,
                                             chunk_size = NULL) {
  if (!is.matrix(assignments) ||
    ncol(assignments) != nrow(edge_matrix)) {
    stop("Assignment rows must label every observed graph.", call. = FALSE)
  }
  if (is.null(chunk_size)) {
    target_temporary_bytes <- 64 * 1024^2
    bytes_per_assignment <- max(1, ncol(edge_matrix)) * 8 * 4
    chunk_size <- max(
      1L,
      min(256L, floor(target_temporary_bytes / bytes_per_assignment))
    )
  }
  chunk_size <- .assert_whole_number(chunk_size, "chunk_size")

  group_sizes <- as.numeric(group_sizes)
  n_groups <- length(group_sizes)
  n_total <- sum(group_sizes)
  pooled_proportions <- colSums(edge_matrix) / n_total
  within_weights <- sqrt(group_sizes) *
    group_sizes / (group_sizes - 1)
  pooled_weights <- sqrt(group_sizes) *
    n_total / (n_total - 1)
  output <- numeric(nrow(assignments))

  starts <- seq.int(1L, nrow(assignments), by = chunk_size)
  for (start in starts) {
    finish <- min(start + chunk_size - 1L, nrow(assignments))
    rows <- start:finish
    labels <- assignments[rows, , drop = FALSE]
    chunk_n <- length(rows)
    scores <- numeric(chunk_n)

    for (group_index in seq_len(n_groups)) {
      membership <- t(labels == group_index)
      storage.mode(membership) <- "double"
      counts <- crossprod(edge_matrix, membership)
      proportions <- counts / group_sizes[[group_index]]
      pooled <- matrix(
        pooled_proportions,
        nrow = nrow(proportions),
        ncol = chunk_n
      )
      within <- 2 * proportions * (1 - proportions)
      between <- proportions + pooled - 2 * proportions * pooled
      scores <- scores +
        within_weights[[group_index]] * colSums(within) -
        pooled_weights[[group_index]] * colSums(between)
    }

    output[rows] <- sqrt(n_groups) * scores
  }
  output
}
