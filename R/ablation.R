.brainnet_ablation <- function(observed_counts, observed_score, edge_order,
                               assignments, edge_matrix, group_sizes,
                               exact) {
  observed_contributions <- .brainnet_edge_contributions(
    observed_counts,
    group_sizes
  )
  observed_residual <- c(
    observed_score,
    observed_score - cumsum(observed_contributions[edge_order])
  )
  extreme_counts <- integer(length(observed_residual))

  for (assignment_index in seq_len(nrow(assignments))) {
    counts <- .brainnet_group_counts(
      edge_matrix,
      assignments[assignment_index, ],
      n_groups = length(group_sizes)
    )
    contributions <- .brainnet_edge_contributions(counts, group_sizes)
    score <- sum(contributions)
    residual <- c(
      score,
      score - cumsum(contributions[edge_order])
    )
    tie_tolerance <- sqrt(.Machine$double.eps) *
      pmax(1, abs(residual), abs(observed_residual))
    extreme_counts <- extreme_counts +
      (abs(residual) >= abs(observed_residual) - tie_tolerance)
  }

  descriptive_tail_fraction <- if (exact) {
    extreme_counts / nrow(assignments)
  } else {
    (1 + extreme_counts) / (nrow(assignments) + 1)
  }

  data.frame(
    n_removed = 0:length(edge_order),
    statistic = observed_residual,
    descriptive_tail_fraction = descriptive_tail_fraction
  )
}

#' Extract an Exploratory Edge-Ablation Path
#'
#' Returns an explicitly descriptive path computed using an edge order selected
#' from the same data. `descriptive_tail_fraction` is not a post-selection
#' p-value and supports no equivalence or minimality claim.
#'
#' @param x A `brainnet_result` created with `ablation = TRUE`.
#'
#' @return A data frame with the number of removed edges, residual score, and
#'   descriptive permutation-tail fraction.
#' @export
ablation_path <- function(x) {
  .validate_brainnet_result(x)
  if (is.null(x$ablation)) {
    stop(
      "The ablation path was not computed; rerun with `ablation = TRUE`.",
      call. = FALSE
    )
  }
  x$ablation
}
