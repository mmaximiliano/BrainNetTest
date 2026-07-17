.edge_test_result <- function(present, group_sizes, method) {
  table <- cbind(
    present = present,
    absent = group_sizes - present
  )

  if (any(colSums(table) == 0L)) {
    return(list(p_value = 1, small_expected = FALSE))
  }

  if (method == "fisher") {
    p_value <- tryCatch(
      stats::fisher.test(table, alternative = "two.sided")$p.value,
      error = function(error) {
        stop(
          "Fisher's exact test failed for edge counts ",
          paste(present, collapse = ", "),
          ": ", conditionMessage(error),
          call. = FALSE
        )
      }
    )
    return(list(p_value = p_value, small_expected = FALSE))
  }

  test <- suppressWarnings(stats::chisq.test(
    table,
    correct = nrow(table) == 2L
  ))
  list(
    p_value = test$p.value,
    small_expected = any(test$expected < 5)
  )
}

.brainnet_edge_inference <- function(counts, group_sizes, edge_indices,
                                     group_names, alpha, method, adjust) {
  group_sizes <- as.integer(group_sizes)
  proportions <- sweep(counts, 2L, group_sizes, "/")
  count_keys <- apply(counts, 1L, paste, collapse = ":")
  unique_keys <- unique(count_keys)
  first_indices <- match(unique_keys, count_keys)
  unique_results <- lapply(
    first_indices,
    function(edge_index) {
      .edge_test_result(counts[edge_index, ], group_sizes, method)
    }
  )
  key_index <- match(count_keys, unique_keys)
  p_values <- vapply(
    unique_results,
    `[[`,
    numeric(1L),
    "p_value"
  )[key_index]
  if (method == "chisq" &&
    any(vapply(unique_results, `[[`, logical(1L), "small_expected"))) {
    warning(
      "Some chi-squared expected counts are below 5; use ",
      "`edge_method = \"fisher\"` for exact edge tests.",
      call. = FALSE
    )
  }
  p_values <- pmin(pmax(p_values, 0), 1)
  p_adjusted <- stats::p.adjust(p_values, method = adjust)

  min_index <- max.col(-proportions, ties.method = "first")
  max_index <- max.col(proportions, ties.method = "first")
  effect_range <- proportions[cbind(seq_len(nrow(proportions)), max_index)] -
    proportions[cbind(seq_len(nrow(proportions)), min_index)]
  risk_difference <- if (length(group_sizes) == 2L) {
    proportions[, 2L] - proportions[, 1L]
  } else {
    rep(NA_real_, nrow(proportions))
  }
  absolute_effect <- if (length(group_sizes) == 2L) {
    abs(risk_difference)
  } else {
    effect_range
  }

  edges <- data.frame(
    node1 = as.integer(edge_indices[, 1L]),
    node2 = as.integer(edge_indices[, 2L]),
    stringsAsFactors = FALSE
  )
  proportion_names <- paste0(
    "proportion_",
    make.unique(make.names(group_names), sep = "_")
  )
  for (group_index in seq_along(group_names)) {
    edges[[proportion_names[[group_index]]]] <- proportions[, group_index]
  }
  edges$risk_difference <- risk_difference
  edges$effect_range <- effect_range
  edges$min_group <- group_names[min_index]
  edges$max_group <- group_names[max_index]
  edges$p_value <- p_values
  edges$p_adjusted <- p_adjusted
  edges$selected <- p_adjusted <= alpha

  order_index <- order(
    edges$p_value,
    -absolute_effect,
    edges$node1,
    edges$node2,
    method = "radix"
  )
  edges$rank <- integer(nrow(edges))
  edges$rank[order_index] <- seq_along(order_index)

  list(
    edges = edges,
    order = order_index,
    family_size = nrow(edges),
    method = method,
    adjustment = adjust,
    alpha = alpha,
    proportion_columns = stats::setNames(
      proportion_names,
      group_names
    )
  )
}
