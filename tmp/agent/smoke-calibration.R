suppressPackageStartupMessages(pkgload::load_all(".", quiet = TRUE))

run_cell <- function(n_rep = 500L, p_a = 0.3, p_b = p_a) {
  global_reject <- logical(n_rep)
  edge_family_error <- logical(n_rep)

  set.seed(20260713)
  for (replicate_index in seq_len(n_rep)) {
    group_a <- replicate(
      6L,
      generate_random_graph(5L, p_a),
      simplify = FALSE
    )
    group_b <- replicate(
      6L,
      generate_random_graph(5L, p_b),
      simplify = FALSE
    )
    data <- brainnet_data(list(GroupA = group_a, GroupB = group_b))
    result <- brainnet_test(
      data,
      exact = "never",
      n_permutations = 199L,
      seed = 100000L + replicate_index
    )
    global_reject[[replicate_index]] <- result$global$p_value <= 0.05
    edge_family_error[[replicate_index]] <- any(result$edges$selected)
  }

  c(
    global_rejection = mean(global_reject),
    any_edge_selected = mean(edge_family_error)
  )
}

null <- run_cell()
alternative <- run_cell(p_a = 0.15, p_b = 0.65)

print(rbind(null = null, alternative = alternative))
