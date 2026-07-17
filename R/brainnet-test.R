#' Test Populations of Binary Brain Networks
#'
#' Performs a whole-graph label-randomization test for independent groups of
#' aligned, binary, undirected brain networks using the unnormalized Fraiman
#' score. It separately tests every edge for marginal frequency differences
#' with multiplicity control.
#'
#' @param x A [brainnet_data()] object.
#' @param alpha Significance level in `(0, 1)`.
#' @param n_permutations Number of random fixed-size label assignments when
#'   Monte Carlo inference is used.
#' @param exact One of `"auto"`, `"always"`, or `"never"`. Automatic mode
#'   exhaustively enumerates at most `max_exact` assignments.
#' @param max_exact Maximum number of assignments allowed for exhaustive
#'   enumeration.
#' @param edge_method Edge-level test, `"fisher"` (default) or `"chisq"`.
#' @param adjust Multiplicity adjustment: `"holm"` (default), `"BY"`, or
#'   `"BH"`.
#' @param ablation Logical; compute the explicitly descriptive edge-ablation
#'   path.
#' @param seed Optional whole-number random seed. The caller's random-number
#'   state is restored on exit.
#'
#' @return An object of class `brainnet_result`. Its `edges` data frame
#'   contains `node1`, `node2`, one uniquely named proportion column per group,
#'   `risk_difference` (group 2 minus group 1 for two groups),
#'   `effect_range`, `min_group`, `max_group`, `p_value`, `p_adjusted`,
#'   `selected`, and deterministic `rank`. The
#'   `edge_inference$proportion_columns` vector maps original group names to
#'   proportion-column names.
#'
#' @details The finite-sample global null is exchangeability of whole-graph
#'   observations with respect to fixed-size group labels. The score is
#'   targeted to marginal edge-frequency changes and is not an omnibus test
#'   for alternatives that preserve every marginal edge probability.
#'   Monte Carlo assignments are sampled uniformly with replacement from the
#'   complete fixed-size orbit, including the possibility of drawing the
#'   observed assignment.
#'
#'   Edge tests are not gated by, and must not be interpreted as a
#'   decomposition of, the global result. Holm adjustment is valid under
#'   arbitrary dependence among edge tests. BY is the corresponding
#'   arbitrary-dependence false-discovery-rate option; BH requires additional
#'   positive-dependence assumptions.
#'
#'   Numerically equivalent permutation statistics are treated as ties using
#'   tolerance `sqrt(.Machine$double.eps) * max(1, abs(T))`, expanded when a
#'   permutation statistic is larger. The realized tolerance is stored in
#'   `result$global$tie_tolerance`.
#'
#' @references
#' Fraiman D, Fraiman R (2018). "An ANOVA approach for statistical comparisons
#' of brain networks." \emph{Scientific Reports}, 8, 4746.
#' \doi{10.1038/s41598-018-23152-5}.
#' @export
brainnet_test <- function(x,
                          alpha = 0.05,
                          n_permutations = 999L,
                          exact = c("auto", "always", "never"),
                          max_exact = 10000L,
                          edge_method = c("fisher", "chisq"),
                          adjust = c("holm", "BY", "BH"),
                          ablation = FALSE,
                          seed = NULL) {
  started <- proc.time()[["elapsed"]]
  .validate_brainnet_data(x)
  .assert_scalar_number(
    alpha,
    "alpha",
    lower = 0,
    upper = 1,
    lower_open = TRUE,
    upper_open = TRUE
  )
  n_permutations <- .assert_whole_number(
    n_permutations,
    "n_permutations"
  )
  max_exact <- .assert_whole_number(max_exact, "max_exact")
  if (max_exact > 1000000L) {
    stop("`max_exact` cannot exceed 1,000,000 assignments.",
      call. = FALSE
    )
  }
  if (missing(exact)) {
    exact <- "auto"
  } else if (length(exact) != 1L ||
    !exact %in% c("auto", "always", "never")) {
    stop(
      "`exact` must be one of \"auto\", \"always\", or \"never\".",
      call. = FALSE
    )
  }
  edge_method <- match.arg(edge_method)
  adjust <- match.arg(adjust)
  if (!is.logical(ablation) || length(ablation) != 1L || is.na(ablation)) {
    stop("`ablation` must be `TRUE` or `FALSE`.", call. = FALSE)
  }
  if (!is.null(seed)) {
    .assert_whole_number(seed, "seed", minimum = 0L)
  }
  if (adjust == "BH") {
    warning(
      "BH adjustment assumes suitable positive dependence between edge ",
      "tests; use Holm or BY when that assumption is not justified.",
      call. = FALSE
    )
  }

  encoded <- .brainnet_encode(x)
  observed_counts <- .brainnet_group_counts(
    encoded$edge_matrix,
    encoded$group,
    n_groups = length(x$group_sizes)
  )
  observed_score <- .brainnet_score(observed_counts, x$group_sizes)
  n_assignments <- .number_group_assignments(x$group_sizes)
  use_exact <- switch(exact,
    always = TRUE,
    never = FALSE,
    auto = is.finite(n_assignments) && n_assignments <= max_exact
  )
  if (use_exact && n_assignments > max_exact) {
    stop(
      "Exhaustive inference needs ", format(n_assignments, scientific = FALSE),
      " assignments, above `max_exact = ", max_exact, "`.",
      call. = FALSE
    )
  }

  permutation <- .with_preserved_seed(seed, function() {
    assignments <- if (use_exact) {
      .enumerate_group_assignments(x$group_sizes)
    } else {
      .sample_group_assignments(
        x$group_sizes,
        n_permutations
      )
    }
    scores <- .brainnet_scores_for_assignments(
      encoded$edge_matrix,
      assignments,
      x$group_sizes
    )
    list(assignments = assignments, scores = scores)
  })

  observed_extremeness <- abs(observed_score)
  permutation_extremeness <- abs(permutation$scores)
  tie_tolerance <- sqrt(.Machine$double.eps) *
    max(1, observed_extremeness, permutation_extremeness)
  extreme <- permutation_extremeness >=
    observed_extremeness - tie_tolerance
  if (use_exact) {
    p_value <- mean(extreme)
    resolution <- 1 / nrow(permutation$assignments)
    permutation_method <- "exact"
  } else {
    p_value <- (1 + sum(extreme)) / (n_permutations + 1)
    resolution <- 1 / (n_permutations + 1)
    permutation_method <- "Monte Carlo"
  }
  n_ties <- sum(
    abs(permutation_extremeness - observed_extremeness) <= tie_tolerance
  )

  if (resolution > alpha) {
    warning(
      "The permutation design cannot attain `alpha = ", alpha,
      "`; its minimum p-value is ", signif(resolution, 4), ".",
      call. = FALSE
    )
  }

  edge_inference <- .brainnet_edge_inference(
    observed_counts,
    x$group_sizes,
    encoded$edge_indices,
    names(x$group_sizes),
    alpha,
    edge_method,
    adjust
  )

  ablation_result <- if (isTRUE(ablation)) {
    .brainnet_ablation(
      observed_counts = observed_counts,
      observed_score = observed_score,
      edge_order = edge_inference$order,
      assignments = permutation$assignments,
      edge_matrix = encoded$edge_matrix,
      group_sizes = x$group_sizes,
      exact = use_exact
    )
  } else {
    NULL
  }

  central_graphs <- lapply(x$populations, central_graph)
  global <- structure(
    list(
      statistic = c(T = unname(observed_score)),
      extremeness = c(`|T|` = unname(observed_extremeness)),
      p_value = unname(p_value),
      p.value = unname(p_value),
      alternative = "two-sided extremeness in the Fraiman score",
      method = "Whole-graph label-randomization test",
      data.name = deparse1(substitute(x)),
      permutation_method = permutation_method,
      n_assignments = if (use_exact) {
        nrow(permutation$assignments)
      } else {
        n_permutations
      },
      total_distinct_assignments = n_assignments,
      resolution = resolution,
      n_ties = as.integer(n_ties),
      tie_tolerance = unname(tie_tolerance)
    ),
    class = c("brainnet_global_test", "htest")
  )

  structure(
    list(
      global = global,
      edges = edge_inference$edges,
      edge_inference = edge_inference[
        c(
          "family_size", "method", "adjustment", "alpha",
          "proportion_columns"
        )
      ],
      central_graphs = central_graphs,
      ablation = ablation_result,
      group_sizes = x$group_sizes,
      n_nodes = x$n_nodes,
      node_labels = x$node_labels,
      communities = x$communities,
      parameters = list(
        alpha = alpha,
        n_permutations = n_permutations,
        exact = exact,
        max_exact = max_exact,
        edge_method = edge_method,
        adjust = adjust,
        ablation = ablation,
        seed = seed
      ),
      elapsed_seconds = unname(proc.time()[["elapsed"]] - started),
      call = match.call()
    ),
    class = "brainnet_result"
  )
}
