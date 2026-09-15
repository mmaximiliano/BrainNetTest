# R/critical_links.R
# The "critical_links" class returned by identify_critical_links(), and its
# print, summary and plot methods.

#' Construct a critical_links object
#'
#' Both the optimised and the naive implementations return through here, so
#' that the two are interchangeable in every respect.
#'
#' @noRd
#' @keywords internal
.critical_links <- function(critical_edges, edges_removed,
                            modified_populations, p_value, n_edges,
                            settings, call) {
  structure(
    list(critical_edges       = critical_edges,
         edges_removed        = edges_removed,
         modified_populations = modified_populations,
         p_value              = p_value,
         n_edges              = n_edges,
         settings             = settings,
         call                 = call),
    class = "critical_links")
}

#' Format a permutation p value for display
#'
#' The estimate is a proportion out of `n_permutations` draws, so it cannot
#' resolve below 1 / n_permutations; report that bound rather than a
#' misleading zero. The bound is rounded to two significant digits so that,
#' say, 3000 replicates read "< 0.00033" rather than "< 0.0003333333".
#'
#' @noRd
#' @keywords internal
.format_pvalue <- function(p, n_permutations) {
  if (is.na(p)) return("NA")
  eps <- 1 / n_permutations
  if (p < eps) paste0("< ", format(signif(eps, 2), scientific = FALSE))
  else format(round(p, 4), nsmall = 4, scientific = FALSE)
}

#' Number of Critical Edges Found
#'
#' @param x A \code{critical_links} object.
#' @return Integer count, zero when no critical set was identified.
#' @noRd
#' @keywords internal
.n_critical <- function(x) {
  if (is.null(x$critical_edges)) 0L else nrow(x$critical_edges)
}

#' Print a Critical-Edge Analysis
#'
#' Compact report of an \code{\link{identify_critical_links}} analysis: the
#' populations compared, the outcome of the global test, and the size of the
#' critical set together with its most significant edges. Use
#' \code{\link{summary.critical_links}} for the call, the analysis settings
#' and the node-level view as well.
#'
#' @param x An object of class \code{"critical_links"}, as returned by
#'   \code{\link{identify_critical_links}}.
#' @param n Number of critical edges to display. Default \code{6}.
#' @param digits Number of significant digits used for the $p$ values.
#'   Default \code{getOption("digits")}.
#' @param ... Further arguments, currently ignored.
#'
#' @return \code{x}, invisibly.
#' @seealso \code{\link{summary.critical_links}},
#'   \code{\link{plot.critical_links}}.
#' @export
#' @examples
#' \donttest{
#' control <- generate_category_graphs(n_graphs = 15, n_nodes = 10,
#'   n_communities = 2, base_intra_prob = 0.8, base_inter_prob = 0.2, seed = 1)
#' patient <- generate_category_graphs(n_graphs = 15, n_nodes = 10,
#'   n_communities = 2, base_intra_prob = 0.5, base_inter_prob = 0.5, seed = 2)
#' result <- identify_critical_links(list(Control = control, Patient = patient),
#'   n_permutations = 200, seed = 42)
#' result
#' }
print.critical_links <- function(x, n = 6, digits = getOption("digits"), ...) {
  s <- x$settings
  cat("\nCritical edges between populations of brain networks\n\n")

  groups <- names(x$modified_populations)
  if (is.null(groups)) groups <- paste("group", seq_along(x$modified_populations))
  cat(sprintf("Populations:     %d (%s), %s graphs\n",
              length(groups), paste(groups, collapse = ", "),
              paste(lengths(x$modified_populations), collapse = " / ")))
  cat(sprintf("Candidate edges: %d\n", x$n_edges))

  if (is.na(x$p_value)) {
    cat("Global test:     not evaluated (no candidate edges)\n")
  } else {
    cat(sprintf("Global test:     p %s (%d permutations, alpha = %s)\n",
                .format_pvalue(x$p_value, s$n_permutations),
                s$n_permutations, format(s$alpha)))
  }

  k <- .n_critical(x)
  if (k == 0L) {
    cat("Critical edges:  none; the populations are not distinguishable",
        "at this level\n\n")
    return(invisible(x))
  }

  cat(sprintf("Critical edges:  %d of %d (%.1f%%), ranked by %s\n\n",
              k, x$n_edges, 100 * k / x$n_edges, .method_label(s$method)))

  show <- min(n, k)
  cat(sprintf("Most significant %s:\n", if (show < k)
    sprintf("%d of them", show) else "edges"))
  print(format(x$critical_edges[seq_len(show), ], digits = digits))
  if (show < k)
    cat(sprintf("... %d more; see the critical_edges component\n", k - show))
  cat("\n")

  invisible(x)
}

#' @noRd
#' @keywords internal
.method_label <- function(method) {
  switch(method,
         fisher        = "Fisher's exact test",
         chi.squared   = "a chi-squared test",
         prop          = "a test of equal proportions",
         method)
}

#' Summarise a Critical-Edge Analysis
#'
#' Extends \code{\link{print.critical_links}} with the settings the analysis
#' used and a node-level view of the critical set: the nodes that participate
#' in the most critical edges, as returned by \code{\link{get_critical_nodes}}.
#'
#' @param object An object of class \code{"critical_links"}.
#' @param node_labels Optional character vector of node labels, passed to
#'   \code{\link{get_critical_nodes}}.
#' @param ... Further arguments, currently ignored.
#'
#' @return An object of class \code{"summary.critical_links"}, a list with the
#'   components \code{n_critical}, \code{n_edges}, \code{p_value},
#'   \code{settings}, \code{group_sizes}, \code{critical_edges} and
#'   \code{critical_nodes}.
#' @seealso \code{\link{identify_critical_links}},
#'   \code{\link{get_critical_nodes}}.
#' @export
#' @examples
#' \donttest{
#' control <- generate_category_graphs(n_graphs = 15, n_nodes = 10,
#'   n_communities = 2, base_intra_prob = 0.8, base_inter_prob = 0.2, seed = 1)
#' patient <- generate_category_graphs(n_graphs = 15, n_nodes = 10,
#'   n_communities = 2, base_intra_prob = 0.5, base_inter_prob = 0.5, seed = 2)
#' result <- identify_critical_links(list(Control = control, Patient = patient),
#'   n_permutations = 200, seed = 42)
#' summary(result)
#' }
summary.critical_links <- function(object, node_labels = NULL, ...) {
  structure(
    list(call           = object$call,
         n_critical     = .n_critical(object),
         n_edges        = object$n_edges,
         p_value        = object$p_value,
         settings       = object$settings,
         group_sizes    = lengths(object$modified_populations),
         critical_edges = object$critical_edges,
         critical_nodes = get_critical_nodes(object, node_labels)),
    class = "summary.critical_links")
}

#' Print a Critical-Edge Analysis Summary
#'
#' @param x An object of class \code{"summary.critical_links"}.
#' @param n Number of critical edges and nodes to display. Default \code{10}.
#' @param digits Number of significant digits used for the $p$ values.
#'   Default \code{getOption("digits")}.
#' @param ... Further arguments, currently ignored.
#'
#' @return \code{x}, invisibly.
#' @export
print.summary.critical_links <- function(x, n = 10,
                                         digits = getOption("digits"), ...) {
  s <- x$settings
  cat("\nCritical edges between populations of brain networks\n\n")
  if (!is.null(x$call)) {
    cat("Call:\n")
    print(x$call)
    cat("\n")
  }

  cat(sprintf("Populations:     %d (%s), %s graphs\n",
              length(x$group_sizes), paste(names(x$group_sizes),
                                           collapse = ", "),
              paste(x$group_sizes, collapse = " / ")))
  cat(sprintf("Candidate edges: %d\n", x$n_edges))
  if (!is.na(x$p_value))
    cat(sprintf("Global test:     p %s\n",
                .format_pvalue(x$p_value, s$n_permutations)))
  cat(sprintf("Critical edges:  %d (%.1f%% of candidates)\n",
              x$n_critical, 100 * x$n_critical / max(x$n_edges, 1)))

  cat("\nSettings:\n")
  cat(sprintf("  marginal test        %s\n", .method_label(s$method)))
  cat(sprintf("  multiplicity         %s\n", s$adjust_method))
  cat(sprintf("  significance level   %s\n", format(s$alpha)))
  cat(sprintf("  edges removed / step %d\n", s$batch_size))
  cat(sprintf("  permutations         %d\n", s$n_permutations))
  cat(sprintf("  normalisation a      %s\n", format(s$a)))
  cat(sprintf("  seed                 %s\n",
              if (is.null(s$seed)) "not set" else format(s$seed)))

  if (x$n_critical == 0L) {
    cat("\nNo critical edges were identified.\n\n")
    return(invisible(x))
  }

  cat("\nMost significant edges:\n")
  print(format(utils::head(x$critical_edges, n), digits = digits))
  if (x$n_critical > n)
    cat(sprintf("... %d more\n", x$n_critical - n))

  cat("\nNodes ranked by critical degree:\n")
  print(utils::head(x$critical_nodes, n))
  if (nrow(x$critical_nodes) > n)
    cat(sprintf("... %d more\n", nrow(x$critical_nodes) - n))
  cat("\n")

  invisible(x)
}

#' Plot a Critical-Edge Analysis
#'
#' Method interface to \code{\link{plot_critical_edges}}: draws the weighted
#' central graph of each population alongside a panel highlighting the
#' critical edges.
#'
#' The populations must be supplied, because the object stores the
#' \emph{modified} populations, in which the critical edges have been removed,
#' and plotting them would hide the very structure the figure is about.
#'
#' @param x An object of class \code{"critical_links"}.
#' @param populations The named \code{list} of populations that was passed to
#'   \code{\link{identify_critical_links}}.
#' @param ... Further arguments passed to \code{\link{plot_critical_edges}},
#'   such as \code{communities}, \code{layout} or \code{reference}.
#'
#' @return Invisibly \code{NULL}; called for its side effect.
#' @seealso \code{\link{plot_critical_edges}}, which this method wraps and
#'   which can also be called directly.
#' @export
#' @examples
#' \donttest{
#' set.seed(1)
#' community_sizes <- c(4, 3, 3)
#' control <- generate_category_graphs(n_graphs = 30, n_nodes = 10,
#'   n_communities = 3, community_sizes = community_sizes,
#'   base_intra_prob = rep(0.8, 3), base_inter_prob = 0.1, seed = 1)
#' patient <- generate_category_graphs(n_graphs = 30, n_nodes = 10,
#'   n_communities = 3, community_sizes = community_sizes,
#'   base_intra_prob = c(0.3, 0.8, 0.8), base_inter_prob = 0.1, seed = 2)
#' populations <- list(Control = control, Patient = patient)
#' result <- identify_critical_links(populations, n_permutations = 200,
#'   seed = 42)
#' plot(result, populations,
#'   communities = rep(seq_along(community_sizes), times = community_sizes))
#' }
plot.critical_links <- function(x, populations, ...) {
  if (missing(populations))
    stop("`populations` is required: supply the populations that were passed ",
         "to identify_critical_links(). The object stores the modified ",
         "populations, from which the critical edges have been removed.",
         call. = FALSE)
  plot_critical_edges(populations, x, ...)
}
