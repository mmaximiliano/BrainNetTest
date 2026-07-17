#' Compute a Population Central Graph
#'
#' Computes the element-wise mean of aligned binary adjacency matrices.
#'
#' @param graphs A non-empty list of aligned binary, symmetric adjacency
#'   matrices.
#'
#' @return A symmetric numeric matrix whose entries are empirical edge
#'   proportions.
#' @export
central_graph <- function(graphs) {
  if (!is.list(graphs) || length(graphs) == 0L) {
    stop("`graphs` must be a non-empty list of adjacency matrices.",
      call. = FALSE
    )
  }

  first <- graphs[[1L]]
  .validate_graph_matrix(first, "`graphs[[1]]`")
  n_nodes <- nrow(first)
  canonical_names <- .matrix_node_names(first, "The first graph")

  for (index in seq_along(graphs)) {
    label <- sprintf("`graphs[[%d]]`", index)
    .validate_graph_matrix(graphs[[index]], label, n_nodes = n_nodes)
    if (!identical(
      .matrix_node_names(graphs[[index]], label),
      canonical_names
    )) {
      stop("All graphs must use the same node ordering and dimnames.",
        call. = FALSE
      )
    }
  }

  Reduce("+", graphs) / length(graphs)
}

.validate_weighted_graph_matrix <- function(graph, label) {
  .validate_graph_matrix(graph, label, binary = FALSE)
  if (any(graph < 0 | graph > 1)) {
    stop(label, " must contain values in [0, 1].", call. = FALSE)
  }
  invisible(graph)
}

#' Compute L1 Distance Between Undirected Graphs
#'
#' Computes the Manhattan distance over the upper triangle, so each undirected
#' edge is counted exactly once.
#'
#' @param x,y Aligned symmetric adjacency or central-graph matrices with values
#'   in `[0, 1]` and zero diagonals.
#'
#' @return One non-negative numeric distance.
#' @export
graph_distance <- function(x, y) {
  .validate_weighted_graph_matrix(x, "`x`")
  .validate_weighted_graph_matrix(y, "`y`")
  if (!identical(dim(x), dim(y))) {
    stop("`x` and `y` must have the same dimensions.", call. = FALSE)
  }

  x_names <- .matrix_node_names(x, "`x`")
  y_names <- .matrix_node_names(y, "`y`")
  if (!identical(x_names, y_names)) {
    stop("`x` and `y` must use the same node ordering and dimnames.",
      call. = FALSE
    )
  }

  sum(abs(x[upper.tri(x)] - y[upper.tri(y)]))
}
