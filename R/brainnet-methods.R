.validate_brainnet_result <- function(x) {
  required <- c(
    "global", "edges", "edge_inference", "central_graphs",
    "group_sizes", "n_nodes"
  )
  if (!inherits(x, "brainnet_result") ||
    !all(required %in% names(x))) {
    stop("`x` must be a `brainnet_result` from `brainnet_test()`.",
      call. = FALSE
    )
  }
  edge_columns <- c(
    "node1", "node2", "p_value", "p_adjusted", "selected", "rank"
  )
  expected_edges <- as.integer(x$n_nodes * (x$n_nodes - 1L) / 2L)
  valid <- is.list(x$global) &&
    is.numeric(x$global$statistic) &&
    length(x$global$statistic) == 1L &&
    identical(names(x$global$statistic), "T") &&
    is.numeric(x$global$p_value) &&
    length(x$global$p_value) == 1L &&
    is.finite(x$global$p_value) &&
    x$global$p_value >= 0 &&
    x$global$p_value <= 1 &&
    is.data.frame(x$edges) &&
    all(edge_columns %in% names(x$edges)) &&
    nrow(x$edges) == expected_edges &&
    identical(x$edge_inference$family_size, expected_edges) &&
    length(x$central_graphs) == length(x$group_sizes) &&
    all(vapply(
      x$central_graphs,
      function(graph) {
        is.matrix(graph) &&
          identical(dim(graph), c(x$n_nodes, x$n_nodes))
      },
      logical(1L)
    ))
  if (!valid) {
    stop(
      "The `brainnet_result` object is internally inconsistent.",
      call. = FALSE
    )
  }
  invisible(x)
}

#' Print and Summarize a BrainNetTest Result
#'
#' @param x,object A `brainnet_result`.
#' @param digits Number of significant digits used for printed values.
#' @param top Number of highest-ranked edges retained by `summary()`.
#' @param ... Additional arguments reserved for methods.
#'
#' @return `print()` returns `x` invisibly. `summary()` returns a
#'   `summary_brainnet_result` containing the global result, group sizes,
#'   edge-inference settings, selected-edge count, and highest-ranked edges.
#' @rdname brainnet_result_methods
#' @export
print.brainnet_result <- function(x, digits = 4L, ...) {
  .validate_brainnet_result(x)
  cat("<brainnet_result>\n")
  cat(
    "Global permutation test: T = ",
    format(signif(x$global$statistic, digits)),
    ", p = ",
    format.pval(x$global$p_value, digits = digits),
    " (", x$global$permutation_method, ")\n",
    sep = ""
  )
  adjustment_label <- if (x$edge_inference$adjustment == "holm") {
    "Holm"
  } else {
    x$edge_inference$adjustment
  }
  cat(
    "Edge inference: ", x$edge_inference$method,
    " + ", adjustment_label,
    "; ", sum(x$edges$selected), " of ",
    x$edge_inference$family_size, " edges selected\n",
    sep = ""
  )
  invisible(x)
}

#' @rdname brainnet_result_methods
#' @export
summary.brainnet_result <- function(object, top = 6L, ...) {
  .validate_brainnet_result(object)
  top <- .assert_whole_number(top, "top", minimum = 1L)
  selected <- selected_edges(object)

  structure(
    list(
      call = object$call,
      global = object$global,
      group_sizes = object$group_sizes,
      edge_inference = object$edge_inference,
      n_selected = nrow(selected),
      top_edges = utils::head(
        object$edges[order(object$edges$rank), , drop = FALSE],
        top
      ),
      elapsed_seconds = object$elapsed_seconds
    ),
    class = "summary_brainnet_result"
  )
}

#' @rdname brainnet_result_methods
#' @export
print.summary_brainnet_result <- function(x, digits = 4L, ...) {
  cat("BrainNetTest result summary\n")
  cat(
    "  Groups: ",
    paste0(names(x$group_sizes), " (n=", x$group_sizes, ")",
      collapse = ", "
    ),
    "\n",
    sep = ""
  )
  cat(
    "  Global T: ", format(signif(x$global$statistic, digits)),
    "\n  Global p: ",
    format.pval(x$global$p_value, digits = digits),
    "\n  Selected edges: ", x$n_selected,
    " (", x$edge_inference$adjustment, " adjustment)\n",
    sep = ""
  )
  if (nrow(x$top_edges) > 0L) {
    cat("  Highest-ranked edges:\n")
    print(x$top_edges, row.names = FALSE)
  }
  invisible(x)
}

#' Convert a BrainNetTest Result to an Edge Table
#'
#' @param x A `brainnet_result`.
#' @param row.names,optional,... Arguments required by the base
#'   [as.data.frame()] generic.
#'
#' @return The complete upper-triangle edge-inference table.
#' @export
as.data.frame.brainnet_result <- function(x, row.names = NULL,
                                          optional = FALSE, ...) {
  .validate_brainnet_result(x)
  x$edges
}

#' Extract Multiplicity-Controlled Selected Edges
#'
#' @param x A `brainnet_result`.
#'
#' @return A data frame containing only edges selected by the configured
#'   adjusted p-value rule, ordered deterministically by rank.
#' @export
selected_edges <- function(x) {
  .validate_brainnet_result(x)
  selected <- x$edges[x$edges$selected, , drop = FALSE]
  selected <- selected[order(selected$rank), , drop = FALSE]
  rownames(selected) <- NULL
  selected
}

#' Summarize Nodes Incident to Selected Edges
#'
#' This is a descriptive incidence summary. It does not provide node-level
#' error control.
#'
#' @param x A `brainnet_result`.
#'
#' @return A data frame with node index, optional label, and selected degree.
#' @export
selected_nodes <- function(x) {
  .validate_brainnet_result(x)
  edges <- selected_edges(x)
  if (nrow(edges) == 0L) {
    output <- data.frame(
      node = integer(),
      selected_degree = integer()
    )
    if (!is.null(x$node_labels)) {
      output$label <- character()
      output <- output[, c("node", "label", "selected_degree")]
    }
    return(output)
  }

  incidence <- table(c(edges$node1, edges$node2))
  output <- data.frame(
    node = as.integer(names(incidence)),
    selected_degree = as.integer(incidence)
  )
  output <- output[
    order(-output$selected_degree, output$node), ,
    drop = FALSE
  ]
  rownames(output) <- NULL

  if (!is.null(x$node_labels)) {
    output$label <- x$node_labels[output$node]
    output <- output[, c("node", "label", "selected_degree")]
  }
  output
}

#' Plot a BrainNetTest Result
#'
#' @param x A `brainnet_result`.
#' @param reference Group name or index used as the selected-edge background.
#' @param threshold Threshold for showing background central-graph edges.
#' @param layout Optional igraph layout function or coordinate matrix.
#' @param ... Additional arguments passed to [igraph::plot.igraph()].
#'
#' @return `x`, invisibly.
#' @export
plot.brainnet_result <- function(x, reference = 1L, threshold = 0.5,
                                 layout = NULL, ...) {
  .validate_brainnet_result(x)
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package 'igraph' is required for plotting.", call. = FALSE)
  }
  .assert_scalar_number(threshold, "threshold", lower = 0, upper = 1)

  group_names <- names(x$central_graphs)
  if (is.character(reference)) {
    reference <- match(reference, group_names)
  }
  reference <- .assert_whole_number(reference, "reference")
  if (reference > length(x$central_graphs)) {
    stop("`reference` is outside the available group range.",
      call. = FALSE
    )
  }
  if (is.null(layout)) {
    layout <- igraph::layout_in_circle
  }

  n_panels <- length(x$central_graphs) + 1L
  n_columns <- ceiling(sqrt(n_panels))
  n_rows <- ceiling(n_panels / n_columns)
  old_par <- graphics::par(
    mfrow = c(n_rows, n_columns),
    mar = c(1, 1, 2, 1)
  )
  on.exit(graphics::par(old_par), add = TRUE)

  vertex_color <- if (is.null(x$communities)) {
    "skyblue"
  } else {
    as.integer(factor(x$communities))
  }
  vertex_label <- if (is.null(x$node_labels)) {
    seq_len(x$n_nodes)
  } else {
    x$node_labels
  }

  for (group_index in seq_along(x$central_graphs)) {
    graph <- igraph::graph_from_adjacency_matrix(
      x$central_graphs[[group_index]],
      mode = "undirected",
      weighted = TRUE,
      diag = FALSE
    )
    weights <- igraph::E(graph)$weight
    igraph::plot.igraph(
      graph,
      layout = layout,
      main = paste(group_names[[group_index]], "central graph"),
      vertex.color = vertex_color,
      vertex.label = vertex_label,
      edge.width = 1 + 4 * weights,
      ...
    )
  }

  background <- x$central_graphs[[reference]] > threshold
  diag(background) <- FALSE
  selected <- selected_edges(x)
  if (nrow(selected) > 0L) {
    background[cbind(selected$node1, selected$node2)] <- TRUE
    background[cbind(selected$node2, selected$node1)] <- TRUE
  }
  graph <- igraph::graph_from_adjacency_matrix(
    background,
    mode = "undirected",
    diag = FALSE
  )
  edge_list <- igraph::as_edgelist(graph, names = FALSE)
  is_selected <- rep(FALSE, nrow(edge_list))
  if (nrow(selected) > 0L && nrow(edge_list) > 0L) {
    selected_keys <- paste(
      pmin(selected$node1, selected$node2),
      pmax(selected$node1, selected$node2),
      sep = "-"
    )
    edge_keys <- paste(
      pmin(edge_list[, 1L], edge_list[, 2L]),
      pmax(edge_list[, 1L], edge_list[, 2L]),
      sep = "-"
    )
    is_selected <- edge_keys %in% selected_keys
  }
  igraph::plot.igraph(
    graph,
    layout = layout,
    main = "Multiplicity-controlled selected edges",
    vertex.color = vertex_color,
    vertex.label = vertex_label,
    edge.color = ifelse(is_selected, "red", "grey80"),
    edge.width = ifelse(is_selected, 3, 1),
    ...
  )

  invisible(x)
}
