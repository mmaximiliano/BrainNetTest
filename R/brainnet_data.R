#' Construct Validated Brain-Network Population Data
#'
#' Creates the input object used by [brainnet_test()]. Each observation must be
#' a binary, symmetric adjacency matrix with a zero diagonal, and all
#' observations must use the same node order. Matrix dimnames are validated
#' when present; otherwise matrix positions are treated as implicit labels.
#'
#' @param populations A named list with at least two groups. Each group is a
#'   list containing at least two adjacency matrices.
#' @param node_labels Optional unique character labels for the nodes. When
#'   omitted, consistent matrix dimnames are used if present.
#' @param communities Optional vector of community memberships, one per node.
#'
#' @return An object of class `brainnet_data`.
#' @export
brainnet_data <- function(populations, node_labels = NULL,
                          communities = NULL) {
  if (!is.list(populations) || length(populations) < 2L) {
    stop("`populations` must be a list with at least 2 groups.",
      call. = FALSE
    )
  }

  group_names <- names(populations)
  if (is.null(group_names) || anyNA(group_names) ||
    any(group_names == "") || anyDuplicated(group_names)) {
    stop("`populations` must be named with unique, non-empty group names.",
      call. = FALSE
    )
  }

  group_sizes <- lengths(populations)
  if (any(!vapply(populations, is.list, logical(1L))) ||
    any(group_sizes < 2L)) {
    stop("Every group must be a list containing at least 2 graphs.",
      call. = FALSE
    )
  }

  first <- populations[[1L]][[1L]]
  .validate_graph_matrix(first, "`populations[[1]][[1]]`")
  n_nodes <- nrow(first)
  canonical_names <- .matrix_node_names(first, "The first graph")

  for (group_index in seq_along(populations)) {
    for (graph_index in seq_along(populations[[group_index]])) {
      label <- sprintf(
        "`populations[[%d]][[%d]]`",
        group_index,
        graph_index
      )
      graph <- populations[[group_index]][[graph_index]]
      .validate_graph_matrix(graph, label, n_nodes = n_nodes)
      graph_names <- .matrix_node_names(graph, label)

      names_match <- identical(graph_names, canonical_names)
      if (!names_match) {
        stop(
          "All adjacency matrices must use the same node ordering and ",
          "dimnames.",
          call. = FALSE
        )
      }
    }
  }

  if (is.null(node_labels)) {
    node_labels <- canonical_names
  } else {
    if (!is.character(node_labels) ||
      length(node_labels) != n_nodes ||
      anyNA(node_labels) ||
      any(node_labels == "") ||
      anyDuplicated(node_labels)) {
      stop(
        "`node_labels` must contain one unique, non-missing character ",
        "label per node.",
        call. = FALSE
      )
    }
    if (!is.null(canonical_names) &&
      !identical(node_labels, canonical_names)) {
      stop(
        "`node_labels` must preserve the node ordering in matrix dimnames.",
        call. = FALSE
      )
    }
  }

  if (!is.null(communities) &&
    (length(communities) != n_nodes || anyNA(communities))) {
    stop(
      "`communities` must contain one non-missing value per node.",
      call. = FALSE
    )
  }

  group_sizes <- as.integer(group_sizes)
  names(group_sizes) <- group_names

  structure(
    list(
      populations = populations,
      group_sizes = group_sizes,
      n_nodes = as.integer(n_nodes),
      node_labels = node_labels,
      communities = communities
    ),
    class = "brainnet_data"
  )
}

.validate_brainnet_data <- function(x) {
  if (!inherits(x, "brainnet_data")) {
    stop("`x` must be a `brainnet_data` object; use `brainnet_data()`.",
      call. = FALSE
    )
  }
  validated <- brainnet_data(
    x$populations,
    node_labels = x$node_labels,
    communities = x$communities
  )
  derived_fields_match <- identical(x$group_sizes, validated$group_sizes) &&
    identical(x$n_nodes, validated$n_nodes) &&
    identical(x$node_labels, validated$node_labels) &&
    identical(x$communities, validated$communities)
  if (!derived_fields_match) {
    stop(
      "The `brainnet_data` object is internally inconsistent; reconstruct ",
      "it with `brainnet_data()`.",
      call. = FALSE
    )
  }
  invisible(validated)
}

#' Print and Summarize Brain-Network Data
#'
#' @param x,object A `brainnet_data` object.
#' @param digits Number of digits used for printed density summaries.
#' @param ... Additional arguments reserved for methods.
#'
#' @return `print()` returns `x` invisibly. `summary()` returns a
#'   `summary_brainnet_data` object containing group sizes, node count, mean
#'   densities, and metadata indicators.
#' @rdname brainnet_data_methods
#' @export
print.brainnet_data <- function(x, ...) {
  .validate_brainnet_data(x)
  cat(
    "<brainnet_data>\n",
    length(x$group_sizes), " groups; ",
    x$n_nodes, " nodes; ",
    sum(x$group_sizes), " graphs\n",
    sep = ""
  )
  cat(
    paste0("  ", names(x$group_sizes), ": ", x$group_sizes),
    sep = "\n"
  )
  cat("\n")
  invisible(x)
}

#' @rdname brainnet_data_methods
#' @export
summary.brainnet_data <- function(object, ...) {
  .validate_brainnet_data(object)
  upper <- upper.tri(matrix(FALSE, object$n_nodes, object$n_nodes))
  mean_density <- vapply(
    object$populations,
    function(group) {
      mean(vapply(group, function(graph) mean(graph[upper]), numeric(1L)))
    },
    numeric(1L)
  )

  structure(
    list(
      group_sizes = object$group_sizes,
      n_nodes = object$n_nodes,
      mean_density = mean_density,
      has_node_labels = !is.null(object$node_labels),
      has_communities = !is.null(object$communities)
    ),
    class = "summary_brainnet_data"
  )
}

#' @rdname brainnet_data_methods
#' @export
print.summary_brainnet_data <- function(x, digits = 3L, ...) {
  cat(
    "Brain-network populations\n",
    "  Nodes: ", x$n_nodes, "\n",
    "  Groups: ", length(x$group_sizes), "\n",
    sep = ""
  )
  for (group in names(x$group_sizes)) {
    cat(
      "  ", group, ": n = ", x$group_sizes[[group]],
      ", mean density = ", format(round(x$mean_density[[group]], digits)),
      "\n",
      sep = ""
    )
  }
  invisible(x)
}
