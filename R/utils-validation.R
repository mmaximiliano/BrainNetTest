.assert_scalar_number <- function(x, arg, lower = -Inf, upper = Inf,
                                  lower_open = FALSE, upper_open = FALSE) {
  valid <- is.numeric(x) && length(x) == 1L && !is.na(x) && is.finite(x)
  if (valid) {
    valid <- if (lower_open) x > lower else x >= lower
  }
  if (valid) {
    valid <- if (upper_open) x < upper else x <= upper
  }
  if (!valid) {
    interval <- paste0(
      if (lower_open) "(" else "[",
      lower,
      ", ",
      upper,
      if (upper_open) ")" else "]"
    )
    stop("`", arg, "` must be one finite number in ", interval, ".",
      call. = FALSE
    )
  }
  invisible(x)
}

.assert_whole_number <- function(x, arg, minimum = 1L) {
  valid <- is.numeric(x) && length(x) == 1L && !is.na(x) &&
    is.finite(x) && x >= minimum && x == floor(x) &&
    x <= .Machine$integer.max
  if (!valid) {
    stop(
      "`", arg, "` must be a whole number from ", minimum,
      " through ", .Machine$integer.max, ".",
      call. = FALSE
    )
  }
  invisible(as.integer(x))
}

.validate_graph_matrix <- function(graph, label, n_nodes = NULL,
                                   binary = TRUE) {
  if (!is.matrix(graph) ||
    !(is.numeric(graph) || is.integer(graph) || is.logical(graph))) {
    stop(label, " must be a numeric, integer, or logical matrix.",
      call. = FALSE
    )
  }
  if (nrow(graph) != ncol(graph)) {
    stop(label, " must be square.", call. = FALSE)
  }
  if (!is.null(n_nodes) && !identical(dim(graph), c(n_nodes, n_nodes))) {
    stop("All adjacency matrices must have the same dimensions.",
      call. = FALSE
    )
  }
  if (nrow(graph) < 2L) {
    stop("Adjacency matrices must contain at least 2 nodes.",
      call. = FALSE
    )
  }
  if (anyNA(graph) || any(!is.finite(graph))) {
    stop(label, " must contain only finite, non-missing values.",
      call. = FALSE
    )
  }
  if (binary && any(!(graph %in% c(0, 1)))) {
    stop(label, " must be binary with values in {0, 1}.", call. = FALSE)
  }
  if (!isTRUE(isSymmetric(graph))) {
    stop(label, " must be symmetric.", call. = FALSE)
  }
  if (any(diag(graph) != 0)) {
    stop(label, " must have a zero diagonal.", call. = FALSE)
  }
  invisible(graph)
}

.matrix_node_names <- function(graph, label) {
  row_names <- rownames(graph)
  col_names <- colnames(graph)

  if (is.null(row_names) && is.null(col_names)) {
    return(NULL)
  }
  if (is.null(row_names) || is.null(col_names) ||
    !identical(row_names, col_names)) {
    stop(label, " must use identical row and column node names.",
      call. = FALSE
    )
  }
  if (anyNA(row_names) || any(row_names == "") || anyDuplicated(row_names)) {
    stop(label, " must use unique, non-missing node names.",
      call. = FALSE
    )
  }
  row_names
}

.normalize_community_sizes <- function(n_nodes, n_communities,
                                       community_sizes) {
  if (is.null(community_sizes)) {
    base_size <- n_nodes %/% n_communities
    remainder <- n_nodes %% n_communities
    sizes <- rep.int(base_size, n_communities)
    if (remainder > 0L) {
      sizes[seq_len(remainder)] <- sizes[seq_len(remainder)] + 1L
    }
    return(as.integer(sizes))
  }

  valid <- is.numeric(community_sizes) &&
    length(community_sizes) == n_communities &&
    all(is.finite(community_sizes)) &&
    all(community_sizes >= 1) &&
    all(community_sizes == floor(community_sizes))
  if (!valid) {
    stop(
      "`community_sizes` must contain one positive whole number per ",
      "community.",
      call. = FALSE
    )
  }
  if (sum(community_sizes) != n_nodes) {
    stop(
      "The sum of `community_sizes` must equal `n_nodes`.",
      call. = FALSE
    )
  }
  as.integer(community_sizes)
}

.validate_probability_vector <- function(x, arg, length_out) {
  valid <- is.numeric(x) &&
    all(is.finite(x)) &&
    length(x) %in% c(1L, length_out) &&
    all(x >= 0 & x <= 1)
  if (!valid) {
    stop(
      "`", arg, "` must be one probability or ", length_out,
      " probabilities in [0, 1].",
      call. = FALSE
    )
  }
  rep(as.numeric(x), length.out = length_out)
}

.with_preserved_seed <- function(seed, code) {
  if (is.null(seed)) {
    return(code())
  }
  .assert_whole_number(seed, "seed", minimum = 0L)

  global <- globalenv()
  had_seed <- exists(".Random.seed", envir = global, inherits = FALSE)
  if (had_seed) {
    old_seed <- get(".Random.seed", envir = global, inherits = FALSE)
  }

  on.exit(
    {
      if (had_seed) {
        assign(".Random.seed", old_seed, envir = global)
      } else if (exists(".Random.seed", envir = global, inherits = FALSE)) {
        rm(".Random.seed", envir = global)
      }
    },
    add = TRUE
  )

  set.seed(as.integer(seed))
  code()
}
