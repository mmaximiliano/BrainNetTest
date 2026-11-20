# R/utils-checks.R
# Internal argument validation helpers, shared by the exported functions.
# They exist so that a mistake in a call is reported in terms of the argument
# that is wrong, instead of surfacing much later from inside the numerical
# code. None of them are exported.

#' Describe a value compactly for use in an error message
#'
#' @noRd
#' @keywords internal
.fmt_value <- function(x) {
  if (is.null(x)) return("NULL")
  if (!is.atomic(x)) return(paste0("an object of class \"", class(x)[1L], "\""))
  txt <- if (is.character(x)) encodeString(x, quote = "\"") else
    format(x, trim = TRUE)
  if (length(txt) == 0L)
    return(paste0("a zero-length ", class(x)[1L], " vector"))
  if (length(txt) > 3L) txt <- c(txt[1:3], "...")
  if (length(txt) == 1L) txt else
    paste0("c(", paste(txt, collapse = ", "), ")")
}

#' Validate a single count argument
#'
#' @noRd
#' @keywords internal
.check_count <- function(x, arg) {
  ok <- is.numeric(x) && length(x) == 1L && is.finite(x) && x >= 1 &&
    x <= .Machine$integer.max && x == trunc(x)
  if (!ok) {
    stop("`", arg, "` must be a single positive integer, not ",
         .fmt_value(x), ".", call. = FALSE)
  }
  as.integer(x)
}

#' Validate a vector of count arguments
#'
#' @noRd
#' @keywords internal
.check_counts <- function(x, arg) {
  ok <- is.numeric(x) && length(x) > 0L && all(is.finite(x)) && all(x >= 1) &&
    all(x <= .Machine$integer.max) && all(x == trunc(x))
  if (!ok) {
    stop("`", arg, "` must be a vector of positive integers, not ",
         .fmt_value(x), ".", call. = FALSE)
  }
  as.integer(x)
}

#' Validate a probability strictly inside the unit interval
#'
#' @noRd
#' @keywords internal
.check_proportion <- function(x, arg) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0 || x >= 1) {
    stop("`", arg, "` must be a single number strictly between 0 and 1, not ",
         .fmt_value(x), ".", call. = FALSE)
  }
  as.numeric(x)
}

#' Validate a single positive number
#'
#' @noRd
#' @keywords internal
.check_positive <- function(x, arg) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0) {
    stop("`", arg, "` must be a single positive number, not ", .fmt_value(x),
         ".", call. = FALSE)
  }
  as.numeric(x)
}

#' Validate the nested list structure used to represent graph populations
#'
#' Checks common dimensions and finite entries. Binary inference additionally
#' requires undirected graphs without self-loops; general helpers can still
#' accept finite weighted matrices.
#'
#' @param populations Candidate populations object.
#' @param binary If `TRUE`, require binary, symmetric, zero-diagonal graphs.
#' @param arg Argument name to use in error messages.
#' @return Invisibly, the common number of nodes.
#' @noRd
#' @keywords internal
.check_populations <- function(populations, binary = FALSE,
                               arg = "populations") {
  if (!is.list(populations) || length(populations) == 0L) {
    stop("`", arg, "` must be a non-empty list of populations, not ",
         .fmt_value(populations), ".", call. = FALSE)
  }

  pop_names <- names(populations)
  n_nodes   <- NA_integer_

  for (k in seq_along(populations)) {
    named <- !is.null(pop_names) && !is.na(pop_names[k]) &&
      nzchar(pop_names[k])
    label <- if (named) paste0("\"", pop_names[k], "\"") else k
    pop   <- populations[[k]]

    if (!is.list(pop) || length(pop) == 0L) {
      stop("Population ", label, " of `", arg,
           "` must be a non-empty list of adjacency matrices, not ",
           .fmt_value(pop), ".", call. = FALSE)
    }

    for (g in seq_along(pop)) {
      G     <- pop[[g]]
      where <- paste0("Graph ", g, " of population ", label)

      if (!is.matrix(G) || !(is.numeric(G) || is.logical(G))) {
        stop(where, " is not a numeric matrix.", call. = FALSE)
      }
      if (nrow(G) != ncol(G)) {
        stop(where, " is not square (", nrow(G), " by ", ncol(G), ").",
             call. = FALSE)
      }
      if (nrow(G) == 0L) {
        stop(where, " must have at least one node.", call. = FALSE)
      }
      if (is.na(n_nodes)) {
        n_nodes <- nrow(G)
      } else if (nrow(G) != n_nodes) {
        stop(where, " has ", nrow(G), " nodes but an earlier graph has ",
             n_nodes, "; all graphs must share the same vertex set.",
             call. = FALSE)
      }
      if (binary && (any(!is.finite(G)) || any(G != 0 & G != 1))) {
        stop(where, " is not binary; adjacency matrices must contain only ",
             "0 and 1.", call. = FALSE)
      }
      if (any(!is.finite(G))) {
        stop(where, " must contain only finite entries.", call. = FALSE)
      }
      if (binary && any(G != t(G))) {
        stop(where, " must be symmetric.", call. = FALSE)
      }
      if (binary && any(diag(G) != 0)) {
        stop(where, " must have a zero diagonal (no self-loops).", call. = FALSE)
      }
    }
  }

  invisible(n_nodes)
}

# Validate the shared input contract of the permutation inference functions.
.check_inference_populations <- function(populations) {
  if (!is.list(populations) || length(populations) < 2L)
    stop("`populations` must be a list with at least 2 groups.", call. = FALSE)
  n_nodes <- .check_populations(populations, binary = TRUE)
  if (any(lengths(populations) < 2L))
    stop("Each population must have at least two graphs.", call. = FALSE)
  n_nodes
}
