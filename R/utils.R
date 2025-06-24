#' Null coalescing operator
#'
#' @param x Left operand
#' @param y Right operand
#' @return x if not NULL, otherwise y
#' @noRd
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

#' Check if package is available
#'
#' @param package Name of the package
#' @param quietly Whether to suppress messages
#' @return TRUE if the package is available, FALSE otherwise
#' @noRd
is_package_available <- function(package, quietly = TRUE) {
  if (quietly) {
    suppressWarnings(suppressMessages(
      requireNamespace(package, quietly = TRUE)
    ))
  } else {
    requireNamespace(package, quietly = TRUE)
  }
}

#' Get number of cores for parallel processing
#'
#' @param max_cores Maximum number of cores to use
#' @return Number of cores to use for parallel processing
#' @noRd
get_usable_cores <- function(max_cores = NULL) {
  available_cores <- parallel::detectCores() - 1
  if (is.null(max_cores)) {
    return(max(1, available_cores))
  } else {
    return(min(max_cores, available_cores))
  }
}
