# Validate marginal-test selection before an edge loop or early return.
.check_edge_method <- function(method, n_populations) {
  choices <- c("fisher", "chi.squared", "prop")
  if (!is.character(method) || length(method) != 1L || is.na(method) ||
      !method %in% choices) {
    stop("Invalid method specified. `method` must be one of ",
         '"fisher", "chi.squared", or "prop".', call. = FALSE)
  }
  if (method == "prop" && n_populations != 2L)
    stop("prop.test is only implemented for two populations.", call. = FALSE)
  method
}
