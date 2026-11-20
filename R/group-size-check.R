# The current lower-tail T test is supported only for balanced designs.
# Keep this separate from graph validation: descriptive edge summaries and
# marginal edge tests can still be used with unequal group sizes.
.check_equal_group_sizes <- function(group_sizes) {
  if (length(unique(group_sizes)) != 1L)
    stop("The current T test requires equal group sizes. Unequal group sizes are not supported; a revised test is needed.",
         call. = FALSE)
  invisible(group_sizes)
}
