# Relative floating-point tolerance for marginal p-value tie comparisons,
# shared by rank_edges() and the permutation calibration helpers.
.PVALUE_REL_TOL <- 64 * .Machine$double.eps
