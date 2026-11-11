# R/combined-test.R
# Internal building blocks of the combined test, which evaluates the global
# statistic T together with the local statistic U (the smallest marginal p
# value of the remaining edges) in one permutation test. Shared by
# global_test() and identify_critical_links(). None of these are exported.

## Two-sided Fisher exact p values for every 2 x 2 table with group sizes n1
## and n2: entry [ct + 1, c1 + 1] is the p value for an edge present in c1
## graphs of the first group and ct graphs in total. Computed from the
## hypergeometric probabilities with the relative tolerance used by
## stats::fisher.test(), so that it agrees with fisher.test() to rounding.
.fisher_pvalue_table <- function(n1, n2) {
  n <- n1 + n2
  H <- matrix(NA_real_, n + 1L, n1 + 1L)
  for (ct in 0:n) {
    x <- max(0L, ct - n2):min(n1, ct)
    d <- stats::dhyper(x, n1, n2, ct)
    H[ct + 1L, x + 1L] <- vapply(seq_along(x), function(i)
      min(1, sum(d[d <= d[i] * (1 + 1e-7)])), numeric(1))
  }
  H
}

## p value of the marginal test of one edge, given its counts in each group
## (same tests as compute_edge_pvalues()).
.edge_test_pvalue <- function(counts, group_sizes, method) {
  tab <- cbind(counts, group_sizes - counts)
  p <- switch(method,
    fisher = stats::fisher.test(tab)$p.value,
    chi.squared = {
      if (any(rowSums(tab) == 0) || any(colSums(tab) == 0))
        stats::fisher.test(tab)$p.value
      else suppressWarnings(stats::chisq.test(tab, correct = TRUE)$p.value)
    },
    prop = {
      if (all(counts == 0) || all(counts == group_sizes)) 1
      else suppressWarnings(stats::prop.test(counts, group_sizes)$p.value)
    },
    stop("Invalid method specified."))
  min(max(p, 0), 1)
}

## Marginal p values of every edge, for the observed counts (an E x m matrix)
## and for every permutation replicate (a list of m E x B count matrices, as
## returned by .permuted_edge_counts()). Returns list(obs = E-vector,
## perm = E x B matrix). For two groups and Fisher's test the p values are read
## from a precomputed table; otherwise each distinct vector of counts is tested
## once.
.marginal_pvalues <- function(obs_counts, perm_counts, group_sizes, method) {
  m <- length(group_sizes)
  E <- nrow(obs_counts)
  B <- ncol(perm_counts[[1L]])
  if (m == 2L && method == "fisher") {
    H  <- .fisher_pvalue_table(group_sizes[1L], group_sizes[2L])
    ct_obs  <- obs_counts[, 1L] + obs_counts[, 2L]
    ct_perm <- perm_counts[[1L]] + perm_counts[[2L]]
    return(list(
      obs  = H[cbind(ct_obs + 1L, obs_counts[, 1L] + 1L)],
      perm = matrix(H[cbind(as.vector(ct_perm) + 1L,
                            as.vector(perm_counts[[1L]]) + 1L)], E, B)))
  }
  all_counts <- rbind(obs_counts,
                      do.call(cbind, lapply(perm_counts, as.vector)))
  key <- do.call(paste, c(as.data.frame(all_counts), sep = "-"))
  ukey <- unique(key)
  first <- match(ukey, key)
  pu <- vapply(first, function(r)
    .edge_test_pvalue(all_counts[r, ], group_sizes, method), numeric(1))
  p_all <- pu[match(key, ukey)]
  list(obs = p_all[seq_len(E)], perm = matrix(p_all[-seq_len(E)], E, B))
}

## Running minimum from the end: element s is min(v[s], ..., v[length(v)]).
.suffix_min <- function(v) rev(cummin(rev(v)))

## Combined p value from the B + 1 values (observed first) of T and U:
## Q_j = min(rank of T_j, rank of U_j) / (B + 1), ties taking the largest rank,
## and p = #{j : Q_j <= Q_0} / (B + 1).
.combined_pvalue <- function(T_all, U_all, tolerance = 0) {
  if (length(T_all) != length(U_all))
    stop("T and U must contain the same observed and permuted assignments.")
  Q <- pmin(.permutation_ranks(T_all, tolerance), .permutation_ranks(U_all))
  mean(Q <= Q[1L])
}

## p value of a statistic for which small values are evidence against H0,
## from its B + 1 values (observed first): #{j : x_j <= x_0} / (B + 1).
.plus_one_pvalue <- function(x_all, tolerance = 0) {
  margin <- .comparison_margin(x_all, tolerance)
  mean(x_all <= x_all[1L] + margin[1L])
}

## Extend each empirical left tail only by its numerical error allowance.
## With zero allowance these are the usual largest ranks for exact ties.
## Applying the SAME rule to all B+1 assignments is permutation equivariant;
## the final exact rank of the integer Q values therefore remains a valid
## randomization calibration. No observed-only rounding is used.
.permutation_ranks <- function(x, tolerance = 0) {
  margin <- .comparison_margin(x, tolerance)
  findInterval(x + margin, sort(x))
}

.comparison_margin <- function(x, tolerance = 0) {
  if (!length(x) || any(!is.finite(x)))
    stop("Permutation statistics must be non-empty and finite.")
  if (!(length(tolerance) %in% c(1L, length(x))) ||
      any(!is.finite(tolerance)) || any(tolerance < 0))
    stop("Numerical tolerances must be finite and non-negative.")
  ## A few ulps for ordinary evaluation/multiplication, with NO absolute
  ## floor: e.g. distinct Fisher p-values of 1e-100 and 2e-100 stay distinct.
  tolerance + 8 * .Machine$double.eps * abs(x)
}
