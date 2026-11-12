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
  balanced <- (n1 == n2)
  for (ct in 0:(n %/% 2L)) {
    x <- max(0L, ct - n2):min(n1, ct)
    d <- stats::dhyper(x, n1, n2, ct)
    row_vals <- vapply(seq_along(x), function(i)
      min(1, sum(d[d <= d[i] * (1 + 1e-7)])), numeric(1))
    ## Enforce exact row-swap (when n1 == n2) and 0 <-> 1 column-complement
    ## symmetries so isomorphic 2 x 2 tables get bit-identical p values.
    if (balanced || 2L * ct == n) {
      half <- seq_len((length(x) + 1L) %/% 2L)
      row_vals[length(x) - half + 1L] <- row_vals[half]
    }
    H[ct + 1L, x + 1L] <- row_vals
    if (2L * ct != n)
      H[n - ct + 1L, n1 - x + 1L] <- row_vals
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

## Canonicalise distinct m-group count vectors under exact row permutations
## (among equal-sized groups) and 0 <-> 1 edge complementation, so isomorphic
## tables are evaluated once and receive bit-identical marginal p values.
.canonical_count_matrix <- function(u_counts, group_sizes) {
  n_vec <- as.integer(group_sizes)
  size_blocks <- split(seq_along(n_vec), n_vec)
  sort_blocks <- function(v) {
    out <- integer(length(v))
    for (idx in size_blocks) {
      out[idx] <- if (length(idx) == 1L) v[idx] else sort.int(v[idx])
    }
    out
  }
  t(vapply(seq_len(nrow(u_counts)), function(r) {
    c_vec <- as.integer(u_counts[r, ])
    v1 <- sort_blocks(c_vec)
    v2 <- sort_blocks(n_vec - c_vec)
    diff_idx <- which(v1 != v2)[1L]
    if (!is.na(diff_idx) && v2[diff_idx] < v1[diff_idx]) v2 else v1
  }, integer(length(n_vec))))
}

## Marginal p values of every edge, for the observed counts (an E x m matrix)
## and for every permutation replicate (a list of m E x B count matrices, as
## returned by .permuted_edge_counts()). Returns list(obs = E-vector,
## perm = E x B matrix). For two groups and Fisher's test the p values are read
## from a precomputed table; otherwise each distinct canonical count vector is
## tested once.
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
  canon <- .canonical_count_matrix(all_counts[first, , drop = FALSE],
                                   group_sizes)
  ckey <- do.call(paste, c(as.data.frame(canon), sep = "-"))
  uckey <- unique(ckey)
  cfirst <- match(uckey, ckey)
  pcu <- vapply(cfirst, function(r)
    .edge_test_pvalue(canon[r, ], group_sizes, method), numeric(1))
  pu <- pcu[match(ckey, uckey)]
  p_all <- pu[match(key, ukey)]
  list(obs = p_all[seq_len(E)], perm = matrix(p_all[-seq_len(E)], E, B))
}

## Running minimum from the end: element s is min(v[s], ..., v[length(v)]).
.suffix_min <- function(v) rev(cummin(rev(v)))

## Combined p value from the B + 1 values (observed first) of T and U:
## Q_j = min(rank of T_j, rank of U_j), ties taking the largest rank within
## the numerical allowance, and p = max(Q_0, #{j : Q_j <= Q_0}) / (B + 1).
## Flooring the count at Q_0 preserves exact randomization calibration while
## guaranteeing p_combined >= min(p_global, p_local) even when non-transitive
## tolerance intervals bridge nearby values.
.combined_pvalue <- function(T_all, U_all, tolerance = 0) {
  if (length(T_all) != length(U_all))
    stop("T and U must contain the same observed and permuted assignments.")
  Q <- pmin(.permutation_ranks(T_all, tolerance), .permutation_ranks(U_all))
  max(Q[1L], sum(Q <= Q[1L])) / length(Q)
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
  ## Relative allowance matching rank_edges() (.PVALUE_REL_TOL = 64 * eps),
  ## with NO absolute floor: distinct tiny Fisher p-values such as 1e-100 and
  ## 2e-100 stay distinct.
  tolerance + .PVALUE_REL_TOL * abs(x)
}
