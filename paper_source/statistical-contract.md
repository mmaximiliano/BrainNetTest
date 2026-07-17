# BrainNetTest 1.0 statistical contract

This document freezes the inferential target used by the package and paper.
Any change to this contract requires corresponding changes to the exact tests,
simulation design, documentation, and manuscript.

## Supported observations

An observation is one participant's binary, undirected brain network without
self-loops, represented by a symmetric square adjacency matrix with a zero
diagonal. All observations use the same labeled and ordered brain atlas.
Groups contain independent subject-level observations and have fixed sample
sizes \(n_1,\ldots,n_m\), with \(m \geq 2\) and \(n_k \geq 2\).

Weighted or directed graphs, missing nodes, non-finite entries, paired or
clustered observations, covariate adjustment, and restricted permutations are
outside the 1.0 contract.

## Global statistic and null

Let \(\mathcal E=\{(i,j):i<j\}\). For edge \(e\), let \(p_{k,e}\) be its
observed proportion in group \(k\), and let

\[
p_e=\frac{\sum_k n_kp_{k,e}}{\sum_k n_k}.
\]

Define

\[
d_{k,e}=2p_{k,e}(1-p_{k,e}), \qquad
D_{k,e}=p_{k,e}+p_e-2p_{k,e}p_e.
\]

The unnormalized score is

\[
T=\sqrt m\sum_{k=1}^m\sqrt{n_k}
\left\{
\frac{n_k}{n_k-1}\sum_{e\in\mathcal E}d_{k,e}
-
\frac{n}{n-1}\sum_{e\in\mathcal E}D_{k,e}
\right\}.
\]

Each undirected edge is counted once. The normalization constant from the
asymptotic version of the Fraiman--Fraiman statistic is omitted because a
positive constant does not change a randomization ordering.

The finite-sample randomization null is that whole-graph observations are
exchangeable with respect to group labels (in the standard independent-groups
setting, equality of the group graph laws). The score is targeted to changes in
marginal edge frequencies. It is not an omnibus detector of alternatives that
change only cross-edge dependence, motifs, or other topology while preserving
all marginal edge probabilities.

Use the predeclared two-sided extremeness statistic \(S=|T|\). This avoids a
directional assumption for unbalanced and multi-group designs.

For exhaustive enumeration of all distinct fixed-size label assignments,
including the observed assignment,

\[
p_{\mathrm{exact}}=
\frac{\#\{S_\pi\geq S_{\mathrm{obs}}\}}{N_{\mathrm{assign}}}.
\]

For \(B\) random fixed-size assignments sampled uniformly with replacement
from the complete orbit (so the observed assignment can be sampled),

\[
p_{\mathrm{MC}}=
\frac{1+\#\{S_b\geq S_{\mathrm{obs}}\}}{B+1}.
\]

Ties are included. Assignment duplicates and equal statistics are never
deduplicated. Floating-point equivalents are treated as ties with tolerance
\(\sqrt{\epsilon}\max(1, |T_{\mathrm{obs}}|, \max_b|T_b|)\), where
\(\epsilon\) is machine double precision; the realized tolerance is stored in
the result. Automatic exhaustive enumeration is used when
\(N_{\mathrm{assign}}\leq 10{,}000\). Results report the method, number of
assignments, attainable p-value resolution, and number of ties. Degenerate
data for which every assignment has the same statistic return \(p=1\).

## Edge-level inference

Edge tests are a separate family of hypotheses and are not gated by, or
interpreted as a decomposition of, the global result. The multiplicity family
contains every upper-triangle edge; no data-dependent prefiltering is applied.

The default raw test is the two-sided Fisher exact test on the \(m\times2\)
present/absent table for each edge. Holm adjustment is the default because it
controls family-wise error under arbitrary dependence between edge tests. BY
is available for false-discovery-rate control under arbitrary dependence. BH
is opt-in and requires the user to accept its positive-dependence assumptions.

For two groups, the reported signed effect is
\(p_{2,e}-p_{1,e}\), using constructor group order. For more than two groups,
the reported effect is \(\max_kp_{k,e}-\min_kp_{k,e}\), together with the group
names attaining the extrema. Omnibus multi-group results make no pairwise
inferential claim.

An edge is `selected` when its adjusted p-value is at most `alpha`. Selected
nodes are a descriptive incidence summary and have no node-level error
guarantee.

## Exploratory ablation

The optional ablation path orders edges deterministically by raw p-value,
decreasing absolute effect, then endpoint indices. Prefix-sum updates report
the residual score and a `descriptive_tail_fraction`.

Because the order is selected from the same labels, that tail fraction is not
a calibrated post-selection p-value. The path returns no first-non-rejecting
set and supports no claim of minimality, equivalence, necessity, causality, or
statistical indistinguishability.

## Validation gate

The paper's complete-null calibration cells use 2,000 independent
replications, and its partial edge-null cells use 1,000, at
\(\alpha=0.05\). The one-sided 95% Wilson upper confidence limit must not
exceed 0.065 for:

1. global rejection under whole-graph exchangeability;
2. Holm family-wise error under complete and partial edge nulls.

BY false-discovery behavior under partial nulls uses a prespecified one-sided
confidence bound for mean false-discovery proportion. These tolerances are
implementation validation criteria, not proofs of exactness; exactness rests
on the randomization argument and exact-oracle tests.

