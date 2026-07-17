# Simulation and benchmark design

The inferential simulation cells and their pass/fail criteria were frozen
before the full results were generated. The prefix-sum timing benchmark was
added during the brain-first manuscript restructuring and uses a separate
deterministic seed stream; it does not change any inferential simulation cell.
Changes to the inferential design require a new design version, a new master
seed, and explicit disclosure in the manuscript.

## Common settings

- Significance level: 0.05.
- Global randomization assignments: 199 Monte Carlo assignments per replicate.
- Global statistic: two-sided extremeness `abs(T)`.
- Edge family: all upper-triangle edges.
- Primary edge procedure: two-sided Fisher tests with Holm adjustment.
- Secondary edge procedure: BY adjustment of the same raw p-values.
- Descriptive comparator: unadjusted raw p-values.
- Random-number generator: L'Ecuyer-CMRG.
- Master seed: 20260714.
- Graph observations are independent between subjects. “Dependent edges”
  always means dependence among edges within one subject-level graph.

## Complete-null calibration

Each cell uses 2,000 independent simulation replicates and six nodes (15 tested
edges).

1. `independent_balanced`: two groups of 10, independent Bernoulli edges with
   probability 0.30.
2. `independent_unbalanced`: group sizes 7 and 13, independent Bernoulli edges
   with probability 0.30.
3. `latent_density`: two groups of 10; each graph draws a density from a beta
   distribution with mean 0.30 and concentration 20, then draws conditionally
   independent edges. This induces positive within-graph edge dependence.
4. `fixed_density`: two groups of 10; every graph contains exactly five
   uniformly sampled edges, inducing negative within-graph edge dependence.
5. `sparse_three_group`: three groups of 8, independent Bernoulli edges with
   probability 0.05.
6. `dense_three_group`: three groups of 8, independent Bernoulli edges with
   probability 0.95.

Primary outcomes are global rejection probability, Holm family-wise error
probability, and BY false-discovery proportion. With no true alternatives, BY
false-discovery proportion equals the indicator of any rejection.

For global rejection and Holm family-wise error, the implementation-validation
criterion is a one-sided 95% Wilson upper confidence limit no greater than
0.065. This tolerance is not a proof of exactness.

## Partial-null and power cells

Each cell uses 1,000 independent simulation replicates, eight nodes (28 tested
edges), and 20 graphs per group unless noted.

1. `localized_positive`: the first six edges have probabilities 0.20 and 0.60
   in Groups A and B; all remaining edges have probability 0.20.
2. `scattered_mixed`: six edges at separated positions in the deterministic
   upper-triangle ordering differ; three shift from 0.20 to 0.60 and three
   shift from 0.60 to 0.20. Null edges have probability 0.20.
3. `latent_dependent`: graph-specific beta density induces dependence among
   null edges; six alternative edges have probabilities 0.20 and 0.60.
4. `three_group`: three groups of 15; six edges have probabilities 0.20, 0.40,
   and 0.65, while all null edges have probability 0.20.

An edge is truly alternative when at least one group marginal differs.

For each adjustment procedure, report:

- true-positive rate `TP / number_of_true_alternatives`;
- precision `TP / max(number_selected, 1)`;
- false-discovery proportion `FP / max(number_selected, 1)`;
- family-wise error indicator `I(FP > 0)`;
- signed set-size error `number_selected - number_of_true_alternatives`; and
- Jaccard index `|selected intersect truth| / |selected union truth|`.

Metrics are reported unconditionally. A separate, explicitly labeled table
reports edge metrics conditional on global rejection; it never replaces the
unconditional table.

Holm partial-null family-wise error uses a one-sided 95% Wilson upper confidence
limit with tolerance 0.065. BY mean false-discovery proportion uses a
one-sided 95% normal-bootstrap upper confidence limit with tolerance 0.065.
Any failed criterion triggers a rerun with a diagnostic seed set and blocks the
corresponding validity claim.

## Weak-null scope stress test

Two groups have the same edge marginals (0.50) but different dependence:
Group A alternates between empty and complete graphs, while Group B draws
independent Bernoulli edges. This violates the global whole-graph
exchangeability null while satisfying every edge-marginal null.

No global calibration criterion applies. Holm edge family-wise error remains a
required pass outcome because all edge-marginal nulls are true.

## Sensitivity

On `localized_positive`, use 500 common-data replicates to compare 99, 199,
999, and 4,999 global assignments. Record global decisions, p-values, elapsed
time, and edge selections. Edge results should be invariant because they do not
use global permutations.

Node relabeling and group-order invariance are exact unit-test properties, not
Monte Carlo study outcomes.

## Performance

For the global permutation engine, use two groups of 20 and 499 assignments at
16, 32, 64, 128, and 256 nodes. Use 10 independent repetitions per size.
Compare the chunked matrix implementation with a direct assignment loop using
the identical assignment matrix.

For the ablation update, use two groups of 20, 99 identical assignments, and
16, 24, 32, 40, 48, 56, and 64 nodes. Use 10 independent repetitions per
size. Compute the complete assignment-by-edge contribution matrix once, then
compare exact prefix-sum updates with repeated direct summation of every
remaining edge tail. This isolates the update optimized by the package and is
not an end-to-end runtime comparison.

For both benchmarks, verify numerical identity before accepting a timing,
warm up each implementation, randomize measurement order, and time calibrated
batches lasting at least 0.25 seconds when feasible. Report per-call medians
and interquartile ranges. For the global engine, also report the measured
input-object footprint; do not describe it as peak memory. Complexity
statements include the number of observations.

## Failure policy

- Do not replace failed cells silently.
- Confirm a failure with a second disjoint seed stream.
- Inspect exact-oracle identities and Monte Carlo resolution first.
- Change a method or claim only by versioning this design and regenerating all
  dependent results, figures, tables, and manuscript macros.
