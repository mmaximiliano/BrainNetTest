# Version 0.2.1 reference

The pre-redesign implementation used in historical comparisons is frozen by:

* Git tag: `v0.2.1`
* Commit: `5fce87a`
* CRAN publication date: 2026-06-05

Replication code that compares version 1.0.0 with the historical adaptive
procedure must retrieve source from that immutable tag. The legacy procedure
is not copied into or exported by version 1.0.0 because its strict-tail
permutation estimate and data-adaptive stopping path do not satisfy the new
inferential contract.
