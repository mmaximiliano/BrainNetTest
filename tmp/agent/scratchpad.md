# BrainNetTest JSS rebuild scratchpad

## State

- Branch: `feat/brainnettest-jss-rebuild`
- Approved direction: breaking 1.0 redesign.
- Evidence: simulations only.
- Global inference: whole-graph exchangeability permutation test, targeted to
  marginal edge-frequency shifts.
- Edge inference: separate all-edge family, Fisher exact tests with Holm by
  default; selected nodes are descriptive.
- Legacy ranked removal: exploratory ablation only.
- Canonical graph distance: upper triangle only.
- Canonical randomization extremeness: `abs(T)`.

## Environment

- Homebrew R 4.6.1 installed on 2026-07-13.
- Development R packages and Pandoc are installed.
- TeX is not installed yet.

## Baseline

- Existing `testthat` suite: PASS.
- `R CMD build .`: PASS after installing Pandoc.
- `R CMD check BrainNetTest_0.2.1.tar.gz --no-manual`: PASS (0 errors,
  0 warnings, 0 notes).
- Manuscript/replication blockers remain the statically verified missing files,
  stale paths, and uninstalled TeX toolchain.

## Current checkpoint

- Contract frozen in `paper_source/statistical-contract.md`.
- New exact global test, Holm edge inference, descriptive ablation, S3 input
  and result classes, methods, validation, migration docs, vignettes, and
  automatic three-OS CI are implemented.
- New and retained generator tests pass.
- A 500-replicate smoke check gave global null rejection 0.044, Holm complete
  null FWER 0.004, and strong-alternative global power 0.994.
- Full local `R CMD check --as-cran`, including PDF manual and HTML
  validation: PASS with 0 errors, 0 warnings, and 0 notes.

Next: pre-register and run the publication validation study, then generate the
new paper from those validated artifacts.

## Final state

- Corrected Monte Carlo assignments sample the complete fixed-size orbit with
  replacement; a finite-orbit size regression test passes.
- Full prespecified validation regenerated with master seed 20260714 using the
  checked and installed 1.0.0 source package.
- All complete- and partial-null error gates pass.
- Full replication runtime: 322.496 seconds.
- Current JSS 3.6 assets are used.
- `article.Rnw` builds to an 11-page `article.pdf`; the cover letter and
  `code.html` also build.
- Read-only submission verification passes 29 required files, three figures,
  18 generated artifacts, and source/design/contract hashes.
- Three independent final audits returned PASS.
- Final local lint, tests, IDE diagnostics, and `R CMD check --as-cran` all
  pass with zero errors, warnings, or notes.
- Exact package archive:
  `9ffdfcff7bda432e593b7f6f434794f6861605cd3d4f4f35065e58f53abbb694`.
- Remaining external actions require explicit authorization: commit/push,
  CRAN submission, and JSS resubmission.

## Audit session 2026-07-13 (fresh review of rebuild)

Verified from scratch: `R CMD check --as-cran` OK (0/0/0), testthat all pass,
`lintr` 0 lints, `code-full.R` reproduces (322s; statistical numbers identical;
only timing-derived speedup macro differs 17.9 vs 17.2, disclaimed), `verify.R`
passes, `article.Rnw` re-knits byte-identical, PDF compiles.

Statistic confirmed against Fraiman & Fraiman (2018) eq. 9: d_{k,e}=2p(1-p),
D_{k,e}=p+p_e-2 p p_e, weights sqrt(n_k) n_k/(n_k-1) and sqrt(n_k) n/(n-1) with
sqrt(m) prefactor. Correct. Two-sided |T| randomization is valid.

Changes made (manuscript/replication only; package R/ untouched):
- Section 5 now demonstrates the S3 system: prints `networks`,
  `summary(networks)`, and adds a `plot(result)` figure (three-group, fills a
  2x2 panel grid). Directly rebuts "no suitable use of classes/methods".
- JSS style: section-3 title wraps `\proglang{R}` with bookmark-safe arg.
- Cited tooling: added `knitr`, `rmarkdown`, `testthat` bib entries + `\citep`
  in Computational details; `\pkg{NBR}`/`\proglang{R}` markup in NBR title.
- New figure `figure/plot-1.pdf` added to `SUBMISSION-FILES.txt`; README and
  response-to-editor updated to describe it.
- article.tex / article.pdf regenerated (13 pages). Cleaned build intermediates.

Top residual acceptance risk (author decision, not fixed): manuscript is
simulation-only with no real-data illustration. Defensible but the single most
likely reviewer ask.

## ABIDE real-data example added 2026-07-13

Implemented the real-data illustration (rebuts "narrow application context").
- New `paper_source/data-application.R`: downloads ABIDE I PCP derivatives
  (NYU site, CPAC filt_noglobal, AAL atlas) from the public S3 bucket, QC by
  mean FD < 0.2, builds one binary connectome per subject by 10% proportional
  threshold on Pearson correlations, validates with `brainnet_data()`, saves
  `application/abide_connectomes.rds` (+ provenance, subjects). Deterministic:
  re-runs are byte-identical; md5 5d632933... . Embedded AAL code->name map.
- Real result: 98 control vs 73 autism, 114 nodes, 6441 edges. Global T=-44.69,
  p=0.001 (rejects). Holm and BY select 0 edges (min Holm p 0.108) -- an honest
  two-layer outcome: global signal detectable, edge localization not, under
  strict multiplicity at this n. Top raw edge Frontal_Sup_Medial_R--Rectus_R.
- New manuscript Section 6 "Illustration on resting-state brain connectomes":
  loads the rds, runs brainnet_test live (network-free knit), shows result,
  top-edges table (anatomical names), and a volcano figure
  (figure/abide-volcano-1.pdf). Abstract, discussion, computational details,
  intro roadmap updated. Manuscript now 15 pages.
- Refs added: DiMartino2014, Craddock2013, TzourioMazoyer2002, vandenHeuvel2017.
- Wiring: verify.R checks application artifacts + connectomes md5 (33 files);
  README, SUBMISSION-FILES, response-to-editor, cover letter, .gitignore
  updated. Derived data shipped under ABIDE CC BY-NC-SA with attribution.
- Manuscript reads only the small derived rds; only data-application.R needs
  network. Full knit+compile+verify all pass; data-application.R lints clean.
- Author decision left open: ship derived rds (current, reproducible) vs
  download-only. Current choice is documented and license-compliant.

## HCP-by-sex request 2026-07-13

Researched: HCP per-subject connectivity is NOT credential-free and HCP Data
Use Terms PROHIBIT redistribution of data/derivatives. So a self-contained,
reproducible in-paper HCP example is impossible (would need every reviewer to
have a ConnectomeDB account, or would breach DUT / fabricate numbers). Did not
ship HCP data or fabricate results.
Delivered instead:
- `data-application-hcp.R`: companion script, downloads NOTHING; reads a user's
  own HCP PTN netmats (netmats1.txt + subjectIDs.txt) + unrestricted behavioral
  (Subject, Gender), builds one binary connectome per subject at 10% density,
  groups healthy adults by sex, runs the same two-layer test. Guards missing
  inputs; tested end-to-end on synthetic HCP-format files; 0 lints.
- Manuscript Section 6 paragraph: states the HCP/Fraiman sex-difference
  motivation, explains why HCP can't be bundled (DUT), points to the companion
  script; cites VanEssen2013 + Smith2013. No fabricated HCP numbers.
- Also confirmed NYU healthy-controls-by-sex is feasible on ABIDE (72M/26F) if
  we ever want an in-paper sex example on openly licensed data.
- Wiring: refs.bib (+VanEssen2013,+Smith2013), README (HCP section), verify.R
  (34 files), SUBMISSION-FILES, response letter updated. Manuscript 16 pages;
  knit+compile+verify all pass; both app scripts lint clean.

### HCP REMOVED per user request 2026-07-13

User: "Since HCP is not accessible, remove it. We should not mention it."
Reverted everything HCP: deleted data-application-hcp.R; removed the Section 6
HCP paragraph; removed VanEssen2013 + Smith2013 from refs.bib; reverted verify.R
(back to 33 files), SUBMISSION-FILES, README (dropped HCP subsection), and the
response letter. Grep across repo (excl. tmp) for HCP/Human Connectome/netmats/
ConnectomeDB/Van Essen/Smith2013 = no matches. Manuscript back to 15 pages;
knit+compile+verify all pass. ABIDE illustration unchanged and intact.

## Brain-first clinical restructuring 2026-07-14

User corrected the earlier strategy: keep the brain-only/domain framing and
real clinical claims, modeled on accepted JSS domain/method/software papers
(FIAR, missSBM, GNAR, rags2ridges, fastnet). Implemented:

- Package/documentation title and DESCRIPTION are brain-first again:
  "Global and Edge-Wise Inference for Brain Network Populations." Root README,
  vignettes, CITATION, NEWS, function title, and paper contract now describe
  subject-level functional/structural connectomes and clinical groups.
- Manuscript title/abstract/intro/contract/discussion restored a clinical
  connectomics narrative. Related-software table now includes FIAR. General
  network applicability remains secondary, not the lead framing.
- Found the original full ABIDE groups were imbalanced on sex (p=.034) and
  motion (p=.00023). Added an outcome-free deterministic coarsened balance
  sensitivity sample: exact sex, 5-year age bins, 0.05 FD bins; 55+55
  participants (49M/6F each), mean age 14.25 vs 14.75, mean FD .06175 vs
  .06183. Full global result T=-44.69,p=.001; balanced T=-27.39,p=.01.
  Both have 0 Holm/BY selected edges. Clinical claim is therefore population
  association/edge-pattern change, not causal, diagnostic, or a biomarker.
- data-application.R now saves full + balanced populations and balance fields;
  deterministic artifact md5 verified. Application figure replaced generic
  volcano with 2-panel anatomical-system heatmap + edge test plot.
- Restored the algorithmic contribution: exact complexity explanation plus
  new prefix-sum benchmark (10 reps, 99 assignments, 16--64 nodes). All paths
  checked numerically identical. Final largest median update-stage speedup:
  953x at 2016 edges (explicitly NOT end-to-end). Global chunk engine largest
  median speedup: 19.1x. Full final replication: 499.66 sec.
- Updated simulation-design disclosure, README, LICENSE, SUBMISSION-FILES,
  verify.R, response-to-editor, cover letter, refs, and generated manifest.
- Final article: 19 pages; no undefined refs or overfull boxes. Final verify:
  44 required files, 4 manuscript figures, 21 generated artifacts, ABIDE hash
  and 55/group balance checks PASS.
- Final R CMD check --as-cran (including PDF/HTML manuals): 0 errors, 0
  warnings, 0 notes. All tests and lints pass.
- Final source archive: 279105 bytes,
  sha256 905266d368c0be914f1d7cc10c06d1811d4f505ccfc5df2c60cd80d82e589f14.
- Follow-up from both readonly reviews: softened ABIDE wording throughout the
  abstract, results, discussion, response, and cover letter. Final claim is an
  extreme statistic under autism-control label randomization and a
  preprocessing-conditional population association, not a general autism
  effect. Kept the clinical case study and brain-first framing.

