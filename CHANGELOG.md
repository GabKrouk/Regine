# Changelog

## 2026-10 — Régine 2.0: TARGET

* Second network: the 65 TARGET experiments of ConnecTF (57 TFs, 30 of them not in DAP-seq; root and shoot protoplasts,
  cycloheximide, FDR < 5%), built by `tools/build_target.mjs` into `data/target_*` and `web/data/target.*`.
* Each experiment is tested on all its regulated genes (kernels, score), on its induced and on its repressed genes:
  a TF *activates* or *represses* the genes of the list. Experiments of the same TF are kept separate.
* TARGET background: the 18,727 genes regulated in at least one experiment.
* Web: DAP-seq, TARGET and side-by-side views; side by side gives, per TF, DAP-seq, the best TARGET experiment, the genes
  bound and regulated, and a combined P (Fisher's method) with its FDR.
* Python version (`python/regine.py`, standard library only) with both networks; tests in `tests/test_regine.py`.
  The R version is unchanged (DAP-seq only).
* Case study: both networks on the root atlas cell types.

## 2026 — Régine 2026

* Case study on the Shahan et al. 2022 root single-cell atlas: Régine recovers the VND (protoxylem), PEAR DOF (protophloem)
  and SOMBRERO (root cap) regulators from cell type marker genes ([`case_studies/shahan2022_root_atlas/`](case_studies/shahan2022_root_atlas/)).

* Exact real / noise / score (no sampling), hypergeometric P and FDR, ranking by score or P.
* Automatic kernel size (maximal AUC between kernels of the list and random kernels).
* TF annotation: symbol, family, UniProt name, Araport11 description, TAIR summary, UniProt function.
* Background: DAP-seq network genes, protein-coding genes or a custom list.
* DAP-seq table deduplicated (836 598 duplicate lines removed).
* R code rewritten in base R without global variables; web version on GitHub Pages.

## 2022 — Régine Pack 2022

* Original R scripts (`R/legacy_2022/`).
