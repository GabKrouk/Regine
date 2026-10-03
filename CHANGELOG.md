# Changelog

## 2026 — Régine 2026

* Bridge with GeneCloud: click a GO term in GeneCloud to rank the TFs most associated with it; send TF targets back to GeneCloud ([INTEGRATION.md](INTEGRATION.md)).
* Exact real / noise / score (no sampling), hypergeometric P and FDR, ranking by score or P.
* Automatic kernel size (maximal AUC between kernels of the list and random kernels).
* TF annotation: symbol, family, UniProt name, Araport11 description, TAIR summary, UniProt function.
* Background: DAP-seq network genes, protein-coding genes or a custom list.
* DAP-seq table deduplicated (836 598 duplicate lines removed).
* R code rewritten in base R without global variables; web version on GitHub Pages.

## 2022 — Régine Pack 2022

* Original R scripts (`R/legacy_2022/`).
