# Régine 2026

**Which transcription factors could regulate your genes?** Give Régine a list of Arabidopsis genes (a cluster,
co-expressed or differentially expressed genes). It tells you whether the genes of your list share DAP-seq regulators
more than random genes do, picks the kernel size at which they differ most, and ranks the 387 transcription factors of
the DAP-seq cistrome (O'Malley et al. 2016) as candidate regulators, with their annotation.

**Online program: <https://gabkrouk.github.io/Regine/>** (runs entirely in your browser; your list is never sent anywhere).

Analyse du code de 2022 et des changements (en français) : [ANALYSE.md](ANALYSE.md).

## How it works

DAP-seq is permissive: a median TF binds 22% of the genes, and two random genes already share about 28 TFs.
Régine therefore looks at **kernels**, groups of *k* genes of your list, and asks which TFs bind *all* of them.
For a TF that binds *k<sub>t</sub>* of the *n* genes of your list and *K<sub>t</sub>* of the *N* background genes:

| | |
|---|---|
| real | C(k<sub>t</sub>, k) / C(n, k): fraction of the kernels of your list entirely bound by the TF |
| noise | C(K<sub>t</sub>, k) / C(N, k): same for kernels of random background genes |
| **Régine score** | real − noise (the *Real − Noise* of Régine 2022, computed exactly instead of sampled) |
| fold, P, FDR | enrichment of the TF's targets in your list, hypergeometric P, Benjamini–Hochberg over the 387 TFs |

**Automatic kernel size.** For *k* = 2…10, random kernels are drawn from your list and from the background and the
number of TFs they share is counted (the boxplots of Régine 2022). Their separation is measured by the AUC, the probability
that a kernel of your list shares more TFs than a random kernel; the *k* with the highest AUC is used.

The score favours TFs that bind almost all of your genes (often broad binders); the P favours TFs whose targets are rare
in the genome. Both rankings are available; the table is ranked by FDR by default. On the 2022 nitrate example, NLP7 binds 10 of the 13 genes and only 829 genes
genome-wide (×30, FDR 1e−11).

## R

Base R only, no package to install.

```r
source("R/regine.R")
reg <- regine_load("data")                                      # a few seconds
res <- regine(reg, "data/example_nitrate_list.txt", name = "nitrate")
# -> results/nitrate.pdf (kernel profile), results/nitrate.regulators.tsv, results/nitrate.profile.tsv
head(res$regulators)

# your own background (recommended: genes expressed in the experiment), fixed kernel, ranked by Régine score
res <- regine(reg, "my_list.txt", background = "expressed_genes.txt", kernel = 6, rank_by = "score", name = "my_list")
```

`Rscript R/run_regine.R` runs the example. The original 2022 scripts are kept in [`R/legacy_2022/`](R/legacy_2022/).

## Repository

| | |
|---|---|
| `R/regine.R` | R implementation |
| `web/` | the online program (static site: `index.html`, `app.js`, `engine.js`, `data/`) |
| `data/dapseq_network.tsv.gz` | DAP-seq network, TF → target, 2 848 928 distinct links (the 2022 table without its 836 598 duplicate lines) |
| `data/tf_annotation.tsv` | the 387 TFs: symbol, family, number of targets, UniProt name, Araport11 description, TAIR summary, UniProt function |
| `data/gene_annotation.tsv.gz` | symbol and description of every Arabidopsis gene |
| `tools/build_data.mjs` | rebuilds `data/` and `web/data/` from the 2022 DAP-seq table and the GeneCloud annotation |
| `tests/` | `node --test tests/engine.test.mjs` and `Rscript tests/test_regine.R` (includes a check that the 2022 procedure converges to the exact values) |

Run the web version locally: `cd web && python3 -m http.server`, then open <http://localhost:8000>.
The site is published from `web/` to the `gh-pages` branch by `.github/workflows/pages.yml` on every push to `main`.

## Data and citation

* DAP-seq: O'Malley RC, Huang SC, Song L, Lewsey MG, Bartlett A, Nery JR, Galli M, Gallavotti A, Ecker JR (2016)
  Cistrome and epicistrome features shape the regulatory DNA landscape. *Cell* 165:1280–1292.
  [doi:10.1016/j.cell.2016.04.038](https://doi.org/10.1016/j.cell.2016.04.038)
* Annotation: TAIR / Araport11 and UniProtKB (CC BY 4.0), as assembled for [GeneCloud](https://github.com/GabKrouk/genecloud).
  TF families are inferred from the annotation text (`tools/families.mjs`).
* Code: MIT licence, G. Krouk (2022, 2026).
