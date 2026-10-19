# Régine 2026

**REGINE**: *Resampling-based Engine for Gene-set Inference and Network Enrichment*. First described and used in
Beckers A, Vantroeyen O, Cazevieille C, Adam H, Krouk G, Gantet P, Laplaze L, Guyomarc'h S (2026) PUCHI controls MYB94
and MYB96 to form the root cap cuticle of developing lateral roots. *Plant Physiol*, in press.
[doi:10.1093/plphys/kiag670](https://doi.org/10.1093/plphys/kiag670)

**Which transcription factors could regulate your genes?** Give Régine a list of Arabidopsis genes (a cluster,
co-expressed or differentially expressed genes). It tells you whether the genes of your list share regulators more than
random genes do, picks the kernel size at which they differ most, and ranks candidate regulators with their annotation,
from two networks:

* **binding**: the DAP-seq cistrome of 387 transcription factors (O'Malley et al. 2016);
* **regulation** (Régine 2.0): the 65 TARGET experiments of ConnecTF, 57 TFs (30 not in DAP-seq), with the direction:
  does the TF *induce* or *repress* the genes of your list?

and side by side, one row per TF, with the genes both bound and regulated and a combined P.

**Online program: <https://gabkrouk.github.io/Regine/>** (runs entirely in your browser; your list is never sent anywhere).

Régine and [GeneCloud](https://gabkrouk.github.io/genecloud/) talk to each other (click a GO term → candidate TFs): see [INTEGRATION.md](INTEGRATION.md).

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

### TARGET (Régine 2.0)

TARGET (Bargmann et al. 2013) moves a TF fused to the glucocorticoid receptor into the nucleus of protoplasts with
dexamethasone, with cycloheximide so that only its direct targets respond; the genes induced or repressed after 3 h
(DESeq2, FDR < 5%) are its targets. The data are the 65 TARGET experiments of the ConnecTF data release v1 (Brooks et al.
2019, 2021; Varala et al. 2018; Alvarez et al. 2020): 59 in root and 6 in shoot protoplasts, pooled.

* Each experiment is a column of the network and is tested on all its regulated genes (kernels, score, P), on its
  induced genes and on its repressed genes. *Activates*: the induced genes are enriched in your list (FDR ≤ 5%);
  *represses*: the repressed genes are.
* A TF tested several times (CDF1, LBD37, LBD38, TGA1, TGA4, CRF4, ERF060, SNZ) keeps one row per experiment, because
  repeated experiments overlap little (Jaccard 0.05 to 0.22).
* The files list the significant genes only, so the TARGET background is the 18,727 genes regulated in at least one
  experiment (or your own background).
* Side by side: per TF, the DAP-seq P and the P of its best TARGET experiment multiplied by its number of experiments,
  combined by Fisher's method (FDR over all TFs); plus the enrichment of the genes both bound in DAP-seq and regulated in
  TARGET (27 TFs are in both).
* A TARGET experiment regulates 2,000 genes on median, so TARGET is less specific than DAP-seq; very strong hits can come
  from a shared module (for example ribosomal protein genes, see the case study).

## Case study: cell identity in the root single-cell atlas

Run on the marker genes of the 14 cell types of the Shahan et al. (2022) root atlas, Régine recovers the identity
regulators from their targets alone: the VNDs in protoxylem (FDR 10⁻³²), the PEAR DOFs in protophloem (10⁻¹⁹) and
SOMBRERO in the root cap. The kernel AUC flags the cell types whose regulators are absent from the DAP-seq set.
With TARGET, NAC096 (protophloem) and MYB61 (procambium, endodermis) are both bound and regulated by their TF and
induce the markers, LBD38 represses the metaphloem markers, and the cytokinin response factors CRF4 and CRF2 explain the
metaxylem markers, which DAP-seq did not.
Details, data, scripts and figure: [`case_studies/shahan2022_root_atlas/`](case_studies/shahan2022_root_atlas/).

## Python (Régine 2.0, both networks)

Python 3.8 or later, standard library only (matplotlib, if installed, draws the AUC profile).

```sh
python3 python/regine.py data/example_nitrate_list.txt --name nitrate
# -> results/nitrate.dapseq.tsv, .target.tsv, .both.tsv, .profile.tsv
python3 python/regine.py my_list.txt --background expressed_genes.txt --network target --kernel 4 --name my_list
```

```python
import sys; sys.path.insert(0, "python")
import regine
data = regine.load("data")
res = regine.run(data, regine.read_genes("data/example_nitrate_list.txt"))
res["both"][0]["symbol"], res["both"][0]["combined_fdr"]     # ('NLP7', 6.9e-15)
```

The random kernels of the profile are drawn with Python's generator, so the AUCs differ slightly from the web version;
everything else (real, noise, score, P, FDR) is exact and identical.

## R (DAP-seq only)

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
| `python/regine.py` | Python implementation, DAP-seq and TARGET |
| `R/regine.R` | R implementation (DAP-seq) |
| `web/` | the online program (static site: `index.html`, `app.js`, `engine.js`, `data/`) |
| `data/dapseq_network.tsv.gz` | DAP-seq network, TF → target, 2 848 928 distinct links (the 2022 table without its 836 598 duplicate lines) |
| `data/tf_annotation.tsv` | the 387 TFs: symbol, family, number of targets, UniProt name, Araport11 description, TAIR summary, UniProt function |
| `data/gene_annotation.tsv.gz` | symbol and description of every Arabidopsis gene |
| `data/target_network.tsv.gz` | TARGET network: experiment, TF, target, log2FC, padj, induced / repressed (139,241 lines) |
| `data/target_experiments.tsv` | the 65 TARGET experiments: TF, symbol, family, in DAP-seq, source, tissue, treatment, cutoff, DOI, counts |
| `tools/build_data.mjs` | rebuilds `data/` and `web/data/` from the 2022 DAP-seq table and the GeneCloud annotation |
| `tools/build_target.mjs` | rebuilds the TARGET files from the ConnecTF data release v1 (`arabidopsis/metadata` and `arabidopsis/data`) |
| `tests/` | `node --test tests/engine.test.mjs`, `python3 tests/test_regine.py` and `Rscript tests/test_regine.R` (includes a check that the 2022 procedure converges to the exact values) |

Run the web version locally: `cd web && python3 -m http.server`, then open <http://localhost:8000>.
The site is published from `web/` to the `gh-pages` branch by `.github/workflows/pages.yml` on every push to `main`.

## Data and citation

* DAP-seq: O'Malley RC, Huang SC, Song L, Lewsey MG, Bartlett A, Nery JR, Galli M, Gallavotti A, Ecker JR (2016)
  Cistrome and epicistrome features shape the regulatory DNA landscape. *Cell* 165:1280–1292.
  [doi:10.1016/j.cell.2016.04.038](https://doi.org/10.1016/j.cell.2016.04.038)
* TARGET: Brooks MD, Juang C-L, Katari MS, Alvarez JM, Pasquino A, Shih H-J, Huang J, Shanks C, Cirrone J, Coruzzi GM
  (2021) ConnecTF: a platform to integrate transcription factor–gene interactions and validate regulatory networks.
  *Plant Physiol* 185:49–66. [doi:10.1093/plphys/kiaa012](https://doi.org/10.1093/plphys/kiaa012), data release v1.
  Experiments from Brooks et al. (2019) *Nat Commun* 10:1569, Varala et al. (2018) *PNAS* 115:6494 and Alvarez et al.
  (2020) *Nat Commun* 11:1157; method: Bargmann et al. (2013) *Mol Plant* 6:978.
* Annotation: TAIR / Araport11 and UniProtKB (CC BY 4.0), as assembled for [GeneCloud](https://github.com/GabKrouk/genecloud).
  TF families are inferred from the annotation text (`tools/families.mjs`).
* Code: MIT licence, G. Krouk (2022, 2026).
