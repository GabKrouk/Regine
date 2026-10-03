# Régine 2026 — example run. Execute from the repository root:  Rscript R/run_regine.R
source("R/regine.R")
reg <- regine_load("data")

# 1. one call: kernel profile (PDF), automatic kernel size, candidate regulators (TSV)
res <- regine(reg, "data/example_nitrate_list.txt", name = "nitrate_2022")
print(head(res$regulators[, c("rank", "agi", "symbol", "family", "your_genes_bound", "real", "noise", "regine_score", "fdr")], 15))

# 2. same table ranked by enrichment significance instead of the Régine score
sig <- regine_regulators(reg, "data/example_nitrate_list.txt", k = res$k, rank_by = "p")
print(head(sig[, c("rank", "agi", "symbol", "family", "your_genes_bound", "background_bound", "fold", "fdr")], 10))

# 3. with your own background (recommended: the genes expressed in the experiment)
# res <- regine(reg, "my_list.txt", background = "expressed_genes.txt", name = "my_list")
