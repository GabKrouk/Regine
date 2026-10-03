# Régine 2026 — R tests. Run from the repository root:  Rscript tests/test_regine.R
source("R/regine.R")
reg <- regine_load("data")
genes <- regine_genes("data/example_nitrate_list.txt")
ok <- function(cond, msg) { if (!isTRUE(cond)) stop("FAILED: ", msg); cat("ok -", msg, "\n") }

ok(length(genes) == 13 && all(grepl("^AT[1-5]G", genes)), "gene list parsed (13 AGIs, upper case)")
ok(ncol(reg$m) == 387 && sum(reg$m) == 2848928, "network: 387 TFs, 2 848 928 distinct links")

# 1. exact "real" = enumeration of every kernel of the list
k <- 4
tab <- regine_regulators(reg, genes, k)
kern <- utils::combn(genes, k)
enum <- rowMeans(apply(kern, 2, function(g) colSums(reg$m[g, , drop = FALSE]) == k))
ok(max(abs(enum[tab$agi] - tab$real)) < 1e-12, "real = fraction of all C(13,4) kernels bound (exact enumeration)")

# 2. the Régine 2022 procedure (Monte Carlo, Find.me.direct.regulators.KERNEL.STAT.gb) converges to the exact values
set.seed(7)
fast <- split(rep(colnames(reg$m), colSums(reg$m)), rownames(reg$m)[which(reg$m, arr.ind = TRUE)[, 1]])  # 2022 memoization
direct <- function(g) Reduce(intersect, fast[g])
all_genes <- reg$network
rep <- 20000; k <- 6
real <- noise <- setNames(numeric(387), colnames(reg$m))
for (i in seq_len(rep)) {
  r <- direct(sample(genes, k)); real[r] <- real[r] + 1
  b <- direct(sample(all_genes, k)); noise[b] <- noise[b] + 1
}
tab6 <- regine_regulators(reg, genes, k)
ok(max(abs(real[tab6$agi] / rep - tab6$real)) < 0.015, "2022 Monte Carlo 'Real' / rep matches the exact real (k = 6)")
ok(max(abs(noise[tab6$agi] / rep - tab6$noise)) < 0.01, "2022 Monte Carlo 'Noise' / rep matches the exact noise (k = 6)")

# 3. statistics
ok(all(tab6$p >= 0 & tab6$p <= 1) && all(tab6$fdr >= tab6$p - 1e-15), "P values in [0, 1], FDR >= P")
nlp7 <- tab6[tab6$symbol == "NLP7", ]
ok(nlp7$your_genes_bound == 10 && nlp7$background_bound == 829, "NLP7 binds 10/13 genes of the nitrate list and 829 genes overall")
prof <- regine_profile(reg, genes, reps = 500)
ok(all(prof$profile$auc > 0.5) && prof$k_auto %in% 3:6, "kernel profile: list departs from random; automatic k in 3..6")
cat("all R tests passed\n")
