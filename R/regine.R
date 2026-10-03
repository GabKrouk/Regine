###############################################################################################################
# Régine 2026 — candidate transcription factors regulating an Arabidopsis gene list, from the DAP-seq cistrome
# (O'Malley et al. 2016, Cell 165:1280). Gabriel Krouk (2022), rewritten 2026.
#
# Base R only (stats, graphics, grDevices, utils). Usage:
#
#   source("R/regine.R")
#   reg <- regine_load("data")                          # a few seconds, once per session
#   res <- regine(reg, "data/example_nitrate_list.txt", name = "nitrate")
#   head(res$regulators)                                # candidate TFs with annotation
#
# What it computes (identical to the web version, https://gabkrouk.github.io/Regine/):
#   for a list L of n genes, a background B of N genes (L is added to B), and each TF t binding k_t genes of L and
#   K_t genes of B:
#     real  = choose(k_t, k) / choose(n, k)    fraction of the kernels of k genes of L entirely bound by t
#     noise = choose(K_t, k) / choose(N, k)    same for kernels of random background genes
#     score = real - noise                     ("Real - Noise" of Régine 2022, exact instead of sampled)
#     p     = P(X >= k_t), X ~ hypergeometric(N, K_t, n);  fdr = Benjamini-Hochberg over all TFs
#   kernel profile: for k in 2..10, `reps` random kernels of L and of B, number of TFs binding all their genes;
#   auc = P(list kernel shares more TFs than a random kernel) (+ ties / 2). The automatic k is the largest auc.
###############################################################################################################

# ---------------------------------------------------------------------------------------------- data
regine_load <- function(dir = "data") {
  net <- utils::read.delim(gzfile(file.path(dir, "dapseq_network.tsv.gz")), colClasses = "character")
  tfs <- utils::read.delim(file.path(dir, "tf_annotation.tsv"), colClasses = "character", quote = "")
  genes <- utils::read.delim(gzfile(file.path(dir, "gene_annotation.tsv.gz")), colClasses = "character", quote = "")
  ids <- genes$agi
  # gene x TF incidence matrix (TRUE = DAP-seq peak), the "memoization" of Régine 2022 done once with indexing
  m <- matrix(FALSE, nrow = length(ids), ncol = nrow(tfs), dimnames = list(ids, tfs$agi))
  m[cbind(match(net$target, ids), match(net$tf, tfs$agi))] <- TRUE
  tfs$n_targets <- as.integer(tfs$n_targets)
  structure(list(m = m, tfs = tfs, genes = genes,
                 protein_coding = ids[genes$protein_coding == "1"], network = ids[genes$in_network == "1"]),
            class = "regine_data")
}

# AGI identifiers from a file, a character vector or free text; upper case, unique, isoform suffix removed
regine_genes <- function(x) {
  if (length(x) == 1 && file.exists(x)) x <- readLines(x, warn = FALSE)
  x <- x[!grepl("^\\s*#", x)]
  ids <- toupper(unlist(regmatches(x, gregexpr("AT[1-5CMcm]G[0-9]{5}", x, ignore.case = TRUE))))
  unique(ids)
}

.regine_sets <- function(reg, genes, background) {
  genes <- regine_genes(genes)
  unknown <- setdiff(genes, rownames(reg$m))
  L <- setdiff(genes, unknown)
  B <- if (identical(background, "network")) reg$network
       else if (identical(background, "protein_coding")) reg$protein_coding
       else intersect(regine_genes(background), rownames(reg$m))
  added <- setdiff(L, B)
  if (length(unknown)) message(length(unknown), " identifier(s) not in Araport11: ", paste(head(unknown, 10), collapse = ", "))
  if (length(added)) message(length(added), " gene(s) of the list added to the background")
  if (length(L) < 2) stop("at least 2 recognised genes are needed")
  list(L = L, B = c(B, added), unknown = unknown, no_peak = L[rowSums(reg$m[L, , drop = FALSE]) == 0])
}

# choose(a, k) / choose(b, k), vectorised over a
.choose_ratio <- function(a, b, k) ifelse(a < k, 0, exp(lchoose(a, k) - lchoose(b, k)))

# ---------------------------------------------------------------------------------------------- kernel profile
regine_profile <- function(reg, genes, background = "network", kernel_range = 2:10, reps = 2000, seed = 1) {
  s <- .regine_sets(reg, genes, background)
  L <- match(s$L, rownames(reg$m)); B <- match(s$B, rownames(reg$m))
  kernel_range <- kernel_range[kernel_range >= 2 & kernel_range <= length(L)]
  cl <- colSums(reg$m[L, , drop = FALSE]); cb <- colSums(reg$m[B, , drop = FALSE])
  set.seed(seed)
  common <- function(pool, k) sum(colSums(reg$m[pool[sample.int(length(pool), k)], , drop = FALSE]) == k)
  draws <- list()
  prof <- do.call(rbind, lapply(kernel_range, function(k) {
    xl <- vapply(seq_len(reps), function(i) common(L, k), 0)
    xb <- vapply(seq_len(reps), function(i) common(B, k), 0)
    draws[[as.character(k)]] <<- data.frame(k = k, set = rep(c("list", "background"), each = reps), n_common = c(xl, xb))
    r <- rank(c(xl, xb))
    data.frame(k = k,
               auc = (sum(r[seq_len(reps)]) - reps * (reps + 1) / 2) / reps^2,
               mean_list = sum(.choose_ratio(cl, length(L), k)),      # exact expectations
               mean_background = sum(.choose_ratio(cb, length(B), k)),
               median_list = stats::median(xl), median_background = stats::median(xb))
  }))
  k_auto <- prof$k[which.max(prof$auc)]
  list(profile = prof, draws = do.call(rbind, draws), k_auto = k_auto, n = length(L), N = length(B),
       unknown = s$unknown, no_peak = s$no_peak, reps = reps)
}

# ---------------------------------------------------------------------------------------------- regulators
regine_regulators <- function(reg, genes, k, background = "network", rank_by = c("score", "p")) {
  rank_by <- match.arg(rank_by)
  s <- .regine_sets(reg, genes, background)
  n <- length(s$L); N <- length(s$B)
  ml <- reg$m[s$L, , drop = FALSE]
  kt <- colSums(ml); Kt <- colSums(reg$m[s$B, , drop = FALSE])
  real <- .choose_ratio(kt, n, k); noise <- .choose_ratio(Kt, N, k)
  p <- stats::phyper(kt - 1, Kt, N - Kt, n, lower.tail = FALSE)
  a <- reg$tfs
  out <- data.frame(agi = a$agi, symbol = a$symbol, family = a$family,
                    your_genes_bound = kt, n_list = n, background_bound = Kt, N_background = N,
                    real = real, noise = noise, regine_score = real - noise,
                    fold = (kt / n) / (Kt / N), p = p, fdr = stats::p.adjust(p, "BH"),
                    tf_targets_genome = a$n_targets, protein_name = a$protein_name,
                    short_description = a$short_description, summary = a$summary,
                    your_genes_bound_ids = apply(ml, 2, function(b) paste(s$L[b], collapse = ",")),
                    stringsAsFactors = FALSE, row.names = NULL)
  o <- if (rank_by == "score") order(-out$regine_score, out$p, -out$your_genes_bound) else order(out$p, -out$regine_score)
  out <- out[o, ]
  out$rank <- seq_len(nrow(out))
  rownames(out) <- NULL
  out[, c("rank", setdiff(names(out), "rank"))]
}

# ---------------------------------------------------------------------------------------------- plot
regine_plot <- function(prof, file = NULL, title = "Régine") {
  if (!is.null(file)) { grDevices::pdf(file, width = 9, height = 7.5); on.exit(grDevices::dev.off()) }
  op <- graphics::par(mfrow = c(2, 1), mar = c(4, 4.5, 3, 1), las = 1); on.exit(graphics::par(op), add = TRUE)
  d <- prof$draws
  d$set <- factor(d$set, levels = c("list", "background"))
  cols <- c("#b58fd8", "#d9dde1")
  graphics::boxplot(log1p(n_common) ~ set + k, data = d, col = cols, border = c("#5e2d8c", "#5d6670"),
                    yaxt = "n", xaxt = "n", outline = FALSE, ylab = "common TFs per kernel", xlab = "kernel size",
                    main = paste0(title, ": kernels of your list (purple) vs random genes (grey)"), cex.main = 1)
  ks <- prof$profile$k
  graphics::axis(1, at = seq(1.5, by = 2, length.out = length(ks)), labels = ks)
  tk <- c(0, 1, 2, 5, 10, 20, 50, 100, 200, 400)
  graphics::axis(2, at = log1p(tk), labels = tk)
  ka <- which(ks == prof$k_auto)
  graphics::rect(2 * ka - 1.5, graphics::par("usr")[3], 2 * ka + 0.5, graphics::par("usr")[4], border = "#5e2d8c", lty = 2)
  graphics::plot(ks, prof$profile$auc, type = "b", pch = 19, col = "#5e2d8c", lwd = 2, xaxt = "n",
                 ylim = range(c(0.5, prof$profile$auc)), xlab = "kernel size", ylab = "AUC",
                 main = sprintf("Separation from random genes (AUC); automatic kernel size = %d", prof$k_auto), cex.main = 1)
  graphics::axis(1, at = ks)
  graphics::abline(h = 0.5, lty = 3)
  graphics::points(prof$k_auto, max(prof$profile$auc), pch = 8, cex = 2, col = "#5e2d8c")
  invisible(NULL)
}

# ---------------------------------------------------------------------------------------------- one call
# kernel: "auto" or a number. Writes <out>.pdf (profile), <out>.regulators.tsv and <out>.profile.tsv when out is given.
regine <- function(reg, genes, background = "network", kernel = "auto", kernel_range = 2:10, reps = 2000, seed = 1,
                   rank_by = c("score", "p"), name = "my_list", out = file.path("results", name)) {
  prof <- regine_profile(reg, genes, background, kernel_range, reps, seed)
  k <- if (identical(kernel, "auto")) prof$k_auto else as.integer(kernel)
  tab <- regine_regulators(reg, genes, k, background, match.arg(rank_by))
  best <- prof$profile[prof$profile$k == prof$k_auto, ]
  message(sprintf("%d genes vs %d background genes. Automatic kernel size %d: AUC %.3f (%.1f vs %.1f common TFs). Table at k = %d.",
                  prof$n, prof$N, prof$k_auto, best$auc, best$mean_list, best$mean_background, k))
  if (!is.null(out)) {
    dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)
    regine_plot(prof, paste0(out, ".pdf"), title = name)
    utils::write.table(tab, paste0(out, ".regulators.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    utils::write.table(prof$profile, paste0(out, ".profile.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    message("written: ", out, ".pdf, .regulators.tsv, .profile.tsv")
  }
  invisible(list(regulators = tab, profile = prof, k = k))
}
