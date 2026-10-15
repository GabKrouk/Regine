"""Régine 2.0 on every cell-type marker list: DAP-seq and TARGET side by side (python/regine.py).
Background = expressed genes, 2000 kernels per size, seed 1.
Usage: python3 06_run_target.py <out_dir>  -> target_all_celltypes.tsv.gz, both_all_celltypes.tsv.gz, target_profiles.tsv"""
import csv
import gzip
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "python"))
import regine  # noqa: E402

D = os.path.join(HERE, "..", "data")
out = sys.argv[1]
data = regine.load(os.path.join(HERE, "..", "..", "..", "data"))
bg = os.path.join(D, "background_expressed_genes.txt")
lists = {}
with open(os.path.join(D, "celltype_markers_top200.tsv")) as f:
    for r in csv.DictReader(f, delimiter="\t"):
        lists.setdefault(r["celltype"], []).append(r["gene"])

T = ["celltype", "n_list", "N_background", "k_auto", "auc", "rank", "experiment", "agi", "symbol", "family", "source", "tissue",
     "kt", "Kt", "fold", "p", "fdr", "k_up", "fdr_up", "k_down", "fdr_down", "effect"]
B = ["celltype", "rank", "agi", "symbol", "family", "dap_genes", "dap_fdr", "target_best", "target_genes", "target_fdr", "target_effect",
     "bound_and_regulated", "bound_and_regulated_fdr", "combined_fdr"]
fmt = lambda v: "%.4g" % v if isinstance(v, float) else str(v)
with gzip.open(os.path.join(out, "target_all_celltypes.tsv.gz"), "wt") as ft, \
     gzip.open(os.path.join(out, "both_all_celltypes.tsv.gz"), "wt") as fb, \
     open(os.path.join(out, "target_profiles.tsv"), "w") as fp:
    ft.write("\t".join(T) + "\n"); fb.write("\t".join(B) + "\n")
    fp.write("celltype\tnetwork\tk\tauc\tmean_list\tmean_background\n")
    for ct, genes in lists.items():
        res = regine.run(data, genes, bg=bg, reps=2000, seed=1)
        for net in ("dap", "target"):
            for p in res["profiles"][net]:
                fp.write("\t".join([ct, net, str(p["k"]), "%.4f" % p["auc"], "%.3f" % p["mean_list"], "%.3f" % p["mean_background"]]) + "\n")
        auc = max(p["auc"] for p in res["profiles"]["target"])
        rows = sorted(res["target"], key=lambda r: min(r["p"], r["p_up"], r["p_down"]))
        for i, r in enumerate(rows, 1):
            ft.write("\t".join(fmt(x) for x in [ct, len(res["genes"]), res["N"]["target"], res["k_auto"]["target"], auc, i] + [r[c] for c in T[6:]]) + "\n")
        for i, r in enumerate(res["both"], 1):
            fb.write("\t".join(fmt(x) for x in [ct, i] + [r[c] for c in B[2:]]) + "\n")
        top = [r for r in rows if min(r["fdr"], r["fdr_up"], r["fdr_down"]) <= 0.05]
        print("%-28s n=%3d  AUC DAP %.3f TARGET %.3f  TARGET exps FDR<=5%%: %2d  top: %s" % (
            ct, len(res["genes"]), max(p["auc"] for p in res["profiles"]["dap"]), auc, len(top),
            ", ".join("%s(%s)%s" % (r["symbol"], r["tissue"][0], {"activates": "↑", "represses": "↓", "both": "↑↓"}.get(r["effect"], "")) for r in top[:6])))
        print("%28s combined top: %s" % ("", ", ".join("%s%s" % (r["symbol"], "*" if r["dap_genes"] != "" and r["target_genes"] != "" else "") for r in res["both"][:8])))
