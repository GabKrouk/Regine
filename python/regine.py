#!/usr/bin/env python3
"""Régine 2026 (Python) — candidate transcription factors regulating an Arabidopsis gene list.

Two networks, the same as the web version (https://gabkrouk.github.io/Regine/):
  * DAP-seq (O'Malley et al. 2016, Cell 165:1280): genes bound in vitro by 387 TFs  -> data/dapseq_network.tsv.gz
  * TARGET (ConnecTF, Brooks et al. 2021, Plant Physiol 185:49): genes induced or repressed by 57 TFs in protoplasts,
    with cycloheximide (direct targets), 65 experiments                               -> data/target_network.tsv.gz

Python 3.8+, standard library only (matplotlib is used for the profile plot if it is installed).

Command line, from the repository root:

    python3 python/regine.py data/example_nitrate_list.txt --name nitrate
    python3 python/regine.py my_list.txt --background expressed_genes.txt --network target --kernel 4

writes results/<name>.dapseq.tsv, results/<name>.target.tsv, results/<name>.both.tsv and results/<name>.profile.tsv
(and results/<name>.profile.pdf with matplotlib). In Python:

    import regine
    data = regine.load("data")
    res = regine.run(data, regine.read_genes("data/example_nitrate_list.txt"))
    res["both"][:5]

For a list L of n genes, a background B of N genes (L is added to B), and each column t of a network (a DAP-seq TF or a
TARGET experiment) with k_t targets in L and K_t in B:
    real  = C(k_t, k) / C(n, k)   fraction of the kernels of k genes of L entirely bound / regulated by t
    noise = C(K_t, k) / C(N, k)   same for random kernels of background genes
    score = real - noise
    p     = P(X >= k_t), X ~ hypergeometric(N, K_t, n);  fdr = Benjamini-Hochberg over the columns of the network
Kernel profile: for k = 2..10, `reps` random kernels of L and of B; number of columns covering every gene of the kernel;
AUC = P(list kernel > random kernel) + ties / 2. The automatic k is the largest AUC.
TARGET is also tested on the induced and on the repressed targets of each experiment. Side by side ("both"): one row per
TF, DAP-seq P and the best TARGET experiment's P times its number of experiments, combined by Fisher's method, FDR over TFs.
"""
import argparse
import csv
import gzip
import math
import os
import random
import re
import sys

AGI = re.compile(r"AT[1-5CM]G\d{5}", re.I)


# ------------------------------------------------------------------------------------------------- data
def _tsv(path):
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt", encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t", quoting=csv.QUOTE_NONE))


def read_genes(x):
    """AGIs from a file name, free text or a list; upper case, unique, in order; isoform suffixes and # lines ignored."""
    if isinstance(x, str) and os.path.exists(x):
        with open(x, encoding="utf-8") as f:
            x = f.read()
    if not isinstance(x, str):
        x = "\n".join(x)
    x = "\n".join(l for l in x.splitlines() if not l.lstrip().startswith("#"))
    return list(dict.fromkeys(g.upper() for g in AGI.findall(x)))


class Network:
    """Columns (TFs or experiments) -> sets of target gene ids, plus per-gene bit masks for the kernels."""

    def __init__(self, name, cols, targets, universe):
        self.name, self.cols, self.targets, self.universe = name, cols, targets, universe
        self.mask = {}
        for j, s in enumerate(targets):
            b = 1 << j
            for g in s:
                self.mask[g] = self.mask.get(g, 0) | b


def load(datadir="data"):
    """Read both networks and the annotation (a few seconds)."""
    genes = _tsv(os.path.join(datadir, "gene_annotation.tsv.gz"))
    tfs = _tsv(os.path.join(datadir, "tf_annotation.tsv"))
    ids = {g["agi"] for g in genes}
    sym = {g["agi"]: g["symbol"] for g in genes}
    protein_coding = {g["agi"] for g in genes if g["protein_coding"] == "1"}

    tcol = {t["agi"]: j for j, t in enumerate(tfs)}
    dsets = [set() for _ in tfs]
    with gzip.open(os.path.join(datadir, "dapseq_network.tsv.gz"), "rt") as f:
        next(f)
        for line in f:
            tf, tg = line.rstrip("\n").split("\t")
            dsets[tcol[tf]].add(tg)
    dap_universe = {g["agi"] for g in genes if g["in_network"] == "1"}
    dap = Network("DAP-seq", [dict(t, key=t["agi"]) for t in tfs], dsets, dap_universe)

    exps = _tsv(os.path.join(datadir, "target_experiments.tsv"))
    ecol = {e["experiment"]: j for j, e in enumerate(exps)}
    up, dn = [set() for _ in exps], [set() for _ in exps]
    with gzip.open(os.path.join(datadir, "target_network.tsv.gz"), "rt") as f:
        next(f)
        for line in f:
            e, tf, tg, fc, padj, d = line.rstrip("\n").split("\t")
            if tg in ids:
                (up if d == "induced" else dn)[ecol[e]].add(tg)
    anyset = [u | d for u, d in zip(up, dn)]
    t_universe = set().union(*anyset)
    tann = {t["agi"]: t for t in tfs}
    cols = []
    for e in exps:
        a = tann.get(e["tf"], {})
        cols.append(dict(e, key=e["experiment"], agi=e["tf"], label="%s · %s · %s" % (e["symbol"], e["source"], e["tissue"]),
                         short_description=a.get("short_description", ""), summary=a.get("summary", "")))
    shared = [j for j, e in enumerate(exps) if e["tf"] in tcol]
    inter = [dsets[tcol[exps[j]["tf"]]] & anyset[j] for j in shared]
    return {
        "ids": ids, "symbol": sym, "protein_coding": protein_coding, "tfs": tfs,
        "dap": dap,
        "target": Network("TARGET", cols, anyset, t_universe),
        "target_up": Network("TARGET induced", cols, up, t_universe),
        "target_down": Network("TARGET repressed", cols, dn, t_universe),
        "inter": Network("bound and regulated", [cols[j] for j in shared], inter, dap_universe & t_universe),
    }


# ------------------------------------------------------------------------------------------------- statistics
def _lchoose(n, k):
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def hyper_upper(x, N, K, n):
    """P(X >= x), X ~ hypergeometric(N population, K successes, n draws)."""
    lo, hi = max(x, 0, n - (N - K)), min(n, K)
    if lo > hi:
        return 0.0
    base = _lchoose(N, n)
    terms = [_lchoose(K, i) + _lchoose(N - K, n - i) - base for i in range(lo, hi + 1)]
    m = max(terms)
    return min(1.0, math.exp(m) * sum(math.exp(t - m) for t in terms))


def bh(p):
    n, q, run = len(p), [0.0] * len(p), 1.0
    for r, i in enumerate(sorted(range(n), key=lambda i: -p[i])):
        run = min(run, p[i] * n / (n - r))
        q[i] = run
    return q


def choose_ratio(a, b, k):
    """C(a, k) / C(b, k)."""
    if a < k:
        return 0.0
    r = 1.0
    for i in range(k):
        r *= (a - i) / (b - i)
    return r


def fisher2(p1, p2):
    """Fisher's method for two P values (chi-square, 4 degrees of freedom)."""
    x = -2 * (math.log(max(p1, 1e-300)) + math.log(max(p2, 1e-300)))
    return min(1.0, math.exp(-x / 2) * (1 + x / 2))


def auc(x, y):
    """P(X > Y) + P(X = Y) / 2."""
    from collections import Counter
    cy, s, below = Counter(y), 0.0, {}
    acc = 0
    for v in sorted(set(x) | set(y)):
        below[v] = acc
        acc += cy.get(v, 0)
    for v in x:
        s += below[v] + cy.get(v, 0) / 2
    return s / (len(x) * len(y))


# ------------------------------------------------------------------------------------------------- one network
def background(data, net, kind):
    if kind == "network":
        return set(net.universe)
    if kind == "protein_coding":
        return set(data["protein_coding"])
    return set(read_genes(kind)) & data["ids"]          # a file of background genes


def enrich(net, L, B, k):
    """One row per column: k_t, K_t, real, noise, score, fold, p, fdr."""
    n, N, Ls = len(L), len(B), set(L)
    rows = []
    for j, s in enumerate(net.targets):
        kt, Kt = len(s & Ls), len(s & B)
        real, noise = choose_ratio(kt, n, k), choose_ratio(Kt, N, k)
        rows.append(dict(col=j, kt=kt, Kt=Kt, real=real, noise=noise, score=real - noise,
                         fold=(kt / n) / (Kt / N) if Kt else float("inf"), p=hyper_upper(kt, N, Kt, n)))
    for r, q in zip(rows, bh([r["p"] for r in rows])):
        r["fdr"] = q
    return rows


def profile(net, L, B, kmin=2, kmax=10, reps=2000, seed=1):
    rnd, Bl = random.Random(seed), sorted(B)
    mask = net.mask
    cl = [len(s & set(L)) for s in net.targets]
    cb = [len(s & B) for s in net.targets]

    def common(genes):
        m = -1
        for g in genes:
            m &= mask.get(g, 0)
        return bin(m).count("1") if m != -1 else 0

    out = []
    for k in range(kmin, min(kmax, len(L)) + 1):
        xl = [common(rnd.sample(L, k)) for _ in range(reps)]
        xb = [common(rnd.sample(Bl, k)) for _ in range(reps)]
        out.append(dict(k=k, auc=auc(xl, xb),
                        mean_list=sum(choose_ratio(c, len(L), k) for c in cl),
                        mean_background=sum(choose_ratio(c, len(B), k) for c in cb),
                        median_list=sorted(xl)[reps // 2], median_background=sorted(xb)[reps // 2]))
    return out


# ------------------------------------------------------------------------------------------------- full analysis
def run(data, genes, bg="network", kernel="auto", reps=2000, seed=1, networks=("dap", "target")):
    """Both networks; returns {'dap': rows, 'target': rows, 'both': rows, 'profiles': {...}, 'k': {...}, ...}."""
    genes = read_genes(genes)
    unknown = [g for g in genes if g not in data["ids"]]
    L = [g for g in genes if g in data["ids"]]
    if len(L) < 2:
        raise ValueError("at least 2 recognised genes are needed")
    res = {"genes": L, "unknown": unknown, "profiles": {}, "k": {}, "N": {}}
    for key in networks:
        net = data[key]
        B = background(data, net, bg) | set(L)
        prof = profile(net, L, B, reps=reps, seed=seed)
        k_auto = max(prof, key=lambda p: p["auc"])["k"]
        k = k_auto if kernel == "auto" else min(int(kernel), len(L))
        rows = enrich(net, L, B, k)
        for r in rows:
            c = net.cols[r["col"]]
            r.update(agi=c["agi"], symbol=c["symbol"], family=c["family"], your_genes=",".join(g for g in L if g in net.targets[r["col"]]))
        if key == "target":
            for name in ("up", "down"):
                sub = enrich(data["target_" + name], L, B, k)
                for r, s in zip(rows, sub):
                    r["k_" + name], r["fdr_" + name], r["p_" + name] = s["kt"], s["fdr"], s["p"]
            for r in rows:
                c = net.cols[r["col"]]
                u, d = r["fdr_up"] <= 0.05, r["fdr_down"] <= 0.05
                r.update(experiment=c["experiment"], source=c["source"], tissue=c["tissue"],
                         effect="both" if u and d else "activates" if u else "represses" if d else "")
            Bi = background(data, data["inter"], bg) | set(L)
            res["inter"] = enrich(data["inter"], L, Bi, 2)
        res[key], res["profiles"][key], res["k"][key], res["N"][key] = rows, prof, k, len(B)
        res.setdefault("k_auto", {})[key] = k_auto
    if "dap" in res and "target" in res:
        res["both"] = side_by_side(data, res)
    return res


def side_by_side(data, res):
    by = {}
    for r in res["dap"]:
        by.setdefault(r["agi"], {"agi": r["agi"], "symbol": r["symbol"], "family": r["family"]})["dap"] = r
    for r in res["target"]:
        by.setdefault(r["agi"], {"agi": r["agi"], "symbol": r["symbol"], "family": r["family"]}).setdefault("exps", []).append(r)
    for j, r in enumerate(res["inter"]):
        by[data["inter"].cols[j]["agi"]].setdefault("inter", []).append(r)
    rows = list(by.values())
    for r in rows:
        d, e, i = r.get("dap"), r.get("exps", []), r.get("inter", [])
        best = min(e, key=lambda x: x["p"]) if e else None
        bi = min(i, key=lambda x: x["p"]) if i else None
        r.update(dap_genes=d["kt"] if d else "", dap_p=d["p"] if d else "", dap_fdr=d["fdr"] if d else "",
                 target_best=best["experiment"] if best else "", target_n_experiments=len(e),
                 target_genes=best["kt"] if best else "", target_p=min(1.0, best["p"] * len(e)) if best else "",
                 target_fdr=best["fdr"] if best else "", target_effect=best["effect"] if best else "",
                 bound_and_regulated=bi["kt"] if bi else "", bound_and_regulated_p=min(1.0, bi["p"] * len(i)) if bi else "",
                 bound_and_regulated_fdr=bi["fdr"] if bi else "")
        r["combined_p"] = fisher2(r["dap_p"], r["target_p"]) if d and best else (r["dap_p"] if d else r["target_p"])
    for r, q in zip(rows, bh([r["combined_p"] for r in rows])):
        r["combined_fdr"] = q
    rows.sort(key=lambda r: (r["combined_p"], r["symbol"]))
    return rows


# ------------------------------------------------------------------------------------------------- output
def _write(path, rows, cols):
    with open(path, "w", encoding="utf-8") as f:
        f.write("\t".join(["rank"] + cols) + "\n")
        for i, r in enumerate(rows, 1):
            f.write("\t".join([str(i)] + ["%.4g" % r[c] if isinstance(r[c], float) else str(r[c]) for c in cols]) + "\n")


def write_results(res, out):
    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)
    written = []
    if "dap" in res:
        _write(out + ".dapseq.tsv", sorted(res["dap"], key=lambda r: (r["p"], -r["score"])),
               ["agi", "symbol", "family", "kt", "Kt", "real", "noise", "score", "fold", "p", "fdr", "your_genes"])
        written.append(out + ".dapseq.tsv")
    if "target" in res:
        _write(out + ".target.tsv", sorted(res["target"], key=lambda r: (r["p"], -r["score"])),
               ["experiment", "agi", "symbol", "family", "source", "tissue", "kt", "Kt", "real", "noise", "score", "fold", "p", "fdr",
                "k_up", "fdr_up", "k_down", "fdr_down", "effect", "your_genes"])
        written.append(out + ".target.tsv")
    if "both" in res:
        _write(out + ".both.tsv", res["both"], ["agi", "symbol", "family", "dap_genes", "dap_fdr", "target_best", "target_n_experiments",
               "target_genes", "target_fdr", "target_effect", "bound_and_regulated", "bound_and_regulated_fdr", "combined_p", "combined_fdr"])
        written.append(out + ".both.tsv")
    with open(out + ".profile.tsv", "w") as f:
        f.write("network\tk\tauc\tmean_list\tmean_background\tmedian_list\tmedian_background\n")
        for net, prof in res["profiles"].items():
            for p in prof:
                f.write("\t".join([net] + [str(round(v, 4)) if isinstance(v, float) else str(v) for v in p.values()]) + "\n")
    written.append(out + ".profile.tsv")
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(7, 3.2))
        for net, color in (("dap", "#7b3fb0"), ("target", "#17806f")):
            if net in res["profiles"]:
                P = res["profiles"][net]
                ax.plot([p["k"] for p in P], [p["auc"] for p in P], "o-", color=color, label="%s (k = %d)" % ("DAP-seq" if net == "dap" else "TARGET", res["k_auto"][net]))
        ax.axhline(0.5, ls=":", color="grey")
        ax.set_xlabel("kernel size")
        ax.set_ylabel("AUC")
        ax.legend(frameon=False)
        fig.tight_layout()
        fig.savefig(out + ".profile.pdf")
        written.append(out + ".profile.pdf")
    except ImportError:
        pass
    return written


def main(argv=None):
    a = argparse.ArgumentParser(description="Régine 2026: candidate regulators of an Arabidopsis gene list (DAP-seq and TARGET).")
    a.add_argument("genes", help="file with AGI identifiers (any separator)")
    a.add_argument("--background", default="network", help="network (default), protein_coding, or a file of background genes")
    a.add_argument("--kernel", default="auto", help="auto (largest AUC) or a kernel size")
    a.add_argument("--network", default="both", choices=["both", "dap", "target"])
    a.add_argument("--reps", type=int, default=2000, help="random kernels per kernel size (profile)")
    a.add_argument("--seed", type=int, default=1)
    a.add_argument("--data", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data"))
    a.add_argument("--name", default="my_list")
    a.add_argument("--out", default=None, help="output prefix (default results/<name>)")
    o = a.parse_args(argv)
    data = load(o.data)
    nets = {"both": ("dap", "target"), "dap": ("dap",), "target": ("target",)}[o.network]
    res = run(data, o.genes, o.background, o.kernel, o.reps, o.seed, nets)
    if res["unknown"]:
        print("%d identifier(s) not in Araport11: %s" % (len(res["unknown"]), ", ".join(res["unknown"][:10])), file=sys.stderr)
    for net in nets:
        best = max(res["profiles"][net], key=lambda p: p["auc"])
        top = sorted(res[net], key=lambda r: r["p"])[:5]
        print("%s: %d genes vs %d background genes, automatic kernel size %d (AUC %.3f); top: %s" % (
            "DAP-seq" if net == "dap" else "TARGET", len(res["genes"]), res["N"][net], best["k"], best["auc"],
            ", ".join("%s %d (FDR %.1e)" % (r["symbol"], r["kt"], r["fdr"]) for r in top)))
    print("written:", ", ".join(write_results(res, o.out or os.path.join("results", o.name))))


if __name__ == "__main__":
    main()
