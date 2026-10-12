// Régine on every cell type marker list (engine of web/engine.js; background = expressed genes, 2000 reps, seed 1).
// Usage: node 03_run_regine.mjs <out_dir>   -> regine_all.tsv, regine_profiles.tsv
import fs from "fs";
import { loadRegine } from "../../../tests/load.mjs";
const U = new URL("../data/", import.meta.url).pathname;
const { reg, tfs } = loadRegine();
const bg = fs.readFileSync(U + "background_expressed_genes.txt", "utf8").split(/\s+/).filter(Boolean);
const lines = fs.readFileSync(U + "celltype_markers_top200.tsv", "utf8").trim().split("\n");
const hdr = lines.shift().split("\t");
const lists = {};
for (const l of lines) { const f = l.split("\t"); (lists[f[0]] ??= []).push(f[1]); }
const out = [], prof = [];
for (const [ct, genes] of Object.entries(lists)) {
  const res = reg.run(genes, { background: "custom", custom: bg, reps: 2000, seed: 1 });
  for (const p of res.profile) prof.push([ct, p.k, p.auc.toFixed(4), p.meanList.toFixed(3), p.meanBg.toFixed(3)].join("\t"));
  const rows = [...res.rows].sort((a, b) => a.p - b.p);
  rows.forEach((r, i) => out.push([ct, res.n, res.N, res.kAuto, i + 1, r.tf.id, r.tf.sym, r.tf.fam, r.kt, r.Kt,
    r.fold.toFixed(3), r.p.toExponential(3), r.q.toExponential(3), r.score.toExponential(3),
    reg.targetsOf(res, r.t).join(",")].join("\t")));
  console.error(ct, res.n, "k=", res.kAuto, "auc=", res.profile.find(p => p.k === res.kAuto).auc.toFixed(3), "nFDR<0.05:", rows.filter(r => r.q < 0.05).length);
}
fs.writeFileSync(process.argv[2] + "/regine_all.tsv", "celltype\tn_list\tN_background\tk_auto\trank_by_p\ttf_agi\tsymbol\tfamily\tlist_bound\tbg_bound\tfold\tp\tfdr\tregine_score\tbound_genes\n" + out.join("\n") + "\n");
fs.writeFileSync(process.argv[2] + "/regine_profiles.tsv", "celltype\tk\tauc\tmean_common_list\tmean_common_bg\n" + prof.join("\n") + "\n");
