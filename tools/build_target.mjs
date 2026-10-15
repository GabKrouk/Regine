// Build the TARGET data files of Régine from the ConnecTF data release.
//
//   node tools/build_target.mjs --connectf <connectf_data_release_v1/arabidopsis> --genecloud <dir with genes.json.gz, texts.json.gz>
//
// Run after tools/build_data.mjs: the gene order of web/data/genes.json.gz and the TF annotation of web/data/tfs.json.gz
// are reused.
//
// Inputs
//   --connectf   the Arabidopsis folder of the ConnecTF data release v1 (Brooks et al. 2021, Plant Physiol 185:49;
//                https://connectf.org): metadata/<experiment>.txt and data/<experiment>.csv. Only the experiments
//                with Edge_Type "Cells:Regulated:TARGET" are used: TF-GR nuclear import with DEX in root or shoot
//                protoplasts, with cycloheximide (direct targets), 180 min, DESeq2 FDR < 0.05 (NLP7: FDR < 0.1).
//                The data files list the significant genes only, with log2FC and adjusted P.
//   --genecloud  the GeneCloud annotation (same as build_data.mjs), for the TFs that are not in the DAP-seq set
//
// Outputs
//   data/target_network.tsv.gz    experiment, TF, target, log2FC, padj, direction (one line per significant gene, for R)
//   data/target_experiments.tsv   the 65 experiments: TF, symbol, family, source, tissue, treatment, cutoff, DOI, counts
//   web/data/target.bin.gz        two gene-major bit matrices (genes in the order of genes.json.gz, W words per gene,
//                                 bit e = experiment e): induced targets, then repressed targets
//   web/data/target.json.gz       {W, n_experiments, n_tfs, n_edges, universe, built, experiments: [...]}
import fs from "fs";
import path from "path";
import zlib from "zlib";
import { familyOf } from "./families.mjs";

const arg = k => { const i = process.argv.indexOf(k); return i > 0 ? process.argv[i + 1] : null; };
const CTF = arg("--connectf"), GC = arg("--genecloud");
if (!CTF || !GC) { console.error("usage: node tools/build_target.mjs --connectf <dir> --genecloud <dir>"); process.exit(1); }
const ROOT = path.resolve(path.dirname(new URL(import.meta.url).pathname), "..");
const out = p => path.join(ROOT, p);
const gunzipJSON = p => JSON.parse(zlib.gunzipSync(fs.readFileSync(p)));
const clean = s => (s || "").replace(/[\t\r\n]+/g, " ").trim();
const AGI = /^AT[1-5CM]G\d{5}$/;
const CONNECTF_DOI = "10.1093/plphys/kiaa012";       // Brooks et al. 2021, Plant Physiol 185:49-66 ("This Work" in the release)

// ---------------------------------------------------------------- experiments
function meta(file) {
  const m = {};
  for (const l of fs.readFileSync(file, "utf8").replace(/^﻿/, "").split(/\r?\n/)) {
    const i = l.indexOf(":");
    if (i > 0) m[l.slice(0, i).replace(/^\*/, "").trim()] = l.slice(i + 1).trim();
  }
  return m;
}
function table(file) {             // first column = gene, then a log2 fold-change column and an adjusted P column
  const lines = fs.readFileSync(file, "utf8").replace(/^﻿/, "").split(/\r?\n/).filter(Boolean);
  const head = lines.shift().split(",").map(s => s.replace(/"/g, "").toLowerCase());
  const fc = head.findIndex(h => h.includes("log2f")), pa = head.findIndex(h => h.includes("padj"));
  if (fc < 0 || pa < 0) throw new Error(`${file}: no log2FC / padj column`);
  return lines.map(l => { const f = l.split(",").map(s => s.replace(/"/g, "").trim()); return { gene: f[0].toUpperCase(), fc: +f[fc], padj: +f[pa] }; });
}
const short = s => s.replace(/^This Work$/, "ConnecTF 2021").replace(/ (Nat Commun|PNAS|Plant Cell|Mol Plant) /, " ");
const tissueOf = s => /root/i.test(s) ? "root" : /shoot/i.test(s) ? "shoot" : s;

const G = gunzipJSON(path.join(GC, "genes.json.gz")), T = gunzipJSON(path.join(GC, "texts.json.gz"));
const gcIdx = new Map(G.ids.map((g, i) => [g, i]));
const genes = gunzipJSON(out("web/data/genes.json.gz")), gIdx = new Map(genes.ids.map((g, i) => [g, i]));
const dapTfs = new Map(gunzipJSON(out("web/data/tfs.json.gz")).map(t => [t.id, t]));

const exps = [];
for (const f of fs.readdirSync(path.join(CTF, "metadata")).sort()) {
  const m = meta(path.join(CTF, "metadata", f));
  if (m.Edge_Type !== "Cells:Regulated:TARGET") continue;
  const name = f.replace(/\.txt$/, ""), tf = m.Transcription_Factor_ID.toUpperCase();
  const rows = table(path.join(CTF, "data", name + ".csv"));
  exps.push({ name, tf, m, rows });
}

function annotate(id) {
  if (dapTfs.has(id)) { const { deg, ...t } = dapTfs.get(id); return { ...t, inDAP: true }; }
  const k = gcIdx.get(id), tx = (k == null ? null : T[k]) || [], sym = k == null ? "" : G.symbols[k] || "";
  return {
    id, sym, fam: familyOf(id, sym, k == null ? "" : G.desc[k], tx),
    desc: clean(tx[0] || (k == null ? "" : G.desc[k])), name: clean((tx[3] || "").split(" — ")[0].split(" | ")[0]),
    summary: clean(tx[1] || ""), uniprot: clean(tx[3] || ""), inDAP: false,
  };
}

// ---------------------------------------------------------------- matrices
const E = exps.length, W = Math.ceil(E / 32), n = genes.ids.length;
const up = new Uint32Array(n * W), dn = new Uint32Array(n * W);
const missing = new Set();
let tsv = "experiment\ttf\ttarget\tlog2FC\tpadj\tdirection\n", nEdges = 0;
const expRows = [];
exps.sort((a, b) => (annotate(a.tf).sym || a.tf).localeCompare(annotate(b.tf).sym || b.tf) || a.name.localeCompare(b.name));
exps.forEach((x, e) => {
  let nu = 0, nd = 0;
  const seen = new Set();
  for (const r of x.rows) {
    if (!AGI.test(r.gene) || seen.has(r.gene)) continue;
    seen.add(r.gene);
    const i = gIdx.get(r.gene);
    if (i == null) { missing.add(r.gene); continue; }
    const m = r.fc > 0 ? up : dn;
    m[i * W + (e >> 5)] |= 1 << (e & 31);
    r.fc > 0 ? nu++ : nd++;
    tsv += `${x.name}\t${x.tf}\t${r.gene}\t${r.fc}\t${r.padj}\t${r.fc > 0 ? "induced" : "repressed"}\n`;
    nEdges++;
  }
  const a = annotate(x.tf), m = x.m;
  if (!a.sym) a.sym = x.name.split("_")[1] || "";          // e.g. AT5G10970_C2H2_DESEQ2 has no TAIR symbol
  expRows.push({
    ...a, key: x.name, n: nu + nd, up: nu, down: nd,
    source: short(m["Data_Source/Publication"] || ""), tissue: tissueOf(m["Tissue/Sample"] || ""),
    treatment: m.Treatments || "", time: m.Treatment_Time || "", cutoff: m.Analysis_Cutoff || "", doi: m.DOI || (m["Data_Source/Publication"] === "This Work" ? CONNECTF_DOI : ""),
  });
});
for (const r of expRows) r.label = `${r.sym || r.id} · ${r.source} · ${r.tissue}`;

let universe = 0;
for (let i = 0; i < n; i++) for (let w = 0; w < W; w++) if (up[i * W + w] | dn[i * W + w]) { universe++; break; }
const tfIds = new Set(expRows.map(r => r.id));
console.log(`${E} TARGET experiments, ${tfIds.size} TFs (${[...tfIds].filter(t => dapTfs.has(t)).length} also in DAP-seq), `
  + `${nEdges} edges, ${universe} genes targeted at least once; ${missing.size} target ids not in the annotation`);
const unk = expRows.filter(r => !r.fam).map(r => r.id);
if (unk.length) console.log("no family:", [...new Set(unk)].join(" "));

// ---------------------------------------------------------------- write
fs.writeFileSync(out("data/target_network.tsv.gz"), zlib.gzipSync(tsv, { level: 9 }));
fs.writeFileSync(out("data/target_experiments.tsv"),
  "experiment\ttf\tsymbol\tfamily\tin_dapseq\tsource\ttissue\ttreatment\ttime\tcutoff\tdoi\tn_targets\tinduced\trepressed\n"
  + expRows.map(r => [r.key, r.id, r.sym, r.fam, r.inDAP ? 1 : 0, r.source, r.tissue, r.treatment, r.time, r.cutoff, r.doi, r.n, r.up, r.down].join("\t")).join("\n") + "\n");
const both = new Uint8Array(2 * n * W * 4);
both.set(new Uint8Array(up.buffer), 0); both.set(new Uint8Array(dn.buffer), n * W * 4);
fs.writeFileSync(out("web/data/target.bin.gz"), zlib.gzipSync(both, { level: 9 }));
fs.writeFileSync(out("web/data/target.json.gz"), zlib.gzipSync(JSON.stringify({
  W, n_experiments: E, n_tfs: tfIds.size, n_edges: nEdges, universe, built: new Date().toISOString().slice(0, 10),
  source: "ConnecTF data release v1 (Brooks et al. 2021, Plant Physiol 185:49-66), TARGET experiments",
  experiments: expRows,
}), { level: 9 }));
