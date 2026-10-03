// Build the Régine data files from the raw inputs.
//
//   node tools/build_data.mjs --dap DAP.seq.database.Dec.2022.txt --genecloud <dir with genes.json.gz, texts.json.gz>
//                             [--tfdesc TF.descriptions.txt]
//
// Inputs
//   --dap        O'Malley et al. 2016 DAP-seq network as distributed with Régine 2022
//                (space-separated, quoted: rowname "tf.at_id" "target.at_id"; 3.68 M lines)
//   --genecloud  the annotation of the GeneCloud web app (gh-pages branch of GabKrouk/genecloud, folder data/):
//                genes.json.gz {ids, symbols, desc, flags} and texts.json.gz [short, curator summary,
//                computational description, UniProt, phenotypes, type] per gene (TAIR / Araport11 / UniProt, CC BY 4.0)
//   --tfdesc     the TF descriptions of Régine 2022, used when TAIR has no curator summary
//
// Outputs
//   data/dapseq_network.tsv.gz   TF <tab> target, one line per distinct edge (for R)
//   data/tf_annotation.tsv       the 387 TFs: AGI, symbol, family, out-degree, UniProt name, short description, summary, UniProt
//   data/gene_annotation.tsv.gz  every gene: AGI, symbol, short description, protein coding (0/1)
//   web/data/genes.json.gz       {ids, sym, desc, flags}   flags: 1 protein coding, 2 in the DAP-seq network
//   web/data/network.bin.gz      gene-major bit matrix: W uint32 words (little-endian) per gene, bit t = bound by TF t
//   web/data/tfs.json.gz         [{id, sym, fam, deg, name, desc, summary, uniprot}] in the bit order of network.bin
//   web/data/meta.json
import fs from "fs";
import path from "path";
import zlib from "zlib";
import { familyOf } from "./families.mjs";

const arg = k => { const i = process.argv.indexOf(k); return i > 0 ? process.argv[i + 1] : null; };
const DAP = arg("--dap"), GC = arg("--genecloud"), TFDESC = arg("--tfdesc");
if (!DAP || !GC) { console.error("usage: node tools/build_data.mjs --dap <file> --genecloud <dir> [--tfdesc <file>]"); process.exit(1); }
const ROOT = path.resolve(path.dirname(new URL(import.meta.url).pathname), "..");
const out = p => path.join(ROOT, p);
const gz = (p, s) => fs.writeFileSync(out(p), zlib.gzipSync(s, { level: 9 }));
const clean = s => (s || "").replace(/[\t\r\n]+/g, " ").trim();

// ---------------------------------------------------------------- network
const edges = new Map();   // target -> Set(TF)
let nLines = 0;
for (const line of fs.readFileSync(DAP, "utf8").split("\n")) {
  const p = line.replace(/"/g, "").trim().split(/\s+/);
  if (p.length < 3 || p[1] === "tf.at_id") continue;
  nLines++;
  const tf = p[1].toUpperCase(), tg = p[2].toUpperCase();
  (edges.get(tg) ?? edges.set(tg, new Set()).get(tg)).add(tf);
}
const tfs = [...new Set([...edges.values()].flatMap(s => [...s]))].sort();
const nEdges = [...edges.values()].reduce((a, s) => a + s.size, 0);
console.log(`${nLines} lines, ${nEdges} distinct edges (${nLines - nEdges} duplicates), ${tfs.length} TFs, ${edges.size} targets`);

// ---------------------------------------------------------------- annotation
const G = JSON.parse(zlib.gunzipSync(fs.readFileSync(path.join(GC, "genes.json.gz"))));
const T = JSON.parse(zlib.gunzipSync(fs.readFileSync(path.join(GC, "texts.json.gz"))));
const gcIdx = new Map(G.ids.map((g, i) => [g, i]));
const desc2022 = new Map();
if (TFDESC) for (const l of fs.readFileSync(TFDESC, "utf8").split("\n")) {
  const m = l.match(/^"([^"]+)"\s+"[^"]*"\s+"(.*)"\s*$/);
  if (m) desc2022.set(m[1].toUpperCase(), m[2]);
}

const ids = [...new Set([...G.ids, ...edges.keys(), ...tfs])].sort();
const W = Math.ceil(tfs.length / 32), tIdx = new Map(tfs.map((t, i) => [t, i]));
const mask = new Uint32Array(ids.length * W), deg = new Array(tfs.length).fill(0);
const sym = [], desc = [], flags = [];
ids.forEach((g, i) => {
  const k = gcIdx.get(g);
  sym.push(k == null ? "" : G.symbols[k] || "");
  desc.push(k == null ? "" : clean(G.desc[k]));
  const s = edges.get(g);
  flags.push((k != null && G.flags[k] & 1 ? 1 : 0) | (s || tIdx.has(g) ? 2 : 0));
  if (s) for (const tf of s) { const t = tIdx.get(tf); mask[i * W + (t >> 5)] |= 1 << (t & 31); deg[t]++; }
});

const tfRows = tfs.map((t, i) => {
  const k = gcIdx.get(t), tx = (k == null ? null : T[k]) || [];
  return {
    id: t, sym: k == null ? "" : G.symbols[k] || "", fam: familyOf(t, k == null ? "" : G.symbols[k], k == null ? "" : G.desc[k], tx),
    deg: deg[i], desc: clean(tx[0] || (k == null ? "" : G.desc[k])),
    name: clean((tx[3] || "").split(" — ")[0].split(" | ")[0]),
    summary: clean(tx[1] || desc2022.get(t) || ""), uniprot: clean(tx[3] || ""),
  };
});
const unk = tfRows.filter(r => !r.fam).map(r => r.id);
if (unk.length) console.log("no family:", unk.join(" "));

// ---------------------------------------------------------------- write
fs.mkdirSync(out("web/data"), { recursive: true });
let tsv = "tf\ttarget\n";
for (const g of [...edges.keys()].sort()) for (const tf of [...edges.get(g)].sort()) tsv += `${tf}\t${g}\n`;
gz("data/dapseq_network.tsv.gz", tsv);
fs.writeFileSync(out("data/tf_annotation.tsv"),
  "agi\tsymbol\tfamily\tn_targets\tprotein_name\tshort_description\tsummary\tuniprot\n"
  + tfRows.map(r => [r.id, r.sym, r.fam, r.deg, r.name, r.desc, r.summary, r.uniprot].join("\t")).join("\n") + "\n");
gz("data/gene_annotation.tsv.gz", "agi\tsymbol\tshort_description\tprotein_coding\tin_network\n"
  + ids.map((g, i) => [g, sym[i], desc[i], flags[i] & 1, flags[i] & 2 ? 1 : 0].join("\t")).join("\n") + "\n");

gz("web/data/genes.json.gz", JSON.stringify({ ids, sym, desc, flags }));
gz("web/data/network.bin.gz", Buffer.from(mask.buffer));
gz("web/data/tfs.json.gz", JSON.stringify(tfRows));
const meta = {
  regine_web: 1, built: new Date().toISOString().slice(0, 10), W, n_tfs: tfs.length, n_genes: ids.length,
  n_network_genes: flags.filter(f => f & 2).length, n_protein_coding: flags.filter(f => f & 1).length,
  n_edges: nEdges, n_duplicate_lines: nLines - nEdges,
  network: "O'Malley RC et al. (2016) Cell 165:1280-1292, DAP-seq (Régine 2022 table, Dec. 2022)",
  annotation: "TAIR / Araport11 / UniProtKB via GeneCloud (" + (G.built || "2026") + ")",
};
fs.writeFileSync(out("web/data/meta.json"), JSON.stringify(meta, null, 1) + "\n");
console.log(meta);
