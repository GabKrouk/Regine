// Régine 2026 — user interface
import { Regine, parseGenes, quantiles } from "./engine.js";

const $ = s => document.querySelector(s), $$ = s => [...document.querySelectorAll(s)];
const esc = s => String(s ?? "").replace(/[&<>"']/g, c => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[c]));
const DATA = "data/";
let REG, GENES, META, RES, VIEW = { k: 0, by: "score", sort: null, dir: -1, all: false, open: new Set() };

// ------------------------------------------------------------------------------------------ theme + banner
const root = document.documentElement;
try { const t = localStorage.getItem("regine-theme"); if (t) root.dataset.theme = t; } catch {}
$("#theme").addEventListener("click", () => {
  const dark = root.dataset.theme ? root.dataset.theme === "dark" : matchMedia("(prefers-color-scheme: dark)").matches;
  root.dataset.theme = dark ? "light" : "dark";
  try { localStorage.setItem("regine-theme", root.dataset.theme); } catch {}
  if (RES) drawCharts();
});
(function banner() {            // a small regulatory network: TFs (bright) wired to targets
  const svg = $("#net"), w = 1400, h = 140, r = (a, b) => a + Math.random() * (b - a);
  svg.setAttribute("viewBox", `0 0 ${w} ${h}`); svg.setAttribute("preserveAspectRatio", "xMidYMid slice");
  const tf = Array.from({ length: 14 }, (_, i) => [80 + i * 95 + r(-20, 20), r(18, 60)]);
  const tg = Array.from({ length: 70 }, () => [r(0, w), r(70, 135)]);
  let s = "";
  for (const t of tg) for (const f of tf) if (Math.random() < 0.12)
    s += `<line x1="${f[0]}" y1="${f[1]}" x2="${t[0]}" y2="${t[1]}" stroke="#c9a8f5" stroke-opacity=".16" stroke-width="1"/>`;
  for (const t of tg) s += `<circle cx="${t[0]}" cy="${t[1]}" r="2.6" fill="#9fb4d9" fill-opacity=".45"/>`;
  for (const f of tf) s += `<circle cx="${f[0]}" cy="${f[1]}" r="5" fill="#d6b8ff" fill-opacity=".75"/>`;
  svg.innerHTML = s;
})();

// ------------------------------------------------------------------------------------------ data
async function fetchBuf(f) {
  const r = await fetch(DATA + f);
  if (!r.ok) throw new Error(`${f}: HTTP ${r.status}`);
  let buf = new Uint8Array(await r.arrayBuffer());
  if (buf[0] === 0x1f && buf[1] === 0x8b) {
    const ds = new Blob([buf]).stream().pipeThrough(new DecompressionStream("gzip"));
    buf = new Uint8Array(await new Response(ds).arrayBuffer());
  }
  return buf;
}
const fetchJSON = async f => JSON.parse(new TextDecoder().decode(await fetchBuf(f)));

async function load() {
  META = await (await fetch(DATA + "meta.json")).json();
  const [genes, tfs, net] = await Promise.all([fetchJSON("genes.json.gz"), fetchJSON("tfs.json.gz"), fetchBuf("network.bin.gz")]);
  GENES = genes;
  const mask = new Uint32Array(net.buffer, net.byteOffset, net.byteLength / 4);
  REG = new Regine({ ids: genes.ids, flags: genes.flags, mask, W: META.W, tfs });
  $("[data-n=network]").textContent = `(${META.n_network_genes.toLocaleString("en")})`;
  $("[data-n=protein_coding]").textContent = `(${META.n_protein_coding.toLocaleString("en")})`;
  const fams = [...new Set(tfs.map(t => t.fam).filter(Boolean))].sort((a, b) => a.localeCompare(b));
  $("#famsel").innerHTML += fams.map(f => `<option>${esc(f)}</option>`).join("");
  $("#annot").textContent = `Network: ${META.n_tfs} TFs, ${META.n_edges.toLocaleString("en")} TF–gene links, ${META.n_network_genes.toLocaleString("en")} genes. `
    + `Annotation: ${META.annotation}; TAIR / Araport11 and UniProtKB are CC BY 4.0. Data built ${META.built}.`;
  $("#status").textContent = "";
  $("#submit").disabled = false;
  fromHash();
}
load().catch(e => { $("#status").innerHTML = `<span class="err">Could not load the data: ${esc(e.message)}</span>`; });

// ------------------------------------------------------------------------------------------ form
function countGenes() {
  const g = parseGenes($("#genes").value);
  $("#gcount").textContent = `${g.length} gene${g.length === 1 ? "" : "s"}`;
  if ($("input[name=bg]:checked").value === "custom") {
    const b = parseGenes($("#bgcustom").value);
    $("#bgcount").textContent = `${b.length} background genes`;
  } else $("#bgcount").textContent = "";
}
$("#genes").addEventListener("input", countGenes);
$("#bgcustom").addEventListener("input", countGenes);
$$("input[name=bg]").forEach(r => r.addEventListener("change", () => {
  $("#bgcustom").style.display = $("input[name=bg]:checked").value === "custom" ? "block" : "none";
  countGenes();
}));
$$(".examples button").forEach(b => b.addEventListener("click", async () => {
  $("#genes").value = await (await fetch(`examples/${b.dataset.ex}.txt`)).text();
  countGenes();
}));
$("#reset").addEventListener("click", () => setTimeout(() => {
  countGenes(); $("#bgcustom").style.display = "none"; $("#results").style.display = "none"; history.replaceState(null, "", location.pathname);
}));

function options() {
  return {
    background: $("input[name=bg]:checked").value,
    custom: parseGenes($("#bgcustom").value),
    kernel: $("input[name=kmode]:checked").value === "auto" ? "auto" : +$("#kfixed").value,
    reps: +$("input[name=reps]:checked").value,
    seed: Math.max(1, +$("#seed").value || 1),
  };
}

$("#form").addEventListener("submit", e => {
  e.preventDefault();
  const genes = parseGenes($("#genes").value), o = options();
  if (genes.length < 2) { $("#status").innerHTML = '<span class="err">Paste at least 2 AGI identifiers.</span>'; return; }
  if (o.background === "custom" && o.custom.length < 50) { $("#status").innerHTML = '<span class="err">Paste your background genes (at least 50).</span>'; return; }
  $("#submit").disabled = true;
  $("#status").innerHTML = '<span class="loader"></span>drawing random kernels…';
  setTimeout(() => {
    try {
      const t0 = performance.now();
      RES = REG.run(genes, o);
      VIEW = { k: RES.k, by: VIEW.by, sort: null, dir: -1, all: false, open: new Set() };
      $("#status").textContent = `done in ${((performance.now() - t0) / 1000).toFixed(1)} s`;
      show();
      writeHash(genes, o);
    } catch (err) {
      $("#status").innerHTML = `<span class="err">${esc(err.message)}</span>`;
    }
    $("#submit").disabled = false;
  }, 30);
});

// ------------------------------------------------------------------------------------------ results
const pct = (x, d = 1) => x === 0 ? "0" : x < 0.001 ? "<0.1%" : (100 * x).toFixed(d) + "%";
const fmtP = p => p >= 0.01 ? p.toFixed(2) : p === 0 ? "<1e-300" : p.toExponential(1).replace("e-", "e−");
const sym = g => GENES.sym[REG.index.get(g)] || "";
const gdesc = g => GENES.desc[REG.index.get(g)] || "";

function show() {
  $("#results").style.display = "block";
  const r = RES, best = r.profile.find(p => p.k === r.kAuto);
  const bgName = { network: "all genes of the DAP-seq network", protein_coding: "protein-coding genes", custom: "your background" }[r.background];
  const strength = best.auc >= 0.8 ? "strongly" : best.auc >= 0.65 ? "clearly" : best.auc >= 0.57 ? "weakly" : "hardly";
  $("#summary").innerHTML = `Your <b>${r.n}</b> genes ${strength} share regulators beyond chance: at the best kernel size,
    <b class="k">k = ${r.kAuto}</b>, a kernel of your list shares more TFs than a kernel of random genes
    (${bgName}, N = ${r.N.toLocaleString("en")}) in <b>${Math.round(best.auc * 100)}%</b> of the draws
    (${best.meanList.toFixed(1)} versus ${best.meanBg.toFixed(1)} common TFs on average).`
    + (r.k !== r.kAuto ? ` The table uses the kernel size you fixed, k = ${r.k}.` : "");
  const w = [];
  if (r.unknown.length) w.push(`${r.unknown.length} identifier${r.unknown.length > 1 ? "s" : ""} not found in Araport11: ${esc(r.unknown.slice(0, 12).join(", "))}${r.unknown.length > 12 ? "…" : ""}.`);
  if (r.notInNetwork.length) w.push(`${r.notInNetwork.length} gene${r.notInNetwork.length > 1 ? "s have" : " has"} no DAP-seq peak (bound by no TF): ${esc(r.notInNetwork.slice(0, 12).join(", "))}${r.notInNetwork.length > 12 ? "…" : ""}.`);
  if (r.added) w.push(`${r.added} gene${r.added > 1 ? "s" : ""} of your list ${r.added > 1 ? "were" : "was"} not in the background and ${r.added > 1 ? "were" : "was"} added to it.`);
  if (best.auc < 0.57) w.push("The separation is close to 0.5 at every kernel size: the candidates below are weak leads.");
  $("#warnings").innerHTML = w.map(s => `<p class="warn">⚠ ${s}</p>`).join("");
  $("#kshow").innerHTML = r.profile.map(p => `<option value="${p.k}"${p.k === VIEW.k ? " selected" : ""}>${p.k}${p.k === r.kAuto ? " (auto)" : ""}</option>`).join("");
  drawCharts();
  table();
  genesTable();
  famTable();
  $("#results").scrollIntoView({ behavior: "smooth", block: "start" });
}

function setK(k) {
  VIEW.k = k;
  $("#kshow").value = k;
  drawCharts();
  table();
  genesTable();
  famTable();
}
$("#kshow").addEventListener("change", e => setK(+e.target.value));

// --------------------------------------------------------------------------- charts (SVG, drawn from CSS tokens)
const css = v => getComputedStyle(root).getPropertyValue(v).trim();
const tip = $("#tip");
function showTip(e, html) { tip.innerHTML = html; tip.style.display = "block"; moveTip(e); }
function moveTip(e) {
  const x = Math.min(e.clientX + 14, innerWidth - tip.offsetWidth - 8), y = e.clientY + 14 + tip.offsetHeight > innerHeight ? e.clientY - tip.offsetHeight - 10 : e.clientY + 14;
  tip.style.left = x + "px"; tip.style.top = y + "px";
}
const hideTip = () => (tip.style.display = "none");

function drawCharts() { profileChart(); aucChart(); }

function profileChart() {
  const P = RES.profile, W = 1060, H = 320, m = { l: 52, r: 12, t: 12, b: 38 };
  const qs = P.map(p => ({ k: p.k, l: quantiles(p.list), b: quantiles(p.bg), p }));
  const ymax = Math.max(1, ...qs.map(q => Math.max(q.l.hi, q.b.hi, q.l.q3, q.b.q3)));
  const f = v => Math.log1p(v), y = v => H - m.b - (f(v) / f(ymax)) * (H - m.t - m.b);
  const cw = (W - m.l - m.r) / P.length, bw = Math.min(26, cw * 0.3);
  const ink = css("--sub"), grid = css("--grid"), cl = css("--list"), cb = css("--bgser"), fl = css("--listfill"), fb = css("--bgfill"), hi = css("--hi");
  let s = `<svg viewBox="0 0 ${W} ${H}" xmlns="http://www.w3.org/2000/svg" font-family="Source Sans 3, Arial, sans-serif" font-size="12">`;
  const ticks = [0, 1, 2, 5, 10, 20, 50, 100, 200, 400].filter(t => t <= ymax * 1.05);
  for (const t of ticks) s += `<line x1="${m.l}" x2="${W - m.r}" y1="${y(t)}" y2="${y(t)}" stroke="${grid}"/><text x="${m.l - 8}" y="${y(t) + 4}" text-anchor="end" fill="${ink}">${t}</text>`;
  s += `<text transform="translate(14 ${(H - m.b + m.t) / 2}) rotate(-90)" text-anchor="middle" fill="${ink}">common TFs per kernel</text>`;
  const box = (x, q, stroke, fill) => {
    const w2 = bw / 2;
    return `<line x1="${x}" x2="${x}" y1="${y(q.hi)}" y2="${y(q.q3)}" stroke="${stroke}" stroke-width="1.5"/>`
      + `<line x1="${x}" x2="${x}" y1="${y(q.q1)}" y2="${y(q.lo)}" stroke="${stroke}" stroke-width="1.5"/>`
      + `<rect x="${x - w2}" y="${y(q.q3)}" width="${bw}" height="${Math.max(1.5, y(q.q1) - y(q.q3))}" rx="3" fill="${fill}" stroke="${stroke}" stroke-width="1.5"/>`
      + `<line x1="${x - w2}" x2="${x + w2}" y1="${y(q.med)}" y2="${y(q.med)}" stroke="${stroke}" stroke-width="2.5"/>`;
  };
  qs.forEach((q, i) => {
    const x0 = m.l + i * cw, on = q.k === VIEW.k;
    s += `<g class="kcol" data-k="${q.k}"><rect class="kbg" x="${x0 + 2}" y="${m.t}" width="${cw - 4}" height="${H - m.t - m.b + 30}" rx="6" fill="${on ? hi : "transparent"}"/>`;
    s += box(x0 + cw / 2 - bw * 0.65, q.l, cl, fl) + box(x0 + cw / 2 + bw * 0.65, q.b, cb, fb);
    s += `<text x="${x0 + cw / 2}" y="${H - m.b + 20}" text-anchor="middle" fill="${on ? css("--fg") : ink}" font-weight="${on ? 700 : 400}">k = ${q.k}${q.k === RES.kAuto ? " ★" : ""}</text></g>`;
  });
  s += "</svg>";
  $("#profile").innerHTML = s;
  $$("#profile .kcol").forEach(g => {
    const q = qs.find(q => q.k === +g.dataset.k);
    g.addEventListener("click", () => setK(q.k));
    g.addEventListener("mousemove", e => showTip(e, `<b>kernel size ${q.k}</b>${q.k === RES.kAuto ? " (automatic)" : ""}<br>`
      + `your list: median ${q.l.med}, mean ${q.p.meanList.toFixed(2)} common TFs<br>random genes: median ${q.b.med}, mean ${q.p.meanBg.toFixed(2)}<br>`
      + `AUC ${q.p.auc.toFixed(3)} · ratio of means ×${(q.p.meanList / Math.max(q.p.meanBg, 1e-9)).toFixed(1)}<br><span class="note">click to use this kernel size</span>`));
    g.addEventListener("mouseleave", hideTip);
  });
}

function aucChart() {
  const P = RES.profile, W = 1060, H = 150, m = { l: 52, r: 12, t: 14, b: 26 };
  const lo = Math.min(0.5, ...P.map(p => p.auc)) - 0.02, hiV = Math.max(...P.map(p => p.auc)) + 0.03;
  const cw = (W - m.l - m.r) / P.length, x = i => m.l + i * cw + cw / 2, y = v => H - m.b - ((v - lo) / (hiV - lo)) * (H - m.t - m.b);
  const ink = css("--sub"), grid = css("--grid"), cl = css("--list"), bg = css("--card");
  let s = `<svg viewBox="0 0 ${W} ${H}" xmlns="http://www.w3.org/2000/svg" font-family="Source Sans 3, Arial, sans-serif" font-size="12">`;
  const step = hiV - lo > 0.3 ? 0.1 : 0.05;
  for (let t = Math.ceil(lo / step) * step; t <= hiV; t += step) s += `<line x1="${m.l}" x2="${W - m.r}" y1="${y(t)}" y2="${y(t)}" stroke="${grid}"/><text x="${m.l - 8}" y="${y(t) + 4}" text-anchor="end" fill="${ink}">${t.toFixed(2)}</text>`;
  s += `<line x1="${m.l}" x2="${W - m.r}" y1="${y(0.5)}" y2="${y(0.5)}" stroke="${ink}" stroke-dasharray="4 4"/><text x="${m.l + 4}" y="${y(0.5) - 5}" fill="${ink}">random (0.5)</text>`;
  s += `<polyline fill="none" stroke="${cl}" stroke-width="2" points="${P.map((p, i) => `${x(i)},${y(p.auc)}`).join(" ")}"/>`;
  P.forEach((p, i) => {
    const on = p.k === VIEW.k, best = p.k === RES.kAuto;
    s += `<g class="kcol" data-k="${p.k}"><rect x="${m.l + i * cw}" y="0" width="${cw}" height="${H}" fill="transparent"/>`
      + `<circle cx="${x(i)}" cy="${y(p.auc)}" r="${on ? 7 : 5}" fill="${best || on ? cl : bg}" stroke="${on ? bg : cl}" stroke-width="2"/>`
      + (best ? `<text x="${x(i)}" y="${y(p.auc) - 12}" text-anchor="middle" fill="${css("--fg")}" font-weight="700">${p.auc.toFixed(3)} ★</text>` : "")
      + `<text x="${x(i)}" y="${H - 6}" text-anchor="middle" fill="${ink}">k = ${p.k}</text></g>`;
  });
  s += "</svg>";
  $("#aucplot").innerHTML = s;
  $$("#aucplot .kcol").forEach(g => {
    const p = P.find(p => p.k === +g.dataset.k);
    g.addEventListener("click", () => setK(p.k));
    g.addEventListener("mousemove", e => showTip(e, `<b>kernel size ${p.k}</b>: AUC ${p.auc.toFixed(3)}${p.k === RES.kAuto ? " (highest: automatic choice)" : ""}`));
    g.addEventListener("mouseleave", hideTip);
  });
}

// --------------------------------------------------------------------------- TF table
$$("#rankby button").forEach(b => b.addEventListener("click", () => {
  $$("#rankby button").forEach(x => x.classList.toggle("on", x === b));
  VIEW.by = b.dataset.r; VIEW.sort = null; table();
}));
["#famsel", "#onlysig"].forEach(s => $(s).addEventListener("change", () => table()));
$("#q").addEventListener("input", () => table());
$("#showall").addEventListener("click", () => { VIEW.all = !VIEW.all; table(); });
$$("#tftab thead th[data-s]").forEach(th => th.addEventListener("click", () => {
  const key = th.dataset.s;
  VIEW.dir = VIEW.sort === key ? -VIEW.dir : (["sym", "fam", "q"].includes(key) ? 1 : -1);
  VIEW.sort = key; table();
}));

function rankedRows() {
  const rows = REG.score(RES, VIEW.k);
  if (VIEW.by === "p") rows.sort((a, b) => a.p - b.p || b.score - a.score);
  rows.forEach((r, i) => (r.rank = i + 1));
  return rows;
}

function table() {
  const rows = rankedRows();
  const fam = $("#famsel").value, sig = $("#onlysig").checked, q = $("#q").value.trim().toLowerCase();
  let shown = rows.filter(r => r.kt > 0 && (!fam || r.tf.fam === fam) && (!sig || r.q <= 0.05)
    && (!q || [r.tf.id, r.tf.sym, r.tf.fam, r.tf.name, r.tf.desc, r.tf.summary].join(" ").toLowerCase().includes(q)));
  if (VIEW.sort) {
    const k = VIEW.sort, get = r => k === "sym" ? (r.tf.sym || r.tf.id) : k === "fam" ? r.tf.fam : r[k];
    shown.sort((a, b) => { const x = get(a), y = get(b); return (typeof x === "string" ? x.localeCompare(y) : x - y) * VIEW.dir; });
  }
  $$("#tftab thead th").forEach(th => th.classList.toggle("sorted", th.dataset.s === VIEW.sort));
  const total = shown.length, limit = VIEW.all ? total : 40;
  shown = shown.slice(0, limit);
  const n = RES.n;
  $("#tftab tbody").innerHTML = shown.map(r => {
    const t = r.tf, open = VIEW.open.has(t.id);
    const row = `<tr class="tfrow${open ? " open" : ""}" data-t="${t.id}"><td class="rank">${r.rank}</td>`
      + `<td class="tf"><b>${esc(t.sym || "—")}</b>${r.q <= 0.05 ? '<span class="sig" title="FDR ≤ 5%">FDR</span>' : ""}<br><a href="https://www.arabidopsis.org/locus?name=${t.id}" target="_blank" rel="noopener" onclick="event.stopPropagation()">${t.id}</a></td>`
      + `<td class="hide-s"><span class="fam">${esc(t.fam || "?")}</span></td>`
      + `<td class="desc hide-s">${tfName(t)}</td>`
      + `<td class="num">${r.kt}/${n}<span class="bar" style="width:${Math.round(40 * r.kt / n)}px"></span></td>`
      + `<td class="num">${pct(r.Kt / RES.N, 0)}<span class="bar barbg" style="width:${Math.round(40 * r.Kt / RES.N)}px"></span></td>`
      + `<td class="num">${pct(r.real)}</td><td class="num">${pct(r.noise)}</td><td class="num"><b>${(100 * r.score).toFixed(1)}</b></td>`
      + `<td class="num">×${r.fold.toFixed(1)}</td><td class="num">${fmtP(r.q)}</td></tr>`;
    return row + (open ? detail(r) : "");
  }).join("") || `<tr><td colspan="11" class="note">No transcription factor matches.</td></tr>`;
  $("#showall").textContent = total > 40 ? (VIEW.all ? "show the first 40 only" : `show all ${total} TFs binding at least one of your genes`) : "";
  $$("#tftab tr.tfrow").forEach(tr => tr.addEventListener("click", () => {
    const id = tr.dataset.t;
    VIEW.open.has(id) ? VIEW.open.delete(id) : VIEW.open.add(id);
    table();
  }));
}

// UniProt name when the Araport11 description is a generic family label
function tfName(t) {
  const generic = /superfamily|family protein|^$/i.test(t.desc);
  if (t.name && generic) return `${esc(t.name)}<br><span class="note">${esc(t.desc)}</span>`;
  return esc(t.desc || t.name);
}

function detail(r) {
  const t = r.tf, bound = new Set(REG.targetsOf(RES, r.t));
  const list = RES.genes.filter(g => REG.index.has(g));
  const chip = g => `<span class="${bound.has(g) ? "" : "off"}" title="${esc(gdesc(g))}">${esc(sym(g) || g)}${sym(g) ? ` <span class="note">${g}</span>` : ""}</span>`;
  const shownGenes = [...list.filter(g => bound.has(g)), ...list.filter(g => !bound.has(g))].slice(0, 300);
  return `<tr class="detail"><td colspan="11"><div class="det">
    <div><h4>${esc(t.sym || t.id)} · ${esc(t.fam || "family unknown")} · ${t.id}</h4>
      <p>${esc(t.summary || t.desc || "No curator summary.")}</p>
      ${t.uniprot ? `<h4 style="margin-top:10px">UniProt</h4><p>${esc(t.uniprot)}</p>` : ""}
      <p class="note" style="margin-top:10px">Binds ${t.deg.toLocaleString("en")} genes genome-wide (${pct(t.deg / META.n_network_genes, 0)} of the network) ·
      ${r.kt} of your ${RES.n} genes · ${r.Kt.toLocaleString("en")} of ${RES.N.toLocaleString("en")} background genes ·
      enrichment ×${r.fold.toFixed(2)}, P = ${fmtP(r.p)}, FDR = ${fmtP(r.q)} ·
      <a href="https://www.arabidopsis.org/locus?name=${t.id}" target="_blank" rel="noopener">TAIR</a> ·
      <a href="https://www.uniprot.org/uniprotkb?query=${t.id}" target="_blank" rel="noopener">UniProt</a></p></div>
    <div><h4>Your genes bound by ${esc(t.sym || t.id)} (${bound.size}/${RES.n}); crossed out: not bound</h4>
      <div class="genes">${shownGenes.map(chip).join("")}${list.length > 300 ? " …" : ""}</div></div>
  </div></td></tr>`;
}

// --------------------------------------------------------------------------- gene and family tables
function genesTable() {
  const top = rankedRows().filter(r => r.kt > 0).slice(0, 20);
  const list = RES.genes.filter(g => REG.index.has(g));
  $("#gtab tbody").innerHTML = list.slice(0, 3000).map(g => {
    const row = REG.index.get(g), w = REG.W, by = top.filter(r => REG.mask[row * w + (r.t >> 5)] & (1 << (r.t & 31)));
    return `<tr><td><a href="https://www.arabidopsis.org/locus?name=${g}" target="_blank" rel="noopener">${g}</a></td><td><b>${esc(sym(g))}</b></td>`
      + `<td class="desc">${esc(gdesc(g))}</td><td class="num">${REG.nTF(row)}</td>`
      + `<td>${by.slice(0, 8).map(r => esc(r.tf.sym || r.tf.id)).join(", ")}${by.length > 8 ? ` +${by.length - 8}` : ""}</td></tr>`;
  }).join("");
}

function famTable() {
  const rows = rankedRows(), by = new Map();
  rows.forEach((r, i) => {
    const f = r.tf.fam || "unknown";
    const e = by.get(f) || by.set(f, { n: 0, sig: 0, top: 0, best: [] }).get(f);
    e.n++; if (r.q <= 0.05) e.sig++; if (i < 20) e.top++;
    if (e.best.length < 4 && r.kt > 0) e.best.push(r.tf.sym || r.tf.id);
  });
  $("#ftab tbody").innerHTML = [...by.entries()].sort((a, b) => b[1].top - a[1].top || b[1].sig - a[1].sig || b[1].n - a[1].n)
    .map(([f, e]) => `<tr><td><span class="fam">${esc(f)}</span></td><td class="num">${e.n}</td><td class="num">${e.sig}</td><td class="num">${e.top || ""}</td><td>${esc(e.best.join(", "))}</td></tr>`).join("");
}

$$(".tabs button").forEach(b => b.addEventListener("click", () => {
  $$(".tabs button").forEach(x => x.classList.toggle("on", x === b));
  $$(".panel").forEach(p => p.classList.toggle("on", p.id === b.dataset.p));
}));

// --------------------------------------------------------------------------- export
function download(name, data, type) {
  const a = document.createElement("a");
  a.href = URL.createObjectURL(new Blob([data], { type }));
  a.download = name; a.click();
  setTimeout(() => URL.revokeObjectURL(a.href), 2000);
}
$("#dltsv").addEventListener("click", () => {
  const rows = rankedRows();
  const head = ["rank", "agi", "symbol", "family", "your_genes_bound", "n_list", "background_bound", "N_background", "real", "noise", "regine_score",
    "fold", "p", "fdr", "tf_targets_genome", "protein_name", "short_description", "summary", "your_genes_bound_ids"];
  const body = rows.map(r => [r.rank, r.tf.id, r.tf.sym, r.tf.fam, r.kt, RES.n, r.Kt, RES.N, r.real, r.noise, r.score, r.fold, r.p, r.q, r.tf.deg,
    r.tf.name, r.tf.desc, r.tf.summary, REG.targetsOf(RES, r.t).join(",")].join("\t"));
  const hdr = `# Régine 2026 · kernel size ${VIEW.k}${VIEW.k === RES.kAuto ? " (automatic)" : ""} · ranked by ${VIEW.by === "p" ? "enrichment P" : "Régine score"} · background ${RES.background} (N=${RES.N})\n`;
  download(`regine_k${VIEW.k}.tsv`, hdr + head.join("\t") + "\n" + body.join("\n") + "\n", "text/tab-separated-values");
});
$("#dlprof").addEventListener("click", () => {
  const head = "kernel_size\tauc\texact_mean_list\texact_mean_background\tlist_median\tbackground_median\tlist_q1\tlist_q3\tbackground_q1\tbackground_q3\tdraws";
  const body = RES.profile.map(p => { const l = quantiles(p.list), b = quantiles(p.bg); return [p.k, p.auc, p.meanList, p.meanBg, l.med, b.med, l.q1, l.q3, b.q1, b.q3, RES.reps].join("\t"); });
  download("regine_kernel_profile.tsv", head + "\n" + body.join("\n") + "\n", "text/tab-separated-values");
});
$("#dlsvg").addEventListener("click", () => {
  const svg = $("#profile svg").outerHTML.replace("<svg ", `<svg style="background:${css("--card")}" `);
  download("regine_kernel_profile.svg", svg, "image/svg+xml");
});
$("#share").addEventListener("click", async () => {
  try { await navigator.clipboard.writeText(location.href); $("#share").textContent = "Link copied"; }
  catch { $("#share").textContent = "Copy the address bar"; }
  setTimeout(() => ($("#share").textContent = "Copy link to this analysis"), 2000);
});

function writeHash(genes, o) {
  if (o.background === "custom") { history.replaceState(null, "", location.pathname); return; }   // too long for a URL
  const p = new URLSearchParams({ g: genes.join(","), bg: o.background, k: o.kernel, r: o.reps, s: o.seed });
  history.replaceState(null, "", "#" + p.toString());
}
function fromHash() {
  if (!location.hash.includes("g=")) return;
  const p = new URLSearchParams(location.hash.slice(1));
  $("#genes").value = (p.get("g") || "").split(",").join("\n");
  const bg = p.get("bg"); if (bg) { const r = $(`input[name=bg][value="${bg}"]`); if (r) r.checked = true; }
  const k = p.get("k");
  if (k && k !== "auto") { $("input[name=kmode][value=fixed]").checked = true; $("#kfixed").value = k; }
  const r = p.get("r"); if (r) { const x = $(`input[name=reps][value="${r}"]`); if (x) x.checked = true; }
  if (p.get("s")) $("#seed").value = p.get("s");
  countGenes();
  $("#form").requestSubmit();
}
