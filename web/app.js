// Régine 2026 — user interface
import { Regine, parseGenes, quantiles, targetNetworks, fisher2, bh } from "./engine.js";

const $ = s => document.querySelector(s), $$ = s => [...document.querySelectorAll(s)];
const esc = s => String(s ?? "").replace(/[&<>"']/g, c => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[c]));
const DATA = "data/";
// networks: dap = DAP-seq (TFs), any = TARGET (experiments), both = side by side (one row per TF)
const NETNAME = { dap: "DAP-seq", any: "TARGET", both: "DAP-seq + TARGET" };
let DAP, T, GENES, META, TMETA, R;
let VIEW = { net: "both", k: {}, by: { dap: "p", any: "p", both: "c" }, sort: null, dir: -1, all: false, open: new Set() };

// ------------------------------------------------------------------------------------------ theme + banner
const root = document.documentElement;
try { const t = localStorage.getItem("regine-theme-2026"); if (t) root.dataset.theme = t; } catch {}
$("#theme").addEventListener("click", () => {
  const dark = root.dataset.theme ? root.dataset.theme === "dark" : matchMedia("(prefers-color-scheme: dark)").matches;
  root.dataset.theme = dark ? "light" : "dark";
  try { localStorage.setItem("regine-theme-2026", root.dataset.theme); } catch {}
  if (R) drawCharts();
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
  const [genes, tfs, net, tmeta, tbuf] = await Promise.all([fetchJSON("genes.json.gz"), fetchJSON("tfs.json.gz"), fetchBuf("network.bin.gz"),
    fetchJSON("target.json.gz"), fetchBuf("target.bin.gz")]);
  GENES = genes; TMETA = tmeta;
  const mask = new Uint32Array(net.buffer, net.byteOffset, net.byteLength / 4);
  DAP = new Regine({ ids: genes.ids, flags: genes.flags, mask, W: META.W, tfs });
  T = targetNetworks(DAP, tbuf, tmeta);
  const fmt = x => x.toLocaleString("en");
  $("[data-n=network]").textContent = `(DAP-seq ${fmt(META.n_network_genes)} · TARGET ${fmt(tmeta.universe)})`;
  $("[data-n=protein_coding]").textContent = `(${fmt(META.n_protein_coding)})`;
  const fams = [...new Set([...tfs, ...tmeta.experiments].map(t => t.fam).filter(Boolean))].sort((a, b) => a.localeCompare(b));
  $("#famsel").innerHTML += fams.map(f => `<option>${esc(f)}</option>`).join("");
  $("#annot").textContent = `DAP-seq: ${META.n_tfs} TFs, ${fmt(META.n_edges)} TF–gene links, ${fmt(META.n_network_genes)} genes. `
    + `TARGET: ${tmeta.n_experiments} experiments, ${tmeta.n_tfs} TFs, ${fmt(tmeta.n_edges)} regulated TF–gene pairs, ${fmt(tmeta.universe)} genes. `
    + `Annotation: ${META.annotation}; TAIR / Araport11 and UniProtKB are CC BY 4.0. Data built ${META.built} (DAP-seq), ${tmeta.built} (TARGET).`;
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
      R = analyse(genes, o);
      VIEW = { net: VIEW.net, k: { dap: R.dap.k, any: R.any.k }, by: VIEW.by, sort: null, dir: -1, all: false, open: new Set() };
      $("#status").textContent = `done in ${((performance.now() - t0) / 1000).toFixed(1)} s`;
      show();
      writeHash(genes, o);
    } catch (err) {
      $("#status").innerHTML = `<span class="err">${esc(err.message)}</span>`;
      console.error(err);
    }
    $("#submit").disabled = false;
  }, 30);
});

// both networks; TARGET also with the induced and the repressed targets, and bound (DAP-seq) ∩ regulated (TARGET)
function analyse(genes, o) {
  const dap = DAP.run(genes, o), any = T.any.run(genes, o);
  const fixed = { ...o, profile: false, kernel: any.k };
  return { genes, o, dap, any, up: T.up.run(genes, fixed), down: T.down.run(genes, fixed), inter: T.inter.run(genes, { ...o, profile: false, kernel: 2 }) };
}

// ------------------------------------------------------------------------------------------ results
const pct = (x, d = 1) => x === 0 ? "0" : x < 0.001 ? "<0.1%" : (100 * x).toFixed(d) + "%";
const fmtP = p => p == null ? "—" : p >= 0.01 ? p.toFixed(2) : p === 0 ? "<1e-300" : p.toExponential(1).replace("e-", "e−");
const sym = g => GENES.sym[DAP.index.get(g)] || "";
const gdesc = g => GENES.desc[DAP.index.get(g)] || "";
const cur = () => VIEW.net === "dap" ? { reg: DAP, res: R.dap } : { reg: T.any, res: R.any };
const bgName = (bg, net) => ({ network: net === "dap" ? "all genes of the DAP-seq network" : "all genes regulated in at least one TARGET experiment",
  protein_coding: "protein-coding genes", custom: "your background" }[bg]);
const strength = a => a >= 0.8 ? "strongly" : a >= 0.65 ? "clearly" : a >= 0.57 ? "weakly" : "hardly";
const effect = (qu, qd) => qu <= 0.05 && qd <= 0.05 ? "both" : qu <= 0.05 ? "activates" : qd <= 0.05 ? "represses" : "";
const effectHtml = e => e ? `<span class="eff ${e}">${{ activates: "↑ activates", represses: "↓ represses", both: "↑↓ both" }[e]}</span>` : "";

function show() {
  $("#results").style.display = "block";
  $$("#netsel button").forEach(b => b.classList.toggle("on", b.dataset.n === VIEW.net));
  const one = (net, res) => {
    const best = res.profile.find(p => p.k === res.kAuto), unit = net === "dap" ? "TFs" : "TARGET experiments";
    return `<b>${NETNAME[net]}:</b> your <b>${res.n}</b> genes ${strength(best.auc)} share regulators beyond chance. At the best kernel size,
      <b class="k">k = ${res.kAuto}</b>, a kernel of your list shares more ${unit} than a kernel of random genes
      (${bgName(res.background, net)}, N = ${res.N.toLocaleString("en")}) in <b>${Math.round(best.auc * 100)}%</b> of the draws
      (${best.meanList.toFixed(1)} versus ${best.meanBg.toFixed(1)} on average).`
      + (res.k !== res.kAuto ? ` The table uses the kernel size you fixed, k = ${res.k}.` : "");
  };
  $("#summary").innerHTML = VIEW.net === "both" ? `<p>${one("dap", R.dap)}</p><p>${one("any", R.any)}</p>`
    : `<p>${one(VIEW.net, VIEW.net === "dap" ? R.dap : R.any)}</p>`;
  const w = [], r = R.dap, list = (a, s) => `${a.length} ${s}: ${esc(a.slice(0, 12).join(", "))}${a.length > 12 ? "…" : ""}.`;
  if (r.unknown.length) w.push(list(r.unknown, `identifier${r.unknown.length > 1 ? "s" : ""} not found in Araport11`));
  if (VIEW.net !== "any" && r.notInNetwork.length) w.push(list(r.notInNetwork, `gene${r.notInNetwork.length > 1 ? "s have" : " has"} no DAP-seq peak (bound by no TF)`));
  const noT = R.any.genes.filter(g => T.any.index.has(g) && !(T.flags[T.any.index.get(g)] & 4));
  if (VIEW.net !== "dap" && noT.length) w.push(list(noT, `gene${noT.length > 1 ? "s are" : " is"} regulated in no TARGET experiment`));
  if (r.added) w.push(`${r.added} gene${r.added > 1 ? "s" : ""} of your list ${r.added > 1 ? "were" : "was"} not in the background and ${r.added > 1 ? "were" : "was"} added to it.`);
  for (const net of VIEW.net === "both" ? ["dap", "any"] : [VIEW.net]) {
    const res = R[net];
    if (res.profile.find(p => p.k === res.kAuto).auc < 0.57) w.push(`${NETNAME[net]}: the separation is close to 0.5 at every kernel size; its candidates are weak leads.`);
  }
  $("#warnings").innerHTML = w.map(s => `<p class="warn">⚠ ${s}</p>`).join("");
  $("#profbox").style.display = VIEW.net === "both" ? "none" : "";
  $("#kwrap").style.display = VIEW.net === "both" ? "none" : "";
  if (VIEW.net !== "both") {
    const res = cur().res;
    $("#kshow").innerHTML = res.profile.map(p => `<option value="${p.k}"${p.k === VIEW.k[VIEW.net] ? " selected" : ""}>${p.k}${p.k === res.kAuto ? " (auto)" : ""}</option>`).join("");
  }
  $("#rankby").innerHTML = RANKS[VIEW.net].map(([k, l]) => `<button type="button" data-r="${k}" class="${VIEW.by[VIEW.net] === k ? "on" : ""}">${l}</button>`).join("");
  $$("#rankby button").forEach(b => b.addEventListener("click", () => { VIEW.by[VIEW.net] = b.dataset.r; VIEW.sort = null; show(); }));
  $("#candnote").innerHTML = NOTES[VIEW.net];
  drawCharts();
  table();
  genesTable();
  famTable();
}
$$("#netsel button").forEach(b => b.addEventListener("click", () => {
  VIEW.net = b.dataset.n; VIEW.sort = null; VIEW.all = false; VIEW.open = new Set(); show(); updateHashNet();
}));

function setK(k) {
  VIEW.k[VIEW.net] = k;
  $("#kshow").value = k;
  if (VIEW.net === "any") { const o = { ...R.o, profile: false, kernel: k }; R.up = T.up.run(R.genes, o); R.down = T.down.run(R.genes, o); }
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

function drawCharts() {
  if (VIEW.net !== "both") profileChart();
  aucChart();
}

function profileChart() {
  const { res } = cur(), P = res.profile, W = 1060, H = 320, m = { l: 52, r: 12, t: 12, b: 38 }, kSel = VIEW.k[VIEW.net];
  const qs = P.map(p => ({ k: p.k, l: quantiles(p.list), b: quantiles(p.bg), p }));
  const ymax = Math.max(1, ...qs.map(q => Math.max(q.l.hi, q.b.hi, q.l.q3, q.b.q3)));
  const f = v => Math.log1p(v), y = v => H - m.b - (f(v) / f(ymax)) * (H - m.t - m.b);
  const cw = (W - m.l - m.r) / P.length, bw = Math.min(26, cw * 0.3);
  const tg = VIEW.net === "any", ink = css("--sub"), grid = css("--grid"), cl = css(tg ? "--tgt" : "--list"), cb = css("--bgser"),
    fl = css(tg ? "--tgtfill" : "--listfill"), fb = css("--bgfill"), hi = css("--hi");
  $("#profttl").textContent = tg ? "Common TARGET experiments per kernel" : "Common regulators per kernel";
  $("#proflegend").innerHTML = `<span><i style="background:${fl};border:1.5px solid ${cl}"></i>kernels of your list</span><span><i style="background:var(--bgfill);border:1.5px solid var(--bgser)"></i>kernels of random genes</span>`;
  let s = `<svg viewBox="0 0 ${W} ${H}" xmlns="http://www.w3.org/2000/svg" font-family="Source Sans 3, Arial, sans-serif" font-size="12">`;
  const ticks = [0, 1, 2, 5, 10, 20, 50, 100, 200, 400].filter(t => t <= ymax * 1.05);
  for (const t of ticks) s += `<line x1="${m.l}" x2="${W - m.r}" y1="${y(t)}" y2="${y(t)}" stroke="${grid}"/><text x="${m.l - 8}" y="${y(t) + 4}" text-anchor="end" fill="${ink}">${t}</text>`;
  s += `<text transform="translate(14 ${(H - m.b + m.t) / 2}) rotate(-90)" text-anchor="middle" fill="${ink}">${tg ? "common experiments" : "common TFs"} per kernel</text>`;
  const box = (x, q, stroke, fill) => {
    const w2 = bw / 2;
    return `<line x1="${x}" x2="${x}" y1="${y(q.hi)}" y2="${y(q.q3)}" stroke="${stroke}" stroke-width="1.5"/>`
      + `<line x1="${x}" x2="${x}" y1="${y(q.q1)}" y2="${y(q.lo)}" stroke="${stroke}" stroke-width="1.5"/>`
      + `<rect x="${x - w2}" y="${y(q.q3)}" width="${bw}" height="${Math.max(1.5, y(q.q1) - y(q.q3))}" rx="3" fill="${fill}" stroke="${stroke}" stroke-width="1.5"/>`
      + `<line x1="${x - w2}" x2="${x + w2}" y1="${y(q.med)}" y2="${y(q.med)}" stroke="${stroke}" stroke-width="2.5"/>`;
  };
  qs.forEach((q, i) => {
    const x0 = m.l + i * cw, on = q.k === kSel;
    s += `<g class="kcol" data-k="${q.k}"><rect class="kbg" x="${x0 + 2}" y="${m.t}" width="${cw - 4}" height="${H - m.t - m.b + 30}" rx="6" fill="${on ? hi : "transparent"}"/>`;
    s += box(x0 + cw / 2 - bw * 0.65, q.l, cl, fl) + box(x0 + cw / 2 + bw * 0.65, q.b, cb, fb);
    s += `<text x="${x0 + cw / 2}" y="${H - m.b + 20}" text-anchor="middle" fill="${on ? css("--fg") : ink}" font-weight="${on ? 700 : 400}">k = ${q.k}${q.k === res.kAuto ? " ★" : ""}</text></g>`;
  });
  s += "</svg>";
  $("#profile").innerHTML = s;
  $$("#profile .kcol").forEach(g => {
    const q = qs.find(q => q.k === +g.dataset.k);
    g.addEventListener("click", () => setK(q.k));
    g.addEventListener("mousemove", e => showTip(e, `<b>kernel size ${q.k}</b>${q.k === res.kAuto ? " (automatic)" : ""}<br>`
      + `your list: median ${q.l.med}, mean ${q.p.meanList.toFixed(2)} common ${tg ? "experiments" : "TFs"}<br>random genes: median ${q.b.med}, mean ${q.p.meanBg.toFixed(2)}<br>`
      + `AUC ${q.p.auc.toFixed(3)} · ratio of means ×${(q.p.meanList / Math.max(q.p.meanBg, 1e-9)).toFixed(1)}<br><span class="note">click to use this kernel size</span>`));
    g.addEventListener("mouseleave", hideTip);
  });
}

function aucChart() {
  const nets = VIEW.net === "both" ? ["dap", "any"] : [VIEW.net];
  const series = nets.map(n => ({ n, P: R[n].profile, res: R[n], color: css(n === "dap" ? "--list" : "--tgt") }));
  const P0 = series[0].P, W = 1060, H = 160, m = { l: 52, r: 12, t: 16, b: 30 }, all = series.flatMap(s => s.P.map(p => p.auc));
  const lo = Math.min(0.5, ...all) - 0.02, hiV = Math.max(...all) + 0.04;
  const cw = (W - m.l - m.r) / P0.length, x = i => m.l + i * cw + cw / 2, y = v => H - m.b - ((v - lo) / (hiV - lo)) * (H - m.t - m.b);
  const ink = css("--sub"), grid = css("--grid"), bg = css("--card");
  let s = `<svg viewBox="0 0 ${W} ${H}" xmlns="http://www.w3.org/2000/svg" font-family="Source Sans 3, Arial, sans-serif" font-size="12">`;
  const step = hiV - lo > 0.3 ? 0.1 : 0.05;
  for (let t = Math.ceil(lo / step) * step; t <= hiV; t += step) s += `<line x1="${m.l}" x2="${W - m.r}" y1="${y(t)}" y2="${y(t)}" stroke="${grid}"/><text x="${m.l - 8}" y="${y(t) + 4}" text-anchor="end" fill="${ink}">${t.toFixed(2)}</text>`;
  s += `<line x1="${m.l}" x2="${W - m.r}" y1="${y(0.5)}" y2="${y(0.5)}" stroke="${ink}" stroke-dasharray="4 4"/><text x="${m.l + 4}" y="${y(0.5) - 5}" fill="${ink}">random (0.5)</text>`;
  for (const se of series) s += `<polyline fill="none" stroke="${se.color}" stroke-width="2" points="${se.P.map((p, i) => `${x(i)},${y(p.auc)}`).join(" ")}"/>`;
  P0.forEach((p0, i) => {
    s += `<g class="kcol" data-k="${p0.k}"><rect x="${m.l + i * cw}" y="0" width="${cw}" height="${H}" fill="transparent"/>`;
    for (const se of series) {
      const p = se.P[i]; if (!p) continue;
      const on = VIEW.net !== "both" && p.k === VIEW.k[se.n], best = p.k === se.res.kAuto;
      s += `<circle cx="${x(i)}" cy="${y(p.auc)}" r="${on ? 7 : 5}" fill="${best || on ? se.color : bg}" stroke="${on ? bg : se.color}" stroke-width="2"/>`
        + (best ? `<text x="${x(i)}" y="${se.n === "any" && series.length > 1 ? y(p.auc) + 20 : y(p.auc) - 11}" text-anchor="middle" fill="${se.color}" font-weight="700">${p.auc.toFixed(3)} ★</text>` : "");
    }
    s += `<text x="${x(i)}" y="${H - 6}" text-anchor="middle" fill="${ink}">k = ${p0.k}</text></g>`;
  });
  s += "</svg>";
  $("#aucplot").innerHTML = s;
  $("#auclegend").innerHTML = series.map(se => `<span><i style="background:${se.color}"></i>${NETNAME[se.n]}</span>`).join("");
  $$("#aucplot .kcol").forEach(g => {
    const k = +g.dataset.k;
    if (VIEW.net !== "both") g.addEventListener("click", () => setK(k));
    g.addEventListener("mousemove", e => showTip(e, `<b>kernel size ${k}</b><br>` + series.map(se => {
      const p = se.P.find(p => p.k === k);
      return p ? `${NETNAME[se.n]}: AUC ${p.auc.toFixed(3)}${k === se.res.kAuto ? " (highest: automatic choice)" : ""}` : "";
    }).join("<br>")));
    g.addEventListener("mouseleave", hideTip);
  });
}

// --------------------------------------------------------------------------- candidate tables
const RANKS = {
  dap: [["p", "FDR"], ["score", "Régine score"]],
  any: [["p", "FDR"], ["pu", "FDR induced"], ["pd", "FDR repressed"], ["score", "Régine score"]],
  both: [["c", "combined"], ["dp", "DAP-seq"], ["tp", "TARGET"], ["ip", "bound ∩ regulated"]],
};
const NOTES = {
  dap: "One row per transcription factor of the DAP-seq cistrome (O'Malley et al. 2016): genes with a binding peak in vitro.",
  any: "One row per TARGET experiment (ConnecTF): genes induced (↑) or repressed (↓) by the TF in root or shoot protoplasts, "
    + "with cycloheximide (direct targets). A TF tested several times has several rows. FDRs are corrected over the experiments, for each of the three tests.",
  both: "One row per transcription factor, DAP-seq and TARGET side by side. TARGET: the best experiment of the TF (its P is multiplied by the number of "
    + "experiments of the TF). <i>Bound ∩ regulated</i>: genes both bound in DAP-seq and regulated in TARGET (TFs present in both). "
    + "Combined: Fisher's method on the DAP-seq and TARGET P values when the TF is in both, otherwise the P of its only network; FDR over all TFs.",
};
const tfCell = (t, q, badge = "FDR") => `<td class="tf"><b>${esc(t.sym || "—")}</b>${q <= 0.05 ? `<span class="sig" title="FDR ≤ 5%">${badge}</span>` : ""}<br>`
  + `<a href="https://www.arabidopsis.org/locus?name=${t.id}" target="_blank" rel="noopener" onclick="event.stopPropagation()">${t.id}</a></td>`;
const bar = (x, of, cls = "") => `<span class="bar ${cls}" style="width:${Math.round(40 * x / Math.max(of, 1))}px"></span>`;

// rows of the current view, ranked; every row has key, tf (annotation), fam, q (for filters) and the sort keys
function rankedRows(net = VIEW.net) {
  let rows;
  if (net === "dap") {
    rows = DAP.score(R.dap, VIEW.k.dap).map(r => ({ ...r, key: r.tf.id, fam: r.tf.fam }));
    rows.sort(VIEW.by.dap === "p" ? (a, b) => a.p - b.p || b.score - a.score : (a, b) => b.score - a.score || a.p - b.p);
  } else if (net === "any") {
    rows = T.any.score(R.any, VIEW.k.any).map(r => ({ ...r, key: r.tf.key, fam: r.tf.fam,
      ku: R.up.cl[r.t], Ku: R.up.cb[r.t], pu: R.up.p[r.t], qu: R.up.q[r.t],
      kd: R.down.cl[r.t], Kd: R.down.cb[r.t], pd: R.down.p[r.t], qd: R.down.q[r.t] }));
    rows.forEach(r => (r.eff = effect(r.qu, r.qd)));
    const by = VIEW.by.any;
    rows.sort(by === "score" ? (a, b) => b.score - a.score || a.p - b.p : (a, b) => a[by] - b[by] || a.p - b.p);
  } else rows = bothRows();
  rows.forEach((r, i) => (r.rank = i + 1));
  return rows;
}

function bothRows() {
  const byTF = new Map(), get = (id, tf) => byTF.get(id) || byTF.set(id, { key: id, tf, fam: tf.fam, exps: [], inter: [] }).get(id);
  DAP.score(R.dap, R.dap.k).forEach(r => Object.assign(get(r.tf.id, r.tf), { d: r }));
  const t = rankedRows("any"), ti = R.inter;
  for (const r of t) get(r.tf.id, r.tf).exps.push(r);
  T.inter.tfs.forEach((x, j) => get(x.id, x).inter.push({ x, kt: ti.cl[j], Kt: ti.cb[j], p: ti.p[j], q: ti.q[j] }));
  const rows = [...byTF.values()];
  for (const r of rows) {
    if (r.d) { r.dp = r.d.p; r.dq = r.d.q; }
    if (r.exps.length) {
      r.best = r.exps.reduce((a, b) => (b.p < a.p ? b : a));
      r.tp = Math.min(1, r.best.p * r.exps.length); r.tq = r.best.q; r.eff = r.best.eff;
    }
    if (r.inter.length) { r.bi = r.inter.reduce((a, b) => (b.p < a.p ? b : a)); r.ip = Math.min(1, r.bi.p * r.inter.length); r.iq = r.bi.q; }
    r.c = r.d && r.best ? fisher2(r.dp, r.tp) : r.d ? r.dp : r.tp;
  }
  const q = bh(rows.map(r => r.c));
  rows.forEach((r, i) => (r.q = r.cq = q[i]));
  const by = VIEW.by.both, v = r => r[by] ?? 2;
  rows.sort((a, b) => v(a) - v(b) || a.c - b.c);
  return rows;
}

const COLS = {
  dap: [
    ["sym", "transcription factor"], ["fam", "family", "hide-s"], [null, "description", "hide-s"],
    ["kt", "your genes", "num", "genes of your list bound by the TF"], ["Kt", "background", "num", "fraction of the background bound by the TF"],
    ["real", "real", "num", "fraction of kernels of your list entirely bound by the TF"], ["noise", "noise", "num", "fraction of random kernels entirely bound by the TF"],
    ["score", "score", "num", "real − noise"], ["fold", "fold", "num", "enrichment of the TF's targets in your list"], ["q", "FDR", "num", "hypergeometric test, Benjamini–Hochberg"],
  ],
  any: [
    ["sym", "transcription factor"], [null, "experiment", "hide-s"], ["fam", "family", "hide-s"],
    ["kt", "your genes", "num", "genes of your list regulated in this experiment"], ["ku", "↑ induced", "num", "genes of your list induced by the TF (FDR)"],
    ["kd", "↓ repressed", "num", "genes of your list repressed by the TF (FDR)"], ["Kt", "background", "num", "fraction of the background regulated in this experiment"],
    ["score", "score", "num", "real − noise, on all regulated genes"], ["fold", "fold", "num", "enrichment of the regulated genes in your list"],
    ["q", "FDR", "num", "hypergeometric test on all regulated genes, Benjamini–Hochberg"], ["eff", "effect", "", "FDR ≤ 5% for the induced or the repressed genes"],
  ],
  both: [
    ["sym", "transcription factor"], ["fam", "family", "hide-s"],
    ["dp", "DAP-seq", "num", "genes of your list bound in DAP-seq; FDR"], ["tp", "TARGET", "num", "best TARGET experiment of the TF: genes of your list regulated; FDR"],
    ["eff", "effect", "", "best TARGET experiment: FDR ≤ 5% for the induced or the repressed genes"],
    ["ip", "bound ∩ regulated", "num", "genes of your list bound in DAP-seq and regulated in TARGET (best experiment); FDR"],
    ["c", "combined FDR", "num", "Fisher's method on the DAP-seq and TARGET P values; Benjamini–Hochberg over all TFs"],
  ],
};

$("#famsel").addEventListener("change", () => table());
$("#onlysig").addEventListener("change", () => table());
$("#q").addEventListener("input", () => table());
$("#showall").addEventListener("click", () => { VIEW.all = !VIEW.all; table(); });

function rowCells(r, net) {
  const n = R.dap.n;
  if (net === "dap") return tfCell(r.tf, r.q) + `<td class="hide-s"><span class="fam">${esc(r.tf.fam || "?")}</span></td><td class="desc hide-s">${tfName(r.tf)}</td>`
    + `<td class="num">${r.kt}/${R.dap.n}${bar(r.kt, R.dap.n)}</td><td class="num">${pct(r.Kt / R.dap.N, 0)}${bar(r.Kt, R.dap.N, "barbg")}</td>`
    + `<td class="num">${pct(r.real)}</td><td class="num">${pct(r.noise)}</td><td class="num"><b>${(100 * r.score).toFixed(1)}</b></td>`
    + `<td class="num">×${r.fold.toFixed(1)}</td><td class="num">${fmtP(r.q)}</td>`;
  if (net === "any") {
    const x = r.tf, a = R.any;
    return tfCell(x, Math.min(r.q, r.qu, r.qd)) + `<td class="hide-s exp">${esc(x.source)} · ${esc(x.tissue)}<br><span class="note">${x.n.toLocaleString("en")} targets (${x.up.toLocaleString("en")} ↑, ${x.down.toLocaleString("en")} ↓)</span></td>`
      + `<td class="hide-s"><span class="fam">${esc(x.fam || "?")}</span></td>`
      + `<td class="num">${r.kt}/${a.n}${bar(r.kt, a.n, "bartgt")}</td>`
      + `<td class="num up">${r.ku}<br><span class="note">${fmtP(r.qu)}</span></td><td class="num dn">${r.kd}<br><span class="note">${fmtP(r.qd)}</span></td>`
      + `<td class="num">${pct(r.Kt / a.N, 0)}${bar(r.Kt, a.N, "barbg")}</td><td class="num"><b>${(100 * r.score).toFixed(1)}</b></td>`
      + `<td class="num">×${r.fold.toFixed(1)}</td><td class="num">${fmtP(r.q)}</td><td>${effectHtml(r.eff)}</td>`;
  }
  const d = r.d, b = r.best, i = r.bi;
  return tfCell(r.tf, r.cq, d && b && r.dq <= 0.05 && r.tq <= 0.05 ? "BOTH" : "FDR") + `<td class="hide-s"><span class="fam">${esc(r.fam || "?")}</span></td>`
    + (d ? `<td class="num">${d.kt}/${n}${bar(d.kt, n)}<br><span class="note">${fmtP(r.dq)}</span></td>` : `<td class="num note">not in DAP-seq</td>`)
    + (b ? `<td class="num">${b.kt}/${R.any.n}${bar(b.kt, R.any.n, "bartgt")}<br><span class="note">${fmtP(r.tq)} · ${esc(b.tf.source)}, ${esc(b.tf.tissue)}${r.exps.length > 1 ? ` (best of ${r.exps.length})` : ""}</span></td>` : `<td class="num note">not in TARGET</td>`)
    + `<td>${effectHtml(r.eff)}</td>`
    + (i ? `<td class="num">${i.kt}/${n}<br><span class="note">${fmtP(r.iq)}</span></td>` : `<td class="num note">—</td>`)
    + `<td class="num"><b>${fmtP(r.cq)}</b></td>`;
}

function table() {
  const net = VIEW.net, cols = COLS[net];
  $("#tftab thead").innerHTML = "<tr><th>#</th>" + cols.map(([k, l, cls = "", title = ""]) =>
    `<th${k ? ` data-s="${k}"` : ""} class="${cls}${k && k === VIEW.sort ? " sorted" : ""}"${title ? ` title="${esc(title)}"` : ""}>${l}</th>`).join("") + "</tr>";
  $$("#tftab thead th[data-s]").forEach(th => th.addEventListener("click", () => {
    const key = th.dataset.s;
    VIEW.dir = VIEW.sort === key ? -VIEW.dir : (["sym", "fam", "q", "dp", "tp", "ip", "c", "eff"].includes(key) ? 1 : -1);
    VIEW.sort = key; table();
  }));
  const rows = rankedRows();
  const fam = $("#famsel").value, sig = $("#onlysig").checked, q = $("#q").value.trim().toLowerCase();
  const hit = r => net === "dap" ? r.kt > 0 : net === "any" ? r.kt > 0 : (r.d?.kt > 0 || r.best?.kt > 0);
  const sigOf = r => net === "any" ? Math.min(r.q, r.qu, r.qd) : r.q;
  let shown = rows.filter(r => hit(r) && (!fam || r.fam === fam) && (!sig || sigOf(r) <= 0.05)
    && (!q || [r.tf.id, r.tf.sym, r.tf.fam, r.tf.name, r.tf.desc, r.tf.summary, r.tf.source].join(" ").toLowerCase().includes(q)));
  if (VIEW.sort) {
    const k = VIEW.sort, get = r => k === "sym" ? (r.tf.sym || r.tf.id) : k === "fam" ? (r.fam || "") : r[k] ?? (VIEW.dir > 0 ? Infinity : -Infinity);
    shown.sort((a, b) => { const x = get(a), y = get(b); return (typeof x === "string" ? x.localeCompare(y) : x - y) * VIEW.dir; });
  }
  const total = shown.length, limit = VIEW.all ? total : 40, ncol = cols.length + 1;
  shown = shown.slice(0, limit);
  $("#tftab tbody").innerHTML = shown.map(r => {
    const open = VIEW.open.has(r.key);
    return `<tr class="tfrow${open ? " open" : ""}" data-t="${esc(r.key)}"><td class="rank">${r.rank}</td>${rowCells(r, net)}</tr>` + (open ? detail(r, ncol) : "");
  }).join("") || `<tr><td colspan="${ncol}" class="note">No transcription factor matches.</td></tr>`;
  const unit = net === "any" ? "experiments regulating" : "TFs binding or regulating";
  $("#showall").textContent = total > 40 ? (VIEW.all ? "show the first 40 only" : `show all ${total} ${unit} at least one of your genes`) : "";
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

const listGenes = () => R.genes.filter(g => DAP.index.has(g));
const chips = (list, cls) => {
  const order = [...list.filter(g => cls(g) !== "off"), ...list.filter(g => cls(g) === "off")].slice(0, 300);
  return `<div class="genes">${order.map(g => `<span class="${cls(g)}" title="${esc(gdesc(g))}">${esc(sym(g) || g)}${sym(g) ? ` <span class="note">${g}</span>` : ""}</span>`).join("")}${list.length > 300 ? " …" : ""}</div>`;
};
const inSet = (reg, res, t) => new Set(reg.targetsOf(res, t));
const links = id => `<a href="https://www.arabidopsis.org/locus?name=${id}" target="_blank" rel="noopener">TAIR</a> ·
  <a href="https://www.uniprot.org/uniprotkb?query=${id}" target="_blank" rel="noopener">UniProt</a>`;
const annot = t => `<h4>${esc(t.sym || t.id)} · ${esc(t.fam || "family unknown")} · ${t.id}</h4>
  <p>${esc(t.summary || t.desc || "No curator summary.")}</p>${t.uniprot ? `<h4 style="margin-top:10px">UniProt</h4><p>${esc(t.uniprot)}</p>` : ""}`;
const expText = x => `${esc(x.source)}, ${esc(x.tissue)} protoplasts, ${esc(x.treatment)}, ${esc(x.time)}, ${esc(x.cutoff)}`
  + (x.doi ? ` · <a href="https://doi.org/${esc(x.doi)}" target="_blank" rel="noopener">doi</a>` : "");

function detail(r, ncol) {
  const list = listGenes();
  let left, right;
  if (VIEW.net === "dap") {
    const bound = inSet(DAP, R.dap, r.t), t = r.tf;
    left = annot(t) + `<p class="note" style="margin-top:10px">Binds ${t.deg.toLocaleString("en")} genes genome-wide (${pct(t.deg / META.n_network_genes, 0)} of the network) ·
      ${r.kt} of your ${R.dap.n} genes · ${r.Kt.toLocaleString("en")} of ${R.dap.N.toLocaleString("en")} background genes ·
      enrichment ×${r.fold.toFixed(2)}, P = ${fmtP(r.p)}, FDR = ${fmtP(r.q)} · ${links(t.id)}</p>`;
    right = `<h4>Your genes bound by ${esc(t.sym || t.id)} (${bound.size}/${R.dap.n}); crossed out: not bound</h4>` + chips(list, g => bound.has(g) ? "" : "off");
  } else if (VIEW.net === "any") {
    const x = r.tf, up = inSet(T.up, R.up, r.t), dn = inSet(T.down, R.down, r.t);
    left = annot(x) + `<p class="note" style="margin-top:10px">TARGET: ${expText(x)} · ${x.n.toLocaleString("en")} regulated genes (${x.up.toLocaleString("en")} induced, ${x.down.toLocaleString("en")} repressed) ·
      your genes: ${r.kt} regulated (P = ${fmtP(r.p)}, FDR = ${fmtP(r.q)}), ${r.ku} induced (FDR ${fmtP(r.qu)}), ${r.kd} repressed (FDR ${fmtP(r.qd)}) · ${links(x.id)}</p>`;
    right = `<h4>Your genes regulated by ${esc(x.sym || x.id)}: <span class="up">↑ induced</span>, <span class="dn">↓ repressed</span>, crossed out: not regulated</h4>`
      + chips(list, g => up.has(g) ? "up" : dn.has(g) ? "dn" : "off");
  } else {
    const t = r.tf, bound = r.d ? inSet(DAP, R.dap, r.d.t) : new Set();
    const up = r.best ? inSet(T.up, R.up, r.best.t) : new Set(), dn = r.best ? inSet(T.down, R.down, r.best.t) : new Set();
    left = annot(t) + `<p class="note" style="margin-top:10px">`
      + (r.d ? `DAP-seq: ${r.d.kt}/${R.dap.n} of your genes bound, FDR ${fmtP(r.dq)}. ` : "Not in the DAP-seq set. ")
      + (r.best ? `TARGET (${r.exps.length} experiment${r.exps.length > 1 ? "s" : ""}): best ${expText(r.best.tf)}, ${r.best.kt}/${R.any.n} regulated, FDR ${fmtP(r.tq)}. ` : "Not in TARGET. ")
      + (r.bi ? `Bound and regulated: ${r.bi.kt}, FDR ${fmtP(r.iq)}. ` : "") + `Combined FDR ${fmtP(r.cq)} · ${links(t.id)}</p>`;
    right = `<h4>Your genes: <span class="bound">bound in DAP-seq</span>; TARGET <span class="up">↑ induced</span> / <span class="dn">↓ repressed</span>; crossed out: neither</h4>`
      + chips(list, g => { const c = [bound.has(g) ? "bound" : "", up.has(g) ? "up" : dn.has(g) ? "dn" : ""].filter(Boolean).join(" "); return c || "off"; });
  }
  return `<tr class="detail"><td colspan="${ncol}"><div class="det"><div>${left}</div><div>${right}</div></div></td></tr>`;
}

// --------------------------------------------------------------------------- gene and family tables
function genesTable() {
  const net = VIEW.net, list = listGenes(), top = rankedRows().slice(0, 20);
  const dapBits = (row, t) => DAP.mask[row * DAP.W + (t >> 5)] & (1 << (t & 31));
  const tBits = (reg, row, e) => reg.mask[row * reg.W + (e >> 5)] & (1 << (e & 31));
  const mark = (row, r) => {
    if (net === "dap") return dapBits(row, r.t) ? esc(r.tf.sym || r.tf.id) : "";
    if (net === "any") return tBits(T.up, row, r.t) ? `${esc(r.tf.sym)}<span class="up">↑</span>` : tBits(T.down, row, r.t) ? `${esc(r.tf.sym)}<span class="dn">↓</span>` : "";
    const b = r.d && dapBits(row, r.d.t), e = r.best?.t, u = e != null && tBits(T.up, row, e), d = e != null && tBits(T.down, row, e);
    return b || u || d ? `${esc(r.tf.sym || r.tf.id)}${b ? '<span class="bound">•</span>' : ""}${u ? '<span class="up">↑</span>' : d ? '<span class="dn">↓</span>' : ""}` : "";
  };
  $("#gtab tbody").innerHTML = list.slice(0, 3000).map(g => {
    const row = DAP.index.get(g), by = top.map(r => mark(row, r)).filter(Boolean);
    return `<tr><td><a href="https://www.arabidopsis.org/locus?name=${g}" target="_blank" rel="noopener">${g}</a></td><td><b>${esc(sym(g))}</b></td>`
      + `<td class="desc">${esc(gdesc(g))}</td><td class="num">${DAP.nTF(row)}</td><td class="num">${T.any.nTF(row)}</td>`
      + `<td>${by.slice(0, 8).join(", ")}${by.length > 8 ? ` +${by.length - 8}` : ""}</td></tr>`;
  }).join("");
}

function famTable() {
  const rows = rankedRows(), by = new Map();
  rows.forEach((r, i) => {
    const f = r.fam || "unknown";
    const e = by.get(f) || by.set(f, { n: 0, sig: 0, top: 0, best: [] }).get(f);
    const q = VIEW.net === "any" ? Math.min(r.q, r.qu, r.qd) : r.q;
    e.n++; if (q <= 0.05) e.sig++; if (i < 20) e.top++;
    const s = r.tf.sym || r.tf.id;
    if (e.best.length < 4 && !e.best.includes(s) && (r.kt > 0 || r.d?.kt > 0 || r.best?.kt > 0)) e.best.push(s);
  });
  $("#ftab thead").innerHTML = `<tr><th>family</th><th class="num">${VIEW.net === "any" ? "experiments" : "TFs"}</th><th class="num">FDR ≤ 5%</th><th class="num">in the top 20</th><th>best TFs</th></tr>`;
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
  const net = VIEW.net, rows = rankedRows(), ids = (reg, res, t) => reg.targetsOf(res, t).join(",");
  let head, body, hdr;
  if (net === "dap") {
    head = ["rank", "agi", "symbol", "family", "your_genes_bound", "n_list", "background_bound", "N_background", "real", "noise", "regine_score",
      "fold", "p", "fdr", "tf_targets_genome", "protein_name", "short_description", "summary", "your_genes_bound_ids"];
    body = rows.map(r => [r.rank, r.tf.id, r.tf.sym, r.tf.fam, r.kt, R.dap.n, r.Kt, R.dap.N, r.real, r.noise, r.score, r.fold, r.p, r.q, r.tf.deg,
      r.tf.name, r.tf.desc, r.tf.summary, ids(DAP, R.dap, r.t)]);
    hdr = `DAP-seq · kernel size ${VIEW.k.dap} · background ${R.dap.background} (N=${R.dap.N})`;
  } else if (net === "any") {
    head = ["rank", "experiment", "agi", "symbol", "family", "source", "tissue", "your_genes_regulated", "n_list", "background_regulated", "N_background",
      "real", "noise", "regine_score", "fold", "p", "fdr", "induced", "fdr_induced", "repressed", "fdr_repressed", "effect", "experiment_targets",
      "your_genes_induced_ids", "your_genes_repressed_ids"];
    body = rows.map(r => [r.rank, r.tf.key, r.tf.id, r.tf.sym, r.tf.fam, r.tf.source, r.tf.tissue, r.kt, R.any.n, r.Kt, R.any.N, r.real, r.noise, r.score,
      r.fold, r.p, r.q, r.ku, r.qu, r.kd, r.qd, r.eff, r.tf.n, ids(T.up, R.up, r.t), ids(T.down, R.down, r.t)]);
    hdr = `TARGET · kernel size ${VIEW.k.any} · background ${R.any.background} (N=${R.any.N})`;
  } else {
    head = ["rank", "agi", "symbol", "family", "dap_your_genes_bound", "dap_p", "dap_fdr", "target_best_experiment", "target_n_experiments",
      "target_your_genes_regulated", "target_p_adjusted", "target_fdr", "target_effect", "bound_and_regulated", "bound_and_regulated_fdr", "combined_p", "combined_fdr"];
    body = rows.map(r => [r.rank, r.tf.id, r.tf.sym, r.fam, r.d?.kt ?? "", r.dp ?? "", r.dq ?? "", r.best?.tf.key ?? "", r.exps.length, r.best?.kt ?? "",
      r.tp ?? "", r.tq ?? "", r.eff ?? "", r.bi?.kt ?? "", r.iq ?? "", r.c, r.cq]);
    hdr = `DAP-seq + TARGET side by side · background ${R.dap.background}`;
  }
  download(`regine_${net === "any" ? "target" : net}.tsv`, `# Régine 2026 · ${hdr}\n` + head.join("\t") + "\n" + body.map(r => r.join("\t")).join("\n") + "\n", "text/tab-separated-values");
});
$("#dlprof").addEventListener("click", () => {
  const head = "network\tkernel_size\tauc\texact_mean_list\texact_mean_background\tlist_median\tbackground_median\tlist_q1\tlist_q3\tbackground_q1\tbackground_q3\tdraws";
  const body = ["dap", "any"].flatMap(n => R[n].profile.map(p => { const l = quantiles(p.list), b = quantiles(p.bg); return [NETNAME[n], p.k, p.auc, p.meanList, p.meanBg, l.med, b.med, l.q1, l.q3, b.q1, b.q3, R[n].reps].join("\t"); }));
  download("regine_kernel_profile.tsv", head + "\n" + body.join("\n") + "\n", "text/tab-separated-values");
});
$("#dlsvg").addEventListener("click", () => {
  const svg = $(VIEW.net === "both" ? "#aucplot svg" : "#profile svg").outerHTML.replace("<svg ", `<svg style="background:${css("--card")}" `);
  download(VIEW.net === "both" ? "regine_auc.svg" : "regine_kernel_profile.svg", svg, "image/svg+xml");
});
$("#share").addEventListener("click", async () => {
  try { await navigator.clipboard.writeText(location.href); $("#share").textContent = "Link copied"; }
  catch { $("#share").textContent = "Copy the address bar"; }
  setTimeout(() => ($("#share").textContent = "Copy link to this analysis"), 2000);
});

function writeHash(genes, o) {
  if (o.background === "custom") { history.replaceState(null, "", location.pathname); return; }   // too long for a URL
  const p = new URLSearchParams({ g: genes.join(","), bg: o.background, k: o.kernel, r: o.reps, s: o.seed, n: VIEW.net });
  history.replaceState(null, "", "#" + p.toString());
}
function updateHashNet() {
  if (!location.hash.includes("g=")) return;
  const p = new URLSearchParams(location.hash.slice(1)); p.set("n", VIEW.net);
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
  if (NETNAME[p.get("n")]) VIEW.net = p.get("n");
  countGenes();
  $("#form").requestSubmit();
}
