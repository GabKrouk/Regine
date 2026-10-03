// Régine 2026 — bridge with GeneCloud (both are static pages of https://gabkrouk.github.io, so no server is needed).
//
//   GeneCloud -> Régine   click a GO term: GeneCloud opens Régine with the genes annotated with this term.
//   Régine -> GeneCloud   send the genes bound by a TF (or the best TFs) to GeneCloud, to see what they have in common.
//
// Two transports, same payload {genes, name, go, from}:
//   * URL hash  #from=genecloud&n=<label>&go=<GO:id>&g=<AGI,AGI,…>   (short lists, works with a plain link)
//   * postMessage handshake  (any list size): the opened page posts "<app>:ready" to its opener, the opener answers
//     "<other app>:genes" with the payload. Only messages from ORIGIN are accepted.
// The protocol is described in INTEGRATION.md.
import { parseGenes } from "./engine.js";

export const ORIGIN = "https://gabkrouk.github.io";
export const GENECLOUD = ORIGIN + "/genecloud/";
export const MAX_HASH_GENES = 400;          // above this, a URL gets too long: use the handshake only

export function parseIncoming(hash) {
  const p = new URLSearchParams(String(hash).replace(/^#/, ""));
  const genes = parseGenes(p.get("g") || "");
  const info = { genes, name: p.get("n") || "", go: p.get("go") || "", from: p.get("from") || "", handshake: p.get("hs") === "1" };
  return genes.length || info.handshake ? info : null;
}

export function buildHash({ genes = [], name = "", go = "", from = "regine", handshake = false }) {
  const p = new URLSearchParams({ from });
  if (name) p.set("n", name);
  if (go) p.set("go", go);
  if (handshake) p.set("hs", "1");
  if (genes.length && genes.length <= MAX_HASH_GENES) p.set("g", genes.join(","));
  return "#" + p.toString();
}

export function validPayload(d) {
  return !!d && typeof d === "object" && d.type === "genecloud:genes" && Array.isArray(d.genes) && parseGenes(d.genes.join(" ")).length > 1;
}

// Regine is the opened page: tell the opener we can receive, then hand the payload to onGenes.
export function listenForGenes(onGenes) {
  addEventListener("message", e => {
    if (e.origin !== ORIGIN || !validPayload(e.data)) return;
    onGenes({ genes: parseGenes(e.data.genes.join(" ")), name: String(e.data.name || ""), go: String(e.data.go || ""), from: "genecloud" });
  });
}
export function announceReady() {
  if (window.opener) window.opener.postMessage({ type: "regine:ready" }, ORIGIN);
}

// Regine is the opener: open GeneCloud and send it the genes (handshake, falling back to the URL hash).
export function sendToGeneCloud(payload) {
  const url = GENECLOUD + buildHash({ ...payload, handshake: true });
  const w = window.open(url, "genecloud");
  if (!w) return false;
  const onMsg = e => {
    if (e.source !== w || e.origin !== ORIGIN || !e.data || e.data.type !== "genecloud:ready") return;
    removeEventListener("message", onMsg);
    w.postMessage({ type: "regine:genes", genes: payload.genes, name: payload.name || "", from: "regine" }, ORIGIN);
  };
  addEventListener("message", onMsg);
  setTimeout(() => removeEventListener("message", onMsg), 20000);
  return true;
}
