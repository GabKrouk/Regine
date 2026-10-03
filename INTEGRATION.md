# Régine ⇄ GeneCloud

Both programs are static pages of `https://gabkrouk.github.io`, so they talk to each other from the browser, with no server.

```
GeneCloud  ──click a GO term──▶  Régine   (genes of the term → ranked DAP-seq transcription factors)
GeneCloud  ◀──send genes───────  Régine   (genes of your list bound by a TF → their GO terms, text, annotation)
```

## Payload

`{ genes: ["AT1G13300", …], name: "response to nitrate", go: "GO:0010167", from: "genecloud" }`

Two transports, the same on both sides (code: [`web/bridge.js`](web/bridge.js)):

| | |
|---|---|
| URL hash | `#from=genecloud&go=GO:0010167&n=response%20to%20nitrate&g=AT1G13300,AT5G67420,…` — plain link, up to 400 genes |
| `postMessage` handshake | any list size. The opener opens the page with `#from=…&hs=1`; the opened page posts `{type:"regine:ready"}` (or `{type:"genecloud:ready"}`) to `window.opener`; the opener answers `{type:"genecloud:genes", genes, name, go}` (or `{type:"regine:genes", …}`). Messages are accepted only from `https://gabkrouk.github.io`. |

Régine side: done (`web/bridge.js`, `web/app.js`). A banner shows where the list comes from, the analysis starts by itself,
and every TF has a *GeneCloud* link plus a *Send to GeneCloud* button for the 20 best TFs.

## What GeneCloud has to add (about 30 lines)

**1. A "Régine" button on each GO term**, sending the genes annotated with the term:

```js
const REGINE = "https://gabkrouk.github.io/Regine/";
function sendToRegine(genes, name, go) {
  const w = window.open(REGINE + "#from=genecloud&hs=1&go=" + encodeURIComponent(go) + "&n=" + encodeURIComponent(name), "regine");
  const on = e => {
    if (e.source !== w || e.origin !== "https://gabkrouk.github.io" || e.data?.type !== "regine:ready") return;
    removeEventListener("message", on);
    w.postMessage({ type: "genecloud:genes", genes, name, go }, "https://gabkrouk.github.io");
  };
  addEventListener("message", on);
}
// e.g. goTermEl.onclick = () => sendToRegine(genesOfTerm(go), goLabel, go);
```

**2. Receive lists from Régine** (the inverse):

```js
addEventListener("message", e => {
  if (e.origin !== "https://gabkrouk.github.io" || e.data?.type !== "regine:genes") return;
  loadGeneList(e.data.genes, e.data.name);          // GeneCloud's own function that takes a gene list
});
if (opener && location.hash.includes("hs=1")) opener.postMessage({ type: "genecloud:ready" }, "https://gabkrouk.github.io");
// short lists also arrive as  #from=regine&n=…&g=AT1G…,AT5G…
```

Notes: use the same AGI identifiers on both sides (Araport11); Régine keeps only the genes of the DAP-seq network
and tells you which ones it dropped. For a GO term the background should ideally be the genes of the GeneCloud annotation
(*My own background* in Régine); a link parameter for that can be added later.
