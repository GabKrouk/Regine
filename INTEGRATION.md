# Régine ⇄ GeneCloud

Both programs are static pages of `https://gabkrouk.github.io`, so they talk to each other from the browser, with no server.

```
GeneCloud  ──search / term────▶  Régine   (genes found by an annotation search, or carrying a term → ranked TFs)
GeneCloud  ◀──send genes───────  Régine   (genes of your list bound by a TF → their GO terms, text, annotation)
```

## Payload

`{ genes: ["AT1G13300", …], name: "response to nitrate", go: "GO:0010167", from: "genecloud" }`

Two transports, the same on both sides (code: [`web/bridge.js`](web/bridge.js)):

| | |
|---|---|
| URL hash | `#from=genecloud&go=GO:0010167&name=response%20to%20nitrate&g=AT1G13300,AT5G67420,…` — plain link, up to 400 genes |
| `postMessage` handshake | any list size. The opener opens the page with `#from=…&hs=1`; the opened page posts `{type:"regine:ready"}` (or `{type:"genecloud:ready"}`) to `window.opener`; the opener answers `{type:"genecloud:genes", genes, name, go}` (or `{type:"regine:genes", …}`). Messages are accepted only from `https://gabkrouk.github.io`. |

The list name travels as `name=` (`n=` is the network view of Régine's own links).

Régine side: done (`web/bridge.js`, `web/app.js`). A banner shows where the list comes from, the analysis starts by itself
and the page scrolls to the TFs; every TF (DAP-seq, TARGET or side by side) has a *→ GeneCloud* link for its genes of the
list, and *Send to GeneCloud* sends the genes of the significant TFs (20 at most).

## GeneCloud side

Done in `GabKrouk/genecloud`, branch `claude/regine-bridge` (the web app lives on `gh-pages`; see `regine.js` there):
*Search the annotation* finds every gene whose annotation (symbol, descriptions, curator summary, UniProt, phenotypes,
GO term names) matches a boolean query — `"nitrate" AND "xylem" OR "auxin"`, `(nitrate OR nitrogen) AND transport*`,
`"root hair" NOT ethylene`; ET / OU / SAUF work too — and *Send to Régine* sends the genes found (any number, `search.js`);
every term of the table gets a *→ Régine* link, the gene panel of a term has a link too (GO terms are sent with their `GO:` identifier),
and *Send the list to Régine* sends the whole list. A list received from Régine runs at once with the GO layer on.

Notes: use the same AGI identifiers on both sides (Araport11); Régine keeps only the genes of the DAP-seq network
and tells you which ones it dropped. For a GO term the background should ideally be the genes of the GeneCloud annotation
(*My own background* in Régine); a link parameter for that can be added later.
