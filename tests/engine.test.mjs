// Régine 2026 — engine tests:  node --test tests/
import test from "node:test";
import assert from "node:assert/strict";
import fs from "fs";
import { parseGenes, chooseRatio, hyperUpper, bh, auc } from "../web/engine.js";
import { loadRegine } from "./load.mjs";

const { reg, meta } = loadRegine();
const nitrate = parseGenes(fs.readFileSync(new URL("../data/example_nitrate_list.txt", import.meta.url), "utf8"));

test("parseGenes: AGIs in any case and separator, isoforms, duplicates", () => {
  assert.deepEqual(parseGenes("At1g01010.1, AT1G01010; at5g67420\nfoo AT2G99999x"), ["AT1G01010", "AT5G67420"]);
  assert.equal(nitrate.length, 13);
});

test("data: 387 TFs, 2 848 928 distinct links", () => {
  assert.equal(meta.n_tfs, 387);
  let links = 0;
  for (let g = 0; g < reg.ids.length; g++) links += reg.nTF(g);
  assert.equal(links, meta.n_edges);
  assert.equal(links, 2848928);
});

test("chooseRatio, hypergeometric tail and BH", () => {
  assert.equal(chooseRatio(3, 5, 4), 0);
  assert.ok(Math.abs(chooseRatio(10, 13, 6) - 210 / 1716) < 1e-15);
  // R: phyper(9, 829, 32606 - 829, 13, lower.tail = FALSE)
  assert.ok(Math.abs(hyperUpper(10, 32606, 829, 13) / 2.856e-14 - 1) < 0.01);
  assert.ok(Math.abs(hyperUpper(0, 100, 10, 5) - 1) < 1e-12);
  assert.deepEqual(bh([0.01, 0.04, 0.03]).map(x => +x.toFixed(6)), [0.03, 0.04, 0.04]);
});

test("auc counts ties as one half", () => {
  assert.equal(auc([1, 1], [1, 1]), 0.5);
  assert.equal(auc([2, 3], [0, 1]), 1);
  assert.equal(auc([1], [0, 1, 2, 1]), 0.5);
});

test("real = fraction of all kernels of the list entirely bound (exact enumeration)", () => {
  const res = reg.run(nitrate, { reps: 200 });
  const k = 4, rows = reg.score(res, k), L = res.L, acc = new Uint32Array(reg.W);
  const counts = new Array(reg.tfs.length).fill(0);
  let total = 0;
  const rec = (start, chosen) => {
    if (chosen.length === k) {
      total++;
      acc.fill(0xffffffff);
      for (const r of chosen) for (let w = 0; w < reg.W; w++) acc[w] &= reg.mask[r * reg.W + w];
      for (let t = 0; t < reg.tfs.length; t++) if (acc[t >> 5] & (1 << (t & 31))) counts[t]++;
      return;
    }
    for (let i = start; i < L.length; i++) rec(i + 1, [...chosen, L[i]]);
  };
  rec(0, []);
  assert.equal(total, 715);
  for (const r of rows) assert.ok(Math.abs(r.real - counts[r.t] / total) < 1e-12);
});

test("nitrate example: NLP7 is the most specific regulator; automatic kernel 3..6", () => {
  const res = reg.run(nitrate, { reps: 3000 });
  const bySig = [...res.rows].sort((a, b) => a.p - b.p);
  assert.equal(bySig[0].tf.sym, "NLP7");
  assert.equal(bySig[0].kt, 10);
  assert.ok(bySig[0].q < 1e-10);
  assert.ok(res.kAuto >= 3 && res.kAuto <= 6);
  assert.ok(res.profile.every(p => p.auc > 0.5));
});

test("a random list does not depart from random genes", () => {
  const ids = reg.background("network"), pick = [];
  for (let i = 0; i < 60; i++) pick.push(reg.ids[ids[(i * 7919) % ids.length]]);
  const res = reg.run(pick, { reps: 3000, seed: 3 });
  assert.ok(Math.max(...res.profile.map(p => p.auc)) < 0.6);
});
