// Régine 2026 — computation engine (browser and Node). No dependency.
//
// Network: gene-major bit matrix, W uint32 words per gene; bit t of gene g = TF t binds g (DAP-seq).
//
// For a gene list L (n genes) and a background B (N genes, L ⊂ B), for each TF t:
//   k_t = genes of L bound by t,  K_t = genes of B bound by t
//   real_t(k)  = C(k_t, k) / C(n, k)   probability that a random kernel of k genes of the list is entirely bound by t
//   noise_t(k) = C(K_t, k) / C(N, k)   same for a random kernel of k background genes
//   score_t(k) = real_t − noise_t      (the "Real − Noise" of Régine 2022, as an exact frequency instead of a count)
// plus the hypergeometric P(X ≥ k_t) of the overlap and its Benjamini–Hochberg FDR over all TFs.
//
// Kernel profile: for every kernel size k, `reps` random kernels of the list and of the background; the number of
// TFs binding all genes of a kernel is recorded. The separation of the two distributions is measured by the AUC
// (probability that a list kernel has more common regulators than a background kernel, ties count one half).
// The automatic kernel size is the one with the largest AUC.

export const AGI_RE = /\bAT[1-5CM]G\d{5}\b/gi;

export function parseGenes(text) {
  const ids = [], seen = new Set();
  for (const m of String(text).matchAll(AGI_RE)) {
    const g = m[0].toUpperCase();
    if (!seen.has(g)) { seen.add(g); ids.push(g); }
  }
  return ids;
}

function popcount(x) {
  x -= (x >>> 1) & 0x55555555;
  x = (x & 0x33333333) + ((x >>> 2) & 0x33333333);
  return (((x + (x >>> 4)) & 0x0f0f0f0f) * 0x01010101) >>> 24;
}

// mulberry32: small seeded generator, so that the same input gives the same profile
export function rng(seed) {
  let s = seed >>> 0;
  return () => {
    s = (s + 0x6d2b79f5) | 0;
    let t = s;
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

// ratio C(a, k) / C(b, k) = prod_{i<k} (a − i) / (b − i)
export function chooseRatio(a, b, k) {
  if (a < k) return 0;
  let r = 1;
  for (let i = 0; i < k; i++) r *= (a - i) / (b - i);
  return r;
}

let LF = new Float64Array([0]);
function logFact(n) {
  if (n >= LF.length) {
    const m = Math.max(n + 1, LF.length * 2), next = new Float64Array(m);
    next.set(LF);
    for (let i = LF.length; i < m; i++) next[i] = next[i - 1] + Math.log(i);
    LF = next;
  }
  return LF[n];
}
const lchoose = (n, k) => logFact(n) - logFact(k) - logFact(n - k);

// P(X ≥ x), X ~ hypergeometric(N population, K successes, n draws)
export function hyperUpper(x, N, K, n) {
  const lo = Math.max(x, 0, n - (N - K)), hi = Math.min(n, K);
  if (lo > hi) return 0;
  const base = lchoose(N, n);
  let m = -Infinity;
  const terms = [];
  for (let i = lo; i <= hi; i++) { const v = lchoose(K, i) + lchoose(N - K, n - i) - base; terms.push(v); if (v > m) m = v; }
  let s = 0;
  for (const v of terms) s += Math.exp(v - m);
  return Math.min(1, Math.exp(m + Math.log(s)));
}

export function bh(p) {
  const n = p.length, o = p.map((v, i) => i).sort((a, b) => p[b] - p[a]), q = new Array(n);
  let min = 1;
  o.forEach((i, r) => { min = Math.min(min, (p[i] * n) / (n - r)); q[i] = min; });
  return q;
}

export class Regine {
  // ids: gene identifiers (row order of mask); flags: 1 protein coding, 2 in network; tfs: [{id, ...}] in bit order
  constructor({ ids, flags, mask, W, tfs }) {
    Object.assign(this, { ids, flags, mask, W, tfs });
    this.index = new Map(ids.map((g, i) => [g, i]));
  }

  background(kind, custom = []) {
    if (kind === "custom") return [...new Set(custom.map(g => this.index.get(g)).filter(i => i != null))];
    const bit = kind === "protein_coding" ? 1 : 2;
    const out = [];
    this.flags.forEach((f, i) => { if (f & bit) out.push(i); });
    return out;
  }

  // number of TFs binding every gene of the kernel (rows)
  common(rows, k, acc) {
    const { mask, W } = this;
    acc.fill(0xffffffff);
    for (let j = 0; j < k; j++) { const o = rows[j] * W; for (let w = 0; w < W; w++) acc[w] &= mask[o + w]; }
    let c = 0;
    for (let w = 0; w < W; w++) c += popcount(acc[w]);
    return c;
  }

  // counts[t] = genes of rows bound by TF t
  coverage(rows) {
    const { mask, W } = this, nT = this.tfs.length, c = new Array(nT).fill(0);
    for (const r of rows) {
      const o = r * W;
      for (let w = 0; w < W; w++) {
        let x = mask[o + w];
        while (x) { const b = 31 - Math.clz32(x & -x); c[w * 32 + b]++; x &= x - 1; }
      }
    }
    return c;
  }

  // genes: AGIs of the list; opts: {background: "network"|"protein_coding"|"custom", custom: [AGI], kmin, kmax, reps, seed}
  run(genes, opts = {}) {
    const { background = "network", custom = [], kmin = 2, reps = 1000, seed = 1 } = opts;
    const unknown = genes.filter(g => !this.index.has(g));
    const known = genes.filter(g => this.index.has(g));
    const L = known.map(g => this.index.get(g));
    const notInNetwork = known.filter(g => !(this.flags[this.index.get(g)] & 2));
    const B0 = this.background(background, custom), inB = new Set(B0);
    const added = L.filter(r => !inB.has(r));
    const B = B0.concat(added);
    const n = L.length, N = B.length;
    if (n < 2) throw new Error("At least 2 recognised genes are needed.");
    const kmax = Math.min(opts.kmax ?? 10, n);

    // kernel profile (Monte Carlo)
    const rand = rng(seed), acc = new Uint32Array(this.W);
    const pl = L.slice(), pb = B.slice();
    const draw = (pool, k) => {             // partial Fisher–Yates; the pool stays a permutation, so draws stay uniform
      for (let i = 0; i < k; i++) { const j = i + Math.floor(rand() * (pool.length - i)); const t = pool[i]; pool[i] = pool[j]; pool[j] = t; }
      return pool;
    };
    const cl = this.coverage(L), cb = this.coverage(B);
    const profile = [];
    for (let k = kmin; k <= kmax; k++) {
      const xl = new Int32Array(reps), xb = new Int32Array(reps);
      for (let r = 0; r < reps; r++) { xl[r] = this.common(draw(pl, k), k, acc); xb[r] = this.common(draw(pb, k), k, acc); }
      profile.push({
        k, list: xl, bg: xb, auc: auc(xl, xb),
        meanList: cl.reduce((s, kt) => s + chooseRatio(kt, n, k), 0),   // exact expectations
        meanBg: cb.reduce((s, Kt) => s + chooseRatio(Kt, N, k), 0),
      });
    }
    let best = profile[0];
    for (const p of profile) if (p.auc > best.auc + 1e-12) best = p;

    const p = cl.map((kt, t) => hyperUpper(kt, N, cb[t], n)), q = bh(p);
    const res = { genes, unknown, notInNetwork, added: added.length, n, N, background, reps, seed, profile, kAuto: best.k, cl, cb, p, q, L };
    res.k = opts.kernel && opts.kernel !== "auto" ? Math.min(+opts.kernel, kmax) : best.k;
    res.rows = this.score(res, res.k);
    return res;
  }

  // TF table for kernel size k (exact; no new sampling)
  score(res, k) {
    const { n, N, cl, cb, p, q } = res;
    const rows = this.tfs.map((tf, t) => {
      const real = chooseRatio(cl[t], n, k), noise = chooseRatio(cb[t], N, k);
      return { t, tf, kt: cl[t], Kt: cb[t], real, noise, score: real - noise,
        fold: cl[t] / n / (cb[t] / N || Infinity), p: p[t], q: q[t] };
    });
    rows.sort((a, b) => b.score - a.score || a.p - b.p || b.kt - a.kt);
    rows.forEach((r, i) => (r.rank = i + 1));
    return rows;
  }

  // genes of the list bound by TF t
  targetsOf(res, t) {
    const w = t >> 5, b = 1 << (t & 31);
    return res.L.filter(r => this.mask[r * this.W + w] & b).map(r => this.ids[r]);
  }

  // number of TFs binding each gene
  nTF(row) {
    let c = 0;
    for (let w = 0; w < this.W; w++) c += popcount(this.mask[row * this.W + w]);
    return c;
  }
}

// AUC = P(X > Y) + P(X = Y) / 2, by counting (values are small integers)
export function auc(x, y) {
  let max = 0;
  for (const v of x) if (v > max) max = v;
  for (const v of y) if (v > max) max = v;
  const hy = new Float64Array(max + 2);
  for (const v of y) hy[v + 1]++;
  const cum = new Float64Array(max + 2);           // cum[v] = #{y < v}
  for (let v = 1; v <= max + 1; v++) cum[v] = cum[v - 1] + hy[v];
  let s = 0;
  for (const v of x) s += cum[v] + hy[v + 1] / 2;
  return s / (x.length * y.length);
}

export function quantiles(a) {
  const s = Array.from(a).sort((x, y) => x - y), q = f => {
    const h = (s.length - 1) * f, l = Math.floor(h);
    return s[l] + (h - l) * ((s[l + 1] ?? s[l]) - s[l]);
  };
  const q1 = q(0.25), q3 = q(0.75), iqr = q3 - q1;
  const lo = s.find(v => v >= q1 - 1.5 * iqr), hi = [...s].reverse().find(v => v <= q3 + 1.5 * iqr);
  return { min: s[0], q1, med: q(0.5), q3, max: s[s.length - 1], lo, hi, mean: s.reduce((a, b) => a + b, 0) / s.length };
}
