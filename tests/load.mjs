import fs from "fs";
import zlib from "zlib";
import { Regine, targetNetworks } from "../web/engine.js";
const D = new URL("../web/data/", import.meta.url);
export function loadRegine() {
  const meta = JSON.parse(fs.readFileSync(new URL("meta.json", D)));
  const genes = JSON.parse(zlib.gunzipSync(fs.readFileSync(new URL("genes.json.gz", D))));
  const tfs = JSON.parse(zlib.gunzipSync(fs.readFileSync(new URL("tfs.json.gz", D))));
  const buf = zlib.gunzipSync(fs.readFileSync(new URL("network.bin.gz", D)));
  const mask = new Uint32Array(buf.buffer, buf.byteOffset, buf.byteLength / 4);
  const reg = new Regine({ ids: genes.ids, flags: genes.flags, mask, W: meta.W, tfs });
  const tmeta = JSON.parse(zlib.gunzipSync(fs.readFileSync(new URL("target.json.gz", D))));
  const tbuf = zlib.gunzipSync(fs.readFileSync(new URL("target.bin.gz", D)));
  return { meta, genes, tfs, reg, tmeta, target: targetNetworks(reg, new Uint8Array(tbuf), tmeta) };
}
