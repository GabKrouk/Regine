import fs from "fs";
import zlib from "zlib";
import { Regine } from "../web/engine.js";
const D = new URL("../web/data/", import.meta.url);
export function loadRegine() {
  const meta = JSON.parse(fs.readFileSync(new URL("meta.json", D)));
  const genes = JSON.parse(zlib.gunzipSync(fs.readFileSync(new URL("genes.json.gz", D))));
  const tfs = JSON.parse(zlib.gunzipSync(fs.readFileSync(new URL("tfs.json.gz", D))));
  const buf = zlib.gunzipSync(fs.readFileSync(new URL("network.bin.gz", D)));
  const mask = new Uint32Array(buf.buffer, buf.byteOffset, buf.byteLength / 4);
  return { meta, genes, tfs, reg: new Regine({ ids: genes.ids, flags: genes.flags, mask, W: meta.W, tfs }) };
}
