# Per-cell-type mean log-normalised expression and fraction of expressing cells in the Shahan et al. 2022 root atlas.
# Input: the atlas exported from GSE152766 (Root_Atlas_spliced_unspliced_raw_counts.rds) as a column-compressed matrix
#   (genes x cells, int32 arrays counts_p / counts_i / counts_x, dim.txt, genes.txt) plus the cell metadata of
#   Supplemental Data 3 (Shahan2022_atlas_cell_metadata.tsv.gz, same cell order). Expression: log1p(counts / total UMIs x 1e4).
# Usage: python3 01_celltype_stats.py <atlas_dir> <out.npz>
import sys, numpy as np, pandas as pd
d, out = sys.argv[1].rstrip('/') + '/', sys.argv[2]
G, C = map(int, open(d + 'derived/dim.txt').read().split())
genes = open(d + 'derived/genes.txt').read().split()
ct = pd.Categorical(pd.read_csv(d + 'metadata/Shahan2022_atlas_cell_metadata.tsv.gz', sep='\t', usecols=['celltype.anno'])['celltype.anno'])
codes = ct.codes.astype(np.int64); T = len(ct.categories)
p = np.fromfile(d + 'derived/counts_p.int32', dtype=np.int32).astype(np.int64)
I = np.memmap(d + 'derived/counts_i.int32', dtype=np.int32, mode='r'); X = np.memmap(d + 'derived/counts_x.int32', dtype=np.int32, mode='r')
S = np.zeros(T * G); N = np.zeros(T * G)
for a in range(0, C, 5000):
    b = min(C, a + 5000); i = I[p[a]:p[b]].astype(np.int64); x = X[p[a]:p[b]].astype(np.float64)
    cell = np.repeat(np.arange(a, b), np.diff(p[a:b + 1]))
    lib = np.bincount(cell - a, weights=x, minlength=b - a)
    key = codes[cell] * G + i
    S += np.bincount(key, weights=np.log1p(x / lib[cell - a] * 1e4), minlength=T * G); N += np.bincount(key, minlength=T * G)
np.savez_compressed(out, S=S.reshape(T, G), N=N.reshape(T, G), ncell=np.bincount(codes, minlength=T),
                    celltypes=np.array(ct.categories), genes=np.array(genes))
