# Marker genes of each cell type (one cell type vs all other cells) and the expressed-gene background.
# Marker: expressed in >= 25% of the cells of the type, mean log difference >= 0.5, pct_in - pct_out >= 0.15; top 200 by
# mean log difference. Background: genes expressed in >= 1% of the cells of at least one cell type.
# Usage: python3 02_markers.py <celltype_stats.npz> <out_dir>
import sys, numpy as np, pandas as pd
z = np.load(sys.argv[1], allow_pickle=True); o = sys.argv[2].rstrip('/') + '/'
S, N, nc, cts, genes = z['S'], z['N'], z['ncell'], z['celltypes'], z['genes']
pct = N / nc[:, None]; St, Nt, ct = S.sum(0), N.sum(0), nc.sum()
a = pd.concat(pd.DataFrame({'celltype': c, 'gene': genes, 'mean_in': S[j] / nc[j], 'mean_out': (St - S[j]) / (ct - nc[j]),
                            'pct_in': pct[j], 'pct_out': (Nt - N[j]) / (ct - nc[j])}) for j, c in enumerate(cts))
a['logFC'] = a.mean_in - a.mean_out
a.to_csv(o + 'celltype_gene_stats.tsv.gz', sep='\t', index=False, float_format='%.4g')
open(o + 'background_expressed_genes.txt', 'w').write('\n'.join(genes[(pct >= 0.01).any(0)]) + '\n')
m = a[(a.pct_in >= 0.25) & (a.logFC >= 0.5) & (a.pct_in - a.pct_out >= 0.15)]
m = m.sort_values(['celltype', 'logFC'], ascending=[True, False]).groupby('celltype').head(200)
m.to_csv(o + 'celltype_markers_top200.tsv', sep='\t', index=False, float_format='%.4g')
