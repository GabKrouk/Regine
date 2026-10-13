# Adds TF expression in each cell type and a cell-type specificity, and selects the candidates.
# specificity = -log10 FDR in this cell type - median -log10 FDR over the 13 other cell types.
# candidates: FDR < 0.05, TF expressed in >= 10% of the cells of the type, specificity > 1.
# Usage: python3 04_annotate.py <dir with regine_all.tsv>  (TF expression from ../data/tf_expression_by_celltype.tsv)
import sys, os, numpy as np, pandas as pd
d = sys.argv[1].rstrip('/') + '/'; here = os.path.dirname(os.path.abspath(__file__))
r = pd.read_csv(d + 'regine_all.tsv', sep='\t')
st = pd.read_csv(os.path.join(here, '../data/tf_expression_by_celltype.tsv'), sep='\t')
mlq = -np.log10(r.fdr.clip(lower=1e-300)); r['_m'] = mlq
piv = r.pivot(index='tf_agi', columns='celltype', values='_m')
r['median_mlq_other'] = [piv.loc[t].drop(c).median() for t, c in zip(r.tf_agi, r.celltype)]
r['specificity'] = r._m - r.median_mlq_other
r = r.drop(columns='_m').merge(st.rename(columns={'gene': 'tf_agi', 'pct_in': 'tf_pct_in_celltype', 'pct_out': 'tf_pct_other',
                                                  'logFC': 'tf_logFC_celltype'}), on=['celltype', 'tf_agi'], how='left')
r.to_csv(d + 'regine_all_celltypes.tsv.gz', sep='\t', index=False, float_format='%.4g')
c = r[(r.fdr < 0.05) & (r.tf_pct_in_celltype >= 0.10) & (r.specificity > 1)].sort_values(['celltype', 'specificity'], ascending=[True, False])
c.to_csv(d + 'candidates.tsv', sep='\t', index=False, float_format='%.4g')
