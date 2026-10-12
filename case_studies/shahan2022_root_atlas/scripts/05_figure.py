import pandas as pd, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.lines import Line2D

import sys
# Dot plot TF x cell type. Usage: python3 05_figure.py <results_dir>
S = O = sys.argv[1].rstrip('/') + '/'
a = pd.read_csv(S + 'regine_all_celltypes.tsv.gz', sep='\t')
prof = pd.read_csv(S + 'kernel_profiles.tsv', sep='\t')
auc = prof.groupby('celltype').auc.max()

cols = ['Columella', 'Lateral Root Cap', 'Atrichoblast', 'Trichoblast', 'Cortex', 'Endodermis', 'Quiescent Center',
        'Xylem Pole Pericycle', 'Phloem Pole Pericycle', 'Procambium', 'Protoxylem', 'Metaxylem', 'Protophloem',
        'Metaphloem & Companion Cell']
short = {'Lateral Root Cap': 'Lateral root cap', 'Quiescent Center': 'Quiescent center', 'Xylem Pole Pericycle': 'XPP',
         'Phloem Pole Pericycle': 'PPP', 'Metaphloem & Companion Cell': 'Metaphloem & CC'}
alias = {'NAC076': 'VND2', 'NAC105': 'VND3', 'NAC007': 'VND4', 'NAC101': 'VND6', 'DOF2.4': 'PEAR1', 'DOF5.1': 'PEAR2'}
known = {'Columella': ['SMB'], 'Protoxylem': ['VND1', 'NAC076', 'NAC105', 'NAC007', 'NAC101'],
         'Protophloem': ['DOF2.4', 'DOF5.1', 'DOF6', 'HCA2']}
reliable = [c for c in cols if auc[c] >= 0.65]

cand = a[(a.fdr < 0.05) & (a.tf_pct_in_celltype >= 0.10) & (a.specificity > 1)]
rows, group = [], []
for c in reliable:
    pick = list(known.get(c, []))
    for s in cand[cand.celltype == c].sort_values('specificity', ascending=False).symbol:
        if len([p for p in pick if p not in known.get(c, [])]) >= 3: break
        if s not in pick and s not in rows: pick.append(s)
    for s in pick:
        if s not in rows: rows.append(s); group.append(c)

M = a[a.symbol.isin(rows)].set_index(['symbol', 'celltype'])
fig = plt.figure(figsize=(8.2, 11.5))
gs = fig.add_gridspec(2, 1, height_ratios=[1.1, len(rows) * 0.33], hspace=0.04)
ax0 = fig.add_subplot(gs[0]); ax = fig.add_subplot(gs[1], sharex=ax0)
ink, muted, grid = '#0b0b0b', '#52514e', '#e4e3df'
blues = LinearSegmentedColormap.from_list('b', ['#cde2fb', '#86b6ef', '#3987e5', '#256abf', '#104281', '#0d366b'])
norm = Normalize(1.3, 20)

x = np.arange(len(cols))
ax0.bar(x, [auc[c] for c in cols], width=0.6, color=['#2a78d6' if auc[c] >= 0.65 else '#c9c8c3' for c in cols])
ax0.axhline(0.5, color=muted, lw=1, ls='--')
ax0.text(6.5, 0.505, 'random', ha='center', va='bottom', fontsize=7.5, color=muted)
ax0.set_ylim(0.4, 0.95); ax0.set_ylabel('Kernel AUC', fontsize=9, color=ink)
ax0.tick_params(axis='x', bottom=False, labelbottom=False); ax0.tick_params(axis='y', labelsize=8, colors=muted)
for s in ['top', 'right']: ax0.spines[s].set_visible(False)
ax0.set_title('Régine on Shahan et al. 2022 root atlas: candidate regulators of cell identity', fontsize=10.5, color=ink, loc='left', pad=10)

for i, r in enumerate(rows):
    for j, c in enumerate(cols):
        if (r, c) not in M.index: continue
        v = M.loc[(r, c)]
        size = 12 + 230 * float(v.tf_pct_in_celltype or 0)
        sig = v.fdr < 0.05
        col = blues(norm(min(-np.log10(v.fdr), 20))) if sig else '#ecebe7'
        ec = ink if (sig and v.tf_logFC_celltype > 0.25) else 'none'
        ax.scatter(j, i, s=size, c=[col], edgecolors=ec, linewidths=1.1, zorder=3)
ax.set_xticks(x); ax.set_xticklabels([short.get(c, c) for c in cols], rotation=45, ha='right', fontsize=8.5)
for t, c in zip(ax.get_xticklabels(), cols): t.set_color(ink if c in reliable else '#a3a29d')
ax.set_yticks(range(len(rows)))
ax.set_yticklabels([f'{alias[r]} ({r})' if r in alias else r for r in rows], fontsize=8)
for t, r in zip(ax.get_yticklabels(), rows):
    if any(r in v for v in known.values()): t.set_fontweight('bold')
ax.set_ylim(len(rows) - 0.5, -0.5); ax.set_xlim(-0.6, len(cols) - 0.4)
# group separators and the selecting cell type highlighted
for c in reliable:
    idx = [i for i, g in enumerate(group) if g == c]
    if not idx: continue
    j = cols.index(c)
    ax.add_patch(plt.Rectangle((j - 0.45, idx[0] - 0.45), 0.9, len(idx) - 0.1, fc='#f3f2ee', ec='none', zorder=1))
    ax.axhline(idx[-1] + 0.5, color=grid, lw=0.8, zorder=0)
for s in ['top', 'right']: ax.spines[s].set_visible(False)
ax.tick_params(colors=muted)

cax = fig.add_axes([0.80, 0.905, 0.12, 0.010])
cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=blues), cax=cax, orientation='horizontal')
cb.set_label('−log10 FDR (targets)', fontsize=7.5, color=ink); cb.ax.tick_params(labelsize=7, colors=muted); cb.outline.set_visible(False)
h = [Line2D([], [], marker='o', ls='', mfc='#86b6ef', mec='none', ms=np.sqrt(12 + 230 * p), label=f'{int(p*100)}%') for p in (0.1, 0.4, 0.8)]
h += [Line2D([], [], marker='o', ls='', mfc='#86b6ef', mec=ink, mew=1.1, ms=8, label='TF enriched\nin cell type'),
      Line2D([], [], marker='o', ls='', mfc='#ecebe7', mec='none', ms=8, label='FDR ≥ 0.05')]
leg = fig.legend(handles=h, title='Cells expressing TF', loc='upper left', bbox_to_anchor=(0.80, 0.855), fontsize=7.5,
                 title_fontsize=7.5, frameon=False, labelspacing=1.1, handletextpad=0.8)
fig.text(0.02, 0.012, 'Markers: top 200 genes per cell type (one vs rest). Régine: DAP-seq (O\'Malley 2016, 387 TFs), expressed-gene\nbackground, automatic kernel. '
         'Grey cell types: kernel AUC < 0.65, the list does not share DAP-seq regulators more than random genes.\nBold: known identity regulators. Shaded: cell type the TF was selected for.',
         fontsize=6.8, color=muted)
fig.subplots_adjust(left=0.2, right=0.77, top=0.95, bottom=0.135)
fig.savefig(O + 'Regine_Shahan2022_dotplot.png', dpi=300); fig.savefig(O + 'Regine_Shahan2022_dotplot.pdf')
print(len(rows), rows)
