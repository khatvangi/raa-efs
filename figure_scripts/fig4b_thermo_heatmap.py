#!/usr/bin/env python3
"""Figure 4B: ThermoMPNN substitution matrix heatmap.

Shows average ΔΔG for each amino acid substitution type in EFS core positions.
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

# load thermompnn data
print("Loading ThermoMPNN data...")
thermo_df = pd.read_csv("/storage/kiran-stuff/RAA/PF00069/thermo_production_full/_family_all_data.csv",
                        usecols=['ddG_pred', 'position', 'wildtype', 'mutation', 'PDB_Name'])
print(f"Loaded {len(thermo_df):,} ThermoMPNN predictions")

# load core scores
print("Loading core scores...")
core_df = pd.read_csv("/storage/kiran-stuff/RAA/PF00069/family_core_scores.csv",
                      usecols=['protein_id', 'res_index', 'core_score'])

# merge
thermo_df['res_index'] = thermo_df['position'] + 1
merged = thermo_df.merge(core_df, left_on=['PDB_Name', 'res_index'],
                          right_on=['protein_id', 'res_index'], how='inner')

# filter to core positions only
core_data = merged[merged['core_score'] >= 20].copy()
print(f"Core position mutations: {len(core_data):,}")

# exclude self-mutations
core_data = core_data[core_data['wildtype'] != core_data['mutation']]

# create substitution matrix
aa_order = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
            'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

# group by wildtype and mutation, compute mean ΔΔG
sub_matrix = core_data.groupby(['wildtype', 'mutation'])['ddG_pred'].mean().unstack()

# reorder to standard AA order (filter to only those present)
present_aas = [aa for aa in aa_order if aa in sub_matrix.index and aa in sub_matrix.columns]
sub_matrix = sub_matrix.loc[present_aas, present_aas]

# create heatmap
fig, ax = plt.subplots(figsize=(12, 10))

# use diverging colormap centered at 0
vmax = max(abs(sub_matrix.min().min()), abs(sub_matrix.max().max()))
vmax = min(vmax, 3.0)  # cap at ±3 for visibility

sns.heatmap(sub_matrix, ax=ax, cmap='RdBu_r', center=0,
            vmin=-vmax, vmax=vmax, annot=True, fmt='.1f',
            annot_kws={'fontsize': 7}, cbar_kws={'label': 'ΔΔG (kcal/mol)'})

ax.set_xlabel('Mutant Amino Acid', fontsize=12)
ax.set_ylabel('Wild-type Amino Acid', fontsize=12)
ax.set_title('ThermoMPNN Substitution Matrix\n(EFS Core Positions Only)', fontsize=14)

# highlight RAA-relevant substitutions
# draw boxes around hydrophobic → hydrophobic (L,I,V,M,F,W,Y)
hydrophobic_aas = ['L', 'I', 'V', 'M', 'F', 'W', 'Y', 'A']

plt.tight_layout()
plt.savefig('/storage/kiran-stuff/RAA/figures/fig4b_thermo_heatmap.png', dpi=300, bbox_inches='tight')
plt.savefig('/storage/kiran-stuff/RAA/figures/fig4b_thermo_heatmap.pdf', bbox_inches='tight')
print("\nSaved fig4b_thermo_heatmap.png/pdf")

# also print RAA-relevant statistics
print("\n=== RAA Substitution Analysis (Core Positions) ===")
# within hydrophobic group
hydro_wt = core_data[core_data['wildtype'].isin(['L', 'I', 'V', 'M'])]
hydro_to_hydro = hydro_wt[hydro_wt['mutation'].isin(['L', 'I', 'V', 'M'])]
hydro_to_charged = hydro_wt[hydro_wt['mutation'].isin(['D', 'E', 'K', 'R'])]

print(f"Hydrophobic → Hydrophobic (I,L,V,M): ΔΔG = {hydro_to_hydro['ddG_pred'].mean():.3f}")
print(f"Hydrophobic → Charged (D,E,K,R): ΔΔG = {hydro_to_charged['ddG_pred'].mean():.3f}")

# cysteine analysis
cys_wt = core_data[core_data['wildtype'] == 'C']
if len(cys_wt) > 0:
    cys_to_hydro = cys_wt[cys_wt['mutation'].isin(['L', 'I', 'V', 'M', 'F', 'A'])]
    cys_to_charged = cys_wt[cys_wt['mutation'].isin(['D', 'E', 'K', 'R'])]
    print(f"Cysteine → Hydrophobic: ΔΔG = {cys_to_hydro['ddG_pred'].mean():.3f}")
    print(f"Cysteine → Charged: ΔΔG = {cys_to_charged['ddG_pred'].mean():.3f}")
