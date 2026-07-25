#!/usr/bin/env python3
"""Figure 4A: ThermoMPNN core vs non-core violin plot.

Shows that EFS core positions are 3× more mutation-sensitive than non-core.
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# load thermompnn data (large file - use chunks for memory efficiency)
print("Loading ThermoMPNN data...")
thermo_df = pd.read_csv("/storage/kiran-stuff/RAA/PF00069/thermo_production_full/_family_all_data.csv",
                        usecols=['ddG_pred', 'position', 'wildtype', 'mutation', 'PDB_Name'])
print(f"Loaded {len(thermo_df):,} ThermoMPNN predictions")

# load core scores
print("Loading core scores...")
core_df = pd.read_csv("/storage/kiran-stuff/RAA/PF00069/family_core_scores.csv",
                      usecols=['protein_id', 'res_index', 'core_score'])
print(f"Loaded {len(core_df):,} core score entries")

# merge on protein_id and position (res_index = position + 1 in ThermoMPNN)
# ThermoMPNN uses 0-based positions, core_scores uses 1-based res_index
thermo_df['res_index'] = thermo_df['position'] + 1
merged = thermo_df.merge(core_df, left_on=['PDB_Name', 'res_index'],
                          right_on=['protein_id', 'res_index'], how='inner')
print(f"Merged: {len(merged):,} entries")

# label core vs non-core (C20 threshold = score >= 20)
merged['is_core'] = merged['core_score'] >= 20

# exclude self-mutations (wt -> wt)
merged = merged[merged['wildtype'] != merged['mutation']]
print(f"After excluding self-mutations: {len(merged):,}")

# calculate statistics
core_ddg = merged[merged['is_core']]['ddG_pred']
noncore_ddg = merged[~merged['is_core']]['ddG_pred']

print(f"\nCore positions: n={len(core_ddg):,}, mean ΔΔG = {core_ddg.mean():.3f}")
print(f"Non-core positions: n={len(noncore_ddg):,}, mean ΔΔG = {noncore_ddg.mean():.3f}")
print(f"Ratio: {core_ddg.mean() / noncore_ddg.mean():.1f}×")

# sample for plotting (too many points otherwise)
np.random.seed(42)
sample_size = 50000
core_sample = core_ddg.sample(min(sample_size, len(core_ddg)))
noncore_sample = noncore_ddg.sample(min(sample_size, len(noncore_ddg)))

# create violin plot
fig, ax = plt.subplots(figsize=(6, 7))

parts = ax.violinplot([noncore_sample, core_sample], positions=[1, 2],
                       showmeans=True, showmedians=True)

# color the violins
colors = ['#3498DB', '#E74C3C']  # blue for non-core, red for core
for i, pc in enumerate(parts['bodies']):
    pc.set_facecolor(colors[i])
    pc.set_alpha(0.7)
    pc.set_edgecolor('black')

# style the other parts
for partname in ('cbars', 'cmins', 'cmaxs', 'cmeans', 'cmedians'):
    if partname in parts:
        parts[partname].set_edgecolor('black')
        parts[partname].set_linewidth(1.5)

ax.set_xticks([1, 2])
ax.set_xticklabels(['Non-Core', 'EFS Core'], fontsize=12)
ax.set_ylabel('ΔΔG (kcal/mol)', fontsize=12)
ax.set_title('ThermoMPNN Mutation Sensitivity\nby Position Type', fontsize=14)

# add horizontal line at 0
ax.axhline(y=0, color='gray', linestyle='--', alpha=0.5)

# add statistics annotations
ax.text(1, noncore_ddg.mean() + 0.3, f'μ = {noncore_ddg.mean():.2f}',
        ha='center', fontsize=11, fontweight='bold', color='#2C3E50')
ax.text(2, core_ddg.mean() + 0.3, f'μ = {core_ddg.mean():.2f}',
        ha='center', fontsize=11, fontweight='bold', color='#C0392B')

# add sample size annotations
ax.text(1, ax.get_ylim()[0] + 0.2, f'n={len(noncore_ddg):,}',
        ha='center', fontsize=9, alpha=0.7)
ax.text(2, ax.get_ylim()[0] + 0.2, f'n={len(core_ddg):,}',
        ha='center', fontsize=9, alpha=0.7)

# add the 3x ratio annotation
ratio = core_ddg.mean() / noncore_ddg.mean()
ax.annotate(f'{ratio:.1f}× more\nsensitive', xy=(1.5, core_ddg.mean()),
            fontsize=12, fontweight='bold', color='#E74C3C', ha='center',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

plt.tight_layout()
plt.savefig('/storage/kiran-stuff/RAA/figures/fig4a_thermo_violin.png', dpi=300, bbox_inches='tight')
plt.savefig('/storage/kiran-stuff/RAA/figures/fig4a_thermo_violin.pdf', bbox_inches='tight')
print("\nSaved fig4a_thermo_violin.png/pdf")
