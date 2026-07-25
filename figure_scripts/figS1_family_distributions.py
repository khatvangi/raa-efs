#!/usr/bin/env python3
"""Figure S1: Per-family TM-score distributions.

Shows distribution of TM-scores across all proteins in each family.
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

# load validation results
df = pd.read_csv("/storage/kiran-stuff/RAA/ultimate_test_set/validation_results.csv")

# get unique families
families = df['family'].unique()
n_families = len(families)

# create subplots
fig, axes = plt.subplots(2, 5, figsize=(15, 6))
axes = axes.flatten()

for i, family in enumerate(sorted(families)):
    ax = axes[i]
    family_data = df[df['family'] == family]

    # plot both EFS and Conservation
    x = np.arange(len(family_data))
    width = 0.35

    efs_vals = family_data['TM_EFS_reduced'].values
    cons_vals = family_data['TM_CONS_reduced'].values

    ax.bar(x - width/2, efs_vals, width, label='EFS', color='#E74C3C', alpha=0.8)
    ax.bar(x + width/2, cons_vals, width, label='Conservation', color='#3498DB', alpha=0.8)

    # add pass threshold line
    ax.axhline(y=0.7, color='gray', linestyle='--', alpha=0.5)

    ax.set_title(family, fontsize=10, fontweight='bold')
    ax.set_ylim(0, 1.1)
    ax.set_xticks([])

    if i == 0:
        ax.legend(fontsize=8)
    if i % 5 == 0:
        ax.set_ylabel('TM-score')

plt.suptitle('Per-Family TM-score Comparison (EFS vs Conservation)', fontsize=14, y=1.02)
plt.tight_layout()
plt.savefig('/storage/kiran-stuff/RAA/figures/figS1_family_distributions.png', dpi=300, bbox_inches='tight')
plt.savefig('/storage/kiran-stuff/RAA/figures/figS1_family_distributions.pdf', bbox_inches='tight')
print("Saved figS1_family_distributions.png/pdf")

# print summary
print("\n=== Family Summary ===")
for family in sorted(families):
    fam_data = df[df['family'] == family]
    efs_mean = fam_data['TM_EFS_reduced'].mean()
    cons_mean = fam_data['TM_CONS_reduced'].mean()
    delta = efs_mean - cons_mean
    winner = "EFS" if delta > 0 else "CONS"
    print(f"{family}: EFS={efs_mean:.3f}, CONS={cons_mean:.3f}, Δ={delta:+.3f} ({winner})")
