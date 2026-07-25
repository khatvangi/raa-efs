#!/usr/bin/env python3
"""Figure 1C: Amino acid composition of EFS-unique positions."""

import pandas as pd
import matplotlib.pyplot as plt
import glob
import os
from collections import Counter

core_dir = "/storage/kiran-stuff/RAA/ultimate_test_set/core_scores/"
csv_files = glob.glob(os.path.join(core_dir, "*_core.csv"))

efs_unique_aas = []

for f in csv_files:
    df = pd.read_csv(f)

    # EFS core but not SASA core
    efs_core = df['core_score'] >= 20
    sasa_core = df['sasa'] < 25
    efs_unique = efs_core & ~sasa_core

    # Get amino acids
    aas = df[efs_unique]['aa'].tolist()
    efs_unique_aas.extend(aas)

# Count
counts = Counter(efs_unique_aas)
total = sum(counts.values())

# Sort by count
sorted_aas = sorted(counts.items(), key=lambda x: -x[1])

# Plot
fig, ax = plt.subplots(figsize=(10, 5))

aas = [x[0] for x in sorted_aas]
values = [x[1] for x in sorted_aas]
percentages = [100*v/total for v in values]

# Color cysteines differently
colors = ['#E74C3C' if aa == 'C' else '#3498DB' for aa in aas]

bars = ax.bar(aas, percentages, color=colors, edgecolor='black', linewidth=0.5)

ax.set_xlabel('Amino Acid', fontsize=12)
ax.set_ylabel('Percentage of EFS-unique positions', fontsize=12)
ax.set_title(f'Amino Acid Composition of EFS-Unique Positions (n={total})', fontsize=14)

# Add percentage labels on bars
for bar, pct in zip(bars, percentages):
    if pct > 2:
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
                f'{pct:.1f}%', ha='center', va='bottom', fontsize=9)

# Highlight cysteine
cys_pct = 100 * counts.get('C', 0) / total
if 'C' in aas:
    cys_idx = aas.index('C')
    ax.annotate(f'Cysteine: {cys_pct:.1f}%\nof EFS-unique',
                xy=(cys_idx, percentages[cys_idx]),
                xytext=(min(cys_idx+3, len(aas)-2), max(percentages)*0.8),
                arrowprops=dict(arrowstyle='->', color='red'),
                fontsize=11, color='#E74C3C')

plt.tight_layout()
plt.savefig('/storage/kiran-stuff/RAA/figures/fig1c_aa_composition.png', dpi=300, bbox_inches='tight')
plt.savefig('/storage/kiran-stuff/RAA/figures/fig1c_aa_composition.pdf', bbox_inches='tight')
print(f"Saved fig1c_aa_composition.png/pdf")
print(f"Total EFS-unique positions: {total}")
print(f"Cysteine percentage: {cys_pct:.1f}%")
