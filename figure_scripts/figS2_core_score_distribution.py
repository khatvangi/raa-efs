#!/usr/bin/env python3
"""Figure S2: Core score distributions across families.

Shows the distribution of EFS core scores and what fraction of positions
are classified as core across different protein families.
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import glob
import os

core_dir = "/storage/kiran-stuff/RAA/ultimate_test_set/core_scores/"
csv_files = glob.glob(os.path.join(core_dir, "*_core.csv"))

# collect data by family
family_data = {}

for f in csv_files:
    filename = os.path.basename(f)
    family = filename.split('_')[0]

    df = pd.read_csv(f)

    if family not in family_data:
        family_data[family] = {'core_scores': [], 'n_core': 0, 'n_total': 0}

    family_data[family]['core_scores'].extend(df['core_score'].tolist())
    family_data[family]['n_core'] += (df['core_score'] >= 20).sum()
    family_data[family]['n_total'] += len(df)

# create figure with two panels
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

# panel 1: core score distributions
families = sorted(family_data.keys())
positions = np.arange(len(families))

for i, family in enumerate(families):
    scores = family_data[family]['core_scores']
    # violin plot
    parts = ax1.violinplot([scores], positions=[i], showmeans=True, showmedians=True)
    for pc in parts['bodies']:
        pc.set_facecolor(plt.cm.tab10(i))
        pc.set_alpha(0.7)

ax1.axhline(y=20, color='red', linestyle='--', linewidth=2, label='Core threshold (≥20)')
ax1.set_xticks(positions)
ax1.set_xticklabels(families, rotation=45, ha='right', fontsize=9)
ax1.set_ylabel('Core Score', fontsize=12)
ax1.set_title('Core Score Distribution by Family', fontsize=12)
ax1.legend(loc='upper right')

# panel 2: fraction core per family
fractions = [family_data[f]['n_core'] / family_data[f]['n_total'] * 100 for f in families]
colors = plt.cm.tab10(np.linspace(0, 1, len(families)))

bars = ax2.bar(positions, fractions, color=colors, edgecolor='black', linewidth=0.5)

# add value labels
for bar, frac in zip(bars, fractions):
    ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1,
             f'{frac:.0f}%', ha='center', va='bottom', fontsize=9)

ax2.set_xticks(positions)
ax2.set_xticklabels(families, rotation=45, ha='right', fontsize=9)
ax2.set_ylabel('Percentage Core Positions', fontsize=12)
ax2.set_title('Fraction of Positions in EFS Core', fontsize=12)
ax2.set_ylim(0, max(fractions) * 1.2)

plt.tight_layout()
plt.savefig('/storage/kiran-stuff/RAA/figures/figS2_core_score_distribution.png', dpi=300, bbox_inches='tight')
plt.savefig('/storage/kiran-stuff/RAA/figures/figS2_core_score_distribution.pdf', bbox_inches='tight')
print("Saved figS2_core_score_distribution.png/pdf")

# print summary
print("\n=== Core Fraction Summary ===")
for family, frac in zip(families, fractions):
    n_core = family_data[family]['n_core']
    n_total = family_data[family]['n_total']
    print(f"{family}: {n_core}/{n_total} core ({frac:.1f}%)")
