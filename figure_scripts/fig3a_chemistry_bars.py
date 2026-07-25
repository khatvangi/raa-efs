#!/usr/bin/env python3
"""Figure 3A: Chemistry inversion experiment bar plot."""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Load results
df = pd.read_csv("/storage/kiran-stuff/RAA/chemistry_inversion/usalign_results.csv")

# Aggregate by condition
condition_stats = df.groupby('condition')['TM_score'].agg(['mean', 'std', 'count']).reset_index()

# Define order and colors
condition_order = ['NATIVE', 'EFS_PRESERVED', 'ANTI_EFS_REDUCED', 'BLIND_REDUCED', 'ANTI_EFS_INVERTED']
colors = ['#2ECC71', '#27AE60', '#F39C12', '#E67E22', '#E74C3C']

# Reorder
condition_stats['condition'] = pd.Categorical(condition_stats['condition'],
                                               categories=condition_order, ordered=True)
condition_stats = condition_stats.sort_values('condition')

fig, ax = plt.subplots(figsize=(10, 6))

x = np.arange(len(condition_stats))
bars = ax.bar(x, condition_stats['mean'],
              yerr=condition_stats['std'], capsize=5,
              color=colors, edgecolor='black', linewidth=1)

ax.set_xticks(x)
ax.set_xticklabels(['Native', 'EFS Preserved', 'Anti-EFS\nReduced', 'Blind\nReduced', 'Anti-EFS\nInverted'],
                   fontsize=10)
ax.set_ylabel('Mean TM-score', fontsize=12)
ax.set_title('Chemistry Inversion Experiment\n(Same positions, different chemistry)', fontsize=14)
ax.set_ylim(0, 1.15)

# Add value labels
for bar, (_, row) in zip(bars, condition_stats.iterrows()):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + row['std'] + 0.02,
            f'{row["mean"]:.2f}', ha='center', va='bottom', fontsize=11, fontweight='bold')

# Add horizontal line at 0.7
ax.axhline(y=0.7, color='gray', linestyle='--', linewidth=1, alpha=0.7)
ax.text(4.5, 0.72, 'Pass threshold', fontsize=9, alpha=0.7)

# Highlight the key comparison
anti_efs_reduced_mean = condition_stats[condition_stats['condition'] == 'ANTI_EFS_REDUCED']['mean'].values[0]
anti_efs_inverted_mean = condition_stats[condition_stats['condition'] == 'ANTI_EFS_INVERTED']['mean'].values[0]
delta = anti_efs_reduced_mean - anti_efs_inverted_mean

ax.annotate('', xy=(2, anti_efs_reduced_mean - 0.05), xytext=(4, anti_efs_inverted_mean + 0.05),
            arrowprops=dict(arrowstyle='<->', color='red', lw=2))
ax.text(3, (anti_efs_reduced_mean + anti_efs_inverted_mean)/2,
        f'Δ = {delta:.2f}\nSame positions,\ndifferent chemistry',
        fontsize=10, ha='center', color='red', fontweight='bold')

plt.tight_layout()
plt.savefig('/storage/kiran-stuff/RAA/figures/fig3a_chemistry_bars.png', dpi=300, bbox_inches='tight')
plt.savefig('/storage/kiran-stuff/RAA/figures/fig3a_chemistry_bars.pdf', bbox_inches='tight')
print("Saved fig3a_chemistry_bars.png/pdf")
