#!/usr/bin/env python3
"""Figure 2B: Per-family bar chart of EFS vs Conservation."""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Load validation results
df = pd.read_csv("/storage/kiran-stuff/RAA/ultimate_test_set/validation_results.csv")

# Aggregate by family
family_stats = df.groupby('family').agg({
    'TM_EFS_reduced': 'mean',
    'TM_CONS_reduced': 'mean'
}).reset_index()

# Sort by EFS TM-score
family_stats = family_stats.sort_values('TM_EFS_reduced', ascending=True)

fig, ax = plt.subplots(figsize=(10, 6))

y = np.arange(len(family_stats))
height = 0.35

bars1 = ax.barh(y - height/2, family_stats['TM_EFS_reduced'], height,
                label='EFS-reduced', color='#E74C3C', edgecolor='black', linewidth=0.5)
bars2 = ax.barh(y + height/2, family_stats['TM_CONS_reduced'], height,
                label='Conservation-reduced', color='#3498DB', edgecolor='black', linewidth=0.5)

ax.set_yticks(y)
ax.set_yticklabels(family_stats['family'])
ax.set_xlabel('Mean TM-score', fontsize=12)
ax.set_ylabel('Protein Family', fontsize=12)
ax.set_title('EFS vs Conservation by Family', fontsize=14)
ax.legend(loc='lower right')
ax.set_xlim(0, 1.1)

# Add delta values
for i, (_, row) in enumerate(family_stats.iterrows()):
    delta = row['TM_EFS_reduced'] - row['TM_CONS_reduced']
    color = 'green' if delta > 0 else 'red'
    ax.text(max(row['TM_EFS_reduced'], row['TM_CONS_reduced']) + 0.02, i,
            f'Δ={delta:+.3f}', va='center', fontsize=9, color=color)

# Add vertical line at 0.7 (pass threshold)
ax.axvline(x=0.7, color='gray', linestyle='--', linewidth=1, alpha=0.7)
ax.text(0.71, len(family_stats)-0.5, 'Pass threshold', fontsize=9, alpha=0.7)

plt.tight_layout()
plt.savefig('/storage/kiran-stuff/RAA/figures/fig2b_family_bars.png', dpi=300, bbox_inches='tight')
plt.savefig('/storage/kiran-stuff/RAA/figures/fig2b_family_bars.pdf', bbox_inches='tight')
print("Saved fig2b_family_bars.png/pdf")
