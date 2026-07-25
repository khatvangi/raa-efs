#!/usr/bin/env python3
"""Figure 2A: Scatter plot EFS vs Conservation TM-scores."""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Load validation results
df = pd.read_csv("/storage/kiran-stuff/RAA/ultimate_test_set/validation_results.csv")

fig, ax = plt.subplots(figsize=(7, 7))

# Scatter plot
ax.scatter(df['TM_CONS_reduced'], df['TM_EFS_reduced'],
           s=80, alpha=0.7, c='#2C3E50', edgecolors='white', linewidth=0.5)

# Unity line
lims = [0, 1]
ax.plot(lims, lims, 'k--', alpha=0.5, linewidth=1, label='y=x')

# Fill regions
ax.fill_between([0, 1], [0, 1], [1, 1], alpha=0.1, color='green', label='EFS wins')
ax.fill_between([0, 1], [0, 0], [0, 1], alpha=0.1, color='red', label='Conservation wins')

# Labels for each point (family prefix)
for idx, row in df.iterrows():
    ax.annotate(row['family'],
                (row['TM_CONS_reduced'], row['TM_EFS_reduced']),
                fontsize=7, alpha=0.6,
                xytext=(3, 3), textcoords='offset points')

ax.set_xlabel('TM-score (Conservation-reduced)', fontsize=12)
ax.set_ylabel('TM-score (EFS-reduced)', fontsize=12)
ax.set_title('Folding Validation: EFS vs Conservation', fontsize=14)
ax.set_xlim(0, 1.05)
ax.set_ylim(0, 1.05)
ax.set_aspect('equal')
ax.legend(loc='lower right')

# Add statistics
n_efs_wins = (df['TM_EFS_reduced'] > df['TM_CONS_reduced']).sum()
n_total = len(df)
mean_delta = (df['TM_EFS_reduced'] - df['TM_CONS_reduced']).mean()
ax.text(0.05, 0.95, f'EFS wins: {n_efs_wins}/{n_total} ({100*n_efs_wins/n_total:.0f}%)\nMean Δ: {mean_delta:+.3f}',
        transform=ax.transAxes, fontsize=11, verticalalignment='top',
        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

plt.tight_layout()
plt.savefig('/storage/kiran-stuff/RAA/figures/fig2a_scatter.png', dpi=300, bbox_inches='tight')
plt.savefig('/storage/kiran-stuff/RAA/figures/fig2a_scatter.pdf', bbox_inches='tight')
print("Saved fig2a_scatter.png/pdf")
