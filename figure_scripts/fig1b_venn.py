#!/usr/bin/env python3
"""Figure 1B: EFS vs SASA core-definition Venn diagram.

Uses the PER-PROTEIN MEAN Jaccard (0.865, Table S6) rather than the pooled
value. Displayed Jaccard is the mean of each protein's own EFS-vs-SASA overlap;
position counts (shared / EFS-unique / SASA-unique) are summed across proteins
and are internally consistent with the per-protein framing.

source of truth: ultimate_test_set/core_definition_comparison.csv
"""

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

df = pd.read_csv("/storage/kiran-stuff/RAA/ultimate_test_set/core_definition_comparison.csv")

# per-protein mean Jaccard (Table S6 framing)
jaccard = df["efs_vs_sasa"].mean()

# position counts, per protein then summed (consistent with per-protein overlap)
efs_unique = int(df["efs_unique"].sum())                 # 246 EFS-unique positions
shared = int((df["efs_n"] - df["efs_unique"]).sum())     # positions in both
sasa_unique = int((df["sasa_n"] - (df["efs_n"] - df["efs_unique"])).sum())

fig, ax = plt.subplots(figsize=(6, 6))

v = venn2(subsets=(efs_unique, sasa_unique, shared),
          set_labels=("EFS Core", "SASA Core"),
          set_colors=("#E74C3C", "#3498DB"),
          alpha=0.7,
          ax=ax)

ax.set_title(f"Core Definition Overlap\nJaccard = {jaccard:.3f}", fontsize=14)

plt.tight_layout()
plt.savefig("/storage/kiran-stuff/RAA/figures/fig1b_venn.png", dpi=300, bbox_inches="tight")
plt.savefig("/storage/kiran-stuff/RAA/figures/fig1b_venn.pdf", bbox_inches="tight")
print("Saved fig1b_venn.png/pdf")
print(f"EFS-unique: {efs_unique}, SASA-unique: {sasa_unique}, Shared: {shared}")
print(f"Jaccard (per-protein mean): {jaccard:.3f}")
