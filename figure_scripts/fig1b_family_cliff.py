#!/usr/bin/env python3
"""Figure 1 panel b: per-family cliff decomposition.

Mean TM vs reduced alphabet (RAA13 -> RAA5) with one line per family.
Point of the panel: the RAA10 -> RAA9 drop is concentrated in PF00089 trypsin
(-0.409) while kinase (PF00069) stays nearly flat.

data: FINAL_RESULTS/aggregate/raa89_cliff_summary.csv (per family x alphabet mean TM)
"""

import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("/storage/kiran-stuff/RAA/FINAL_RESULTS/aggregate/raa89_cliff_summary.csv")

# x-axis: coarser to the right
alphabets = ["RAA13", "RAA10", "RAA9", "RAA8", "RAA7", "RAA5"]
cols = [f"{a}_meanTM" for a in alphabets]
x = list(range(len(alphabets)))

families = {
    "PF00069": ("PF00069 kinase", "#E74C3C"),
    "PF00089": ("PF00089 trypsin", "#3498DB"),
    "PF01145": ("PF01145 ferredoxin", "#27AE60"),
    "PF07686": ("PF07686 spectrin", "#F39C12"),
    "PF12796": ("PF12796 ankyrin", "#8E44AD"),
}

fig, ax = plt.subplots(figsize=(7, 6))

for fam, (label, color) in families.items():
    row = df[df["family"] == fam]
    if row.empty:
        continue
    y = row[cols].values.flatten().astype(float)
    ax.plot(x, y, marker="o", markersize=6, linewidth=2, color=color,
            label=label, markeredgecolor="black", markeredgewidth=0.5)

# pass threshold
ax.axhline(y=0.70, color="gray", linestyle="--", linewidth=1.5, alpha=0.8)
ax.text(len(alphabets) - 1, 0.71, "TM = 0.70", ha="right", va="bottom",
        fontsize=10, color="#2C3E50")

ax.set_xticks(x)
ax.set_xticklabels(alphabets, fontsize=12)
ax.set_xlabel("Reduced alphabet", fontsize=12)
ax.set_ylabel("Mean TM-score", fontsize=12)
ax.set_ylim(0.15, 1.0)
ax.legend(fontsize=10, loc="lower left")

plt.tight_layout()
plt.savefig("/storage/kiran-stuff/RAA/figures/fig1b_family_cliff.png", dpi=300, bbox_inches="tight")
plt.savefig("/storage/kiran-stuff/RAA/figures/fig1b_family_cliff.pdf", bbox_inches="tight")
print("Saved fig1b_family_cliff.png/pdf")
for fam, (label, _) in families.items():
    row = df[df["family"] == fam]
    if not row.empty:
        print(f"{label}: RAA10={row['RAA10_meanTM'].values[0]}  RAA9={row['RAA9_meanTM'].values[0]}")
