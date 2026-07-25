#!/usr/bin/env python3
"""Figure 3: EFS-guided RAA5 rescue of PF00069 kinase folds.

Blind RAA5 collapses (mean TM 0.399, 69% collapse); EFS-guided RAA5 (C20 core
protection) rescues folding (mean TM 0.791, 0% collapse). Distribution-based
overlaid histograms with the TM=0.70 pass threshold marked.

data paths:
  blind      = PF00069/blind-5/usalign_validation_report.csv       (per-structure)
  EFS-guided = PF00069/usalign_EFS_C_results.csv  (RAA5, C_cutoff=20 rows)
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# --- load per-structure TM distributions ---
blind = pd.read_csv("/storage/kiran-stuff/RAA/PF00069/blind-5/usalign_validation_report.csv")
blind_tm = blind["TM-score"].dropna().values

efs = pd.read_csv("/storage/kiran-stuff/RAA/PF00069/usalign_EFS_C_results.csv")
efs["raa"] = efs["variant_pdb"].str.extract(r"RAA(\d+)")
efs5 = efs[(efs["raa"] == "5") & (efs["C_cutoff"] == 20.0)]   # C20 core = manuscript definition
efs_tm = efs5["TM_score"].dropna().values

# --- summary stats (from the distributions themselves) ---
blind_mean, blind_collapse = blind_tm.mean(), (blind_tm < 0.5).mean() * 100
efs_mean, efs_collapse = efs_tm.mean(), (efs_tm < 0.5).mean() * 100

# --- plot: overlaid histograms ---
fig, ax = plt.subplots(figsize=(7, 6))

bins = np.linspace(0, 1, 41)
ax.hist(blind_tm, bins=bins, color="#E74C3C", alpha=0.6, edgecolor="black",
        linewidth=0.4, label=f"Blind RAA5  (n={len(blind_tm)})")
ax.hist(efs_tm, bins=bins, color="#27AE60", alpha=0.6, edgecolor="black",
        linewidth=0.4, label=f"EFS-guided RAA5  (n={len(efs_tm)})")

# pass threshold
ax.axvline(x=0.70, color="gray", linestyle="--", linewidth=1.5, alpha=0.8)
ax.text(0.70, ax.get_ylim()[1] * 0.97, " TM = 0.70\n (pass)", ha="left", va="top",
        fontsize=10, color="#2C3E50")

# mean markers + annotations
ax.axvline(x=blind_mean, color="#C0392B", linestyle="-", linewidth=1.2, alpha=0.9)
ax.axvline(x=efs_mean, color="#1E8449", linestyle="-", linewidth=1.2, alpha=0.9)

ax.annotate(f"blind μ = {blind_mean:.3f}\n{blind_collapse:.0f}% collapse",
            xy=(blind_mean, ax.get_ylim()[1] * 0.72), xytext=(0.05, ax.get_ylim()[1] * 0.72),
            fontsize=10, fontweight="bold", color="#C0392B", va="center",
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.85, edgecolor="#E74C3C"))
ax.annotate(f"EFS μ = {efs_mean:.3f}\n{efs_collapse:.0f}% collapse",
            xy=(efs_mean, ax.get_ylim()[1] * 0.55), xytext=(0.05, ax.get_ylim()[1] * 0.55),
            fontsize=10, fontweight="bold", color="#1E8449", va="center",
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.85, edgecolor="#27AE60"))

ax.set_xlabel("TM-score (fold recovery vs native)", fontsize=12)
ax.set_ylabel("Number of structures", fontsize=12)
ax.set_xlim(0, 1)
ax.legend(fontsize=11, loc="upper center")

plt.tight_layout()
plt.savefig("/storage/kiran-stuff/RAA/figures/fig3_efs_rescue.png", dpi=300, bbox_inches="tight")
plt.savefig("/storage/kiran-stuff/RAA/figures/fig3_efs_rescue.pdf", bbox_inches="tight")
print("Saved fig3_efs_rescue.png/pdf")
print(f"blind: n={len(blind_tm)} mean={blind_mean:.4f} collapse={blind_collapse:.1f}%")
print(f"EFS:   n={len(efs_tm)} mean={efs_mean:.4f} collapse={efs_collapse:.1f}%")
