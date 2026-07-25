#!/usr/bin/env python3
"""
score_efs_vs_sasa.py

USalign each predicted PDB against its native. Produces paired
(EFS, SASA) TM-score comparison at each alphabet level and runs
a paired Wilcoxon test.

predicted PDBs are named like:
  PF00069_O70126_EFScore_RAA5.pdb
  PF00069_O70126_SASAcore_RAA5.pdb
"""

import re
import subprocess
import sys
from pathlib import Path
import pandas as pd

BASE = Path("/storage/kiran-stuff/RAA")
PRED_DIR = BASE / "FINAL_RESULTS/efs_vs_sasa/folding_outputs"
OUT_DIR = BASE / "FINAL_RESULTS/efs_vs_sasa"
USALIGN = "/storage/kiran-stuff/US-align/USalign"

NAME_RE = re.compile(r"(PF\d+)_([A-Z0-9]+)_(EFScore|SASAcore)_(RAA\d+)")


def parse_tm(stdout):
    """parse USalign output -> dict with tm_native, tm_pred, rmsd, seq_id."""
    res = {"tm_native": None, "tm_pred": None, "rmsd": None, "seq_id": None}
    for line in stdout.split("\n"):
        # Structure_1 = predicted (passed first), Structure_2 = native (second)
        m = re.search(r"TM-score=\s*([\d.]+)\s*\(normalized by length of Structure_1", line)
        if m:
            res["tm_pred"] = float(m.group(1))
        m = re.search(r"TM-score=\s*([\d.]+)\s*\(normalized by length of Structure_2", line)
        if m:
            res["tm_native"] = float(m.group(1))
        m = re.search(r"Aligned length=\s*\d+,\s*RMSD=\s*([\d.]+),\s*Seq_ID=n_identical/n_aligned=\s*([\d.]+)", line)
        if m:
            res["rmsd"] = float(m.group(1))
            res["seq_id"] = float(m.group(2))
    return res


def run_usalign(predicted, native):
    try:
        proc = subprocess.run(
            [USALIGN, str(predicted), str(native)],
            capture_output=True, text=True, timeout=120,
        )
        return parse_tm(proc.stdout)
    except Exception as e:
        return {"error": str(e)}


def main():
    pdbs = sorted(PRED_DIR.glob("*.pdb"))
    print(f"found {len(pdbs)} predicted PDBs")
    rows = []
    for i, pdb in enumerate(pdbs, 1):
        m = NAME_RE.match(pdb.stem)
        if not m:
            print(f"  ! cannot parse: {pdb.stem}")
            continue
        family, pid, core_type, raa = m.groups()
        native = BASE / family / "Original" / f"{pid}.pdb"
        if not native.exists():
            print(f"  ! missing native: {native}")
            continue
        metrics = run_usalign(pdb, native)
        row = {
            "family": family, "protein_id": pid,
            "core_type": "EFS" if core_type == "EFScore" else "SASA",
            "raa": raa,
            **metrics,
        }
        rows.append(row)
        if i % 5 == 0 or i == len(pdbs):
            print(f"  [{i}/{len(pdbs)}] {family}/{pid}/{core_type}/{raa}: TM_native={metrics.get('tm_native')}")
    df = pd.DataFrame(rows)
    df.to_csv(OUT_DIR / "efs_vs_sasa_results.csv", index=False)
    print(f"\nresults -> efs_vs_sasa_results.csv ({len(df)} rows)")

    # paired comparison: pivot to wide so EFS and SASA are side by side
    print("\n" + "=" * 70)
    print("PAIRED EFS vs SASA COMPARISON (TM-score normalized by native length)")
    print("=" * 70)

    pivot = df.pivot_table(
        values="tm_native",
        index=["family", "protein_id", "raa"],
        columns="core_type",
        aggfunc="first",
    ).reset_index()
    pivot["delta_efs_minus_sasa"] = pivot["EFS"] - pivot["SASA"]
    pivot.to_csv(OUT_DIR / "efs_vs_sasa_paired.csv", index=False)
    print(f"\npaired data -> efs_vs_sasa_paired.csv\n")
    print(pivot.to_string(index=False))

    # summary per (family, RAA)
    print("\n" + "=" * 70)
    print("SUMMARY BY (family, RAA)")
    print("=" * 70)
    summary = (
        pivot.groupby(["family", "raa"])
        .agg(
            n=("delta_efs_minus_sasa", "count"),
            mean_EFS=("EFS", "mean"),
            mean_SASA=("SASA", "mean"),
            mean_delta=("delta_efs_minus_sasa", "mean"),
            efs_wins=("delta_efs_minus_sasa", lambda x: int((x > 0).sum())),
            ties=("delta_efs_minus_sasa", lambda x: int((x == 0).sum())),
        )
        .round(3)
        .reset_index()
    )
    summary.to_csv(OUT_DIR / "efs_vs_sasa_summary.csv", index=False)
    print(summary.to_string(index=False))

    # overall paired test (try scipy if available)
    print("\n" + "=" * 70)
    print("OVERALL PAIRED COMPARISON")
    print("=" * 70)
    for raa in sorted(pivot["raa"].unique()):
        sub = pivot[pivot["raa"] == raa].dropna(subset=["EFS", "SASA"])
        delta = sub["EFS"].values - sub["SASA"].values
        n_wins = int((delta > 0).sum())
        n_ties = int((delta == 0).sum())
        n_losses = int((delta < 0).sum())
        try:
            from scipy.stats import wilcoxon
            stat, p = wilcoxon(sub["EFS"], sub["SASA"], alternative="greater")
            p_str = f", Wilcoxon (EFS>SASA) p={p:.4f}"
        except Exception:
            p_str = ""
        print(f"  {raa}: EFS wins {n_wins}/{len(sub)} (ties={n_ties}, losses={n_losses}), "
              f"mean ΔTM = {delta.mean():+.3f}{p_str}")


if __name__ == "__main__":
    main()
