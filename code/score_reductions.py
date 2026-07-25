#!/usr/bin/env python3
"""
score_reductions.py - USalign each predicted reduction PDB against its native

Parses predicted PDBs named like PF00069_O70126_RAA8.pdb, finds the native
at /storage/kiran-stuff/RAA/<family>/Original/<protein_id>.pdb, runs USalign,
and compiles results into a CSV + summary table.

Output: FINAL_RESULTS/reps_5_per_family/raa89_validation_results.csv
"""
import subprocess
import re
from pathlib import Path
import pandas as pd
import sys

BASE = Path("/storage/kiran-stuff/RAA")
PRED_DIR = BASE / "FINAL_RESULTS/reps_5_per_family/folding_outputs/full_run"
OUT_DIR = BASE / "FINAL_RESULTS/reps_5_per_family"
USALIGN = "/storage/kiran-stuff/US-align/USalign"


def parse_tm_score(output: str) -> dict:
    """parse USalign stdout. returns dict with tm1, tm2, rmsd, seq_id."""
    result = {"tm1": None, "tm2": None, "rmsd": None, "seq_id": None, "aligned_len": None}
    for line in output.split("\n"):
        # TM-score= 0.92538 (normalized by length of Structure_1: L=345, d0=6.77)
        m = re.search(r"TM-score=\s*([\d.]+)\s*\(normalized by length of Structure_1", line)
        if m:
            result["tm1"] = float(m.group(1))
        m = re.search(r"TM-score=\s*([\d.]+)\s*\(normalized by length of Structure_2", line)
        if m:
            result["tm2"] = float(m.group(1))
        # Aligned length= 340, RMSD=   2.35, Seq_ID=n_identical/n_aligned= 0.785
        m = re.search(r"Aligned length=\s*(\d+),\s*RMSD=\s*([\d.]+),\s*Seq_ID=n_identical/n_aligned=\s*([\d.]+)", line)
        if m:
            result["aligned_len"] = int(m.group(1))
            result["rmsd"] = float(m.group(2))
            result["seq_id"] = float(m.group(3))
    return result


def run_usalign(native_pdb: Path, predicted_pdb: Path) -> dict:
    """run USalign, return parsed metrics."""
    try:
        result = subprocess.run(
            [USALIGN, str(predicted_pdb), str(native_pdb)],
            capture_output=True, text=True, timeout=120,
        )
        return parse_tm_score(result.stdout)
    except Exception as exc:
        return {"error": str(exc)}


def main():
    if not PRED_DIR.exists():
        print(f"error: {PRED_DIR} missing")
        sys.exit(1)

    pdbs = sorted(PRED_DIR.glob("*.pdb"))
    print(f"found {len(pdbs)} predicted PDBs")

    records = []
    for i, pdb in enumerate(pdbs, 1):
        # parse name: {family}_{protein}_RAA{N}.pdb
        stem = pdb.stem
        m = re.match(r"(PF\d+)_([A-Z0-9]+)_RAA(\d+)", stem)
        if not m:
            print(f"  ! cannot parse name: {stem}")
            continue
        family, protein_id, raa = m.group(1), m.group(2), int(m.group(3))

        native = BASE / family / "Original" / f"{protein_id}.pdb"
        if not native.exists():
            print(f"  ! native missing for {protein_id}: {native}")
            continue

        metrics = run_usalign(native, pdb)
        record = {
            "family": family,
            "protein_id": protein_id,
            "raa": raa,
            "predicted_pdb": str(pdb),
            "native_pdb": str(native),
            **metrics,
        }
        records.append(record)
        if i % 10 == 0 or i == len(pdbs):
            tm = metrics.get("tm1", "?")
            print(f"  [{i}/{len(pdbs)}] {family}/{protein_id}/RAA{raa}: TM={tm}")

    df = pd.DataFrame(records)
    out_csv = OUT_DIR / "raa89_validation_results.csv"
    df.to_csv(out_csv, index=False)
    print(f"\nwrote {len(df)} rows -> {out_csv}")

    # summary: mean TM per (family, raa)
    print("\nMean TM-score (normalized by Chain_1 = native) per family x RAA:")
    pivot = df.pivot_table(values="tm1", index="family", columns="raa", aggfunc="mean")
    print(pivot.round(3).to_string())

    # overall summary
    print("\nOverall mean TM per RAA:")
    print(df.groupby("raa")["tm1"].agg(["mean", "std", "count"]).round(3).to_string())

    summary_csv = OUT_DIR / "raa89_validation_summary.csv"
    summary = df.groupby(["family", "raa"]).agg(
        n=("tm1", "count"),
        mean_tm=("tm1", "mean"),
        std_tm=("tm1", "std"),
        min_tm=("tm1", "min"),
        max_tm=("tm1", "max"),
    ).round(3)
    summary.to_csv(summary_csv)
    print(f"\nsummary -> {summary_csv}")


if __name__ == "__main__":
    main()
