#!/usr/bin/env python3
"""
generate_efs_vs_sasa_reductions.py

paired EFS-core vs SASA-core reduction experiment.

For each protein in the 5-rep set (PF00069 + PF00089):
  - extract sequence from family_core_scores.csv (sorted by res_index)
  - identify EFS core: positions with core_score >= 20
  - identify SASA core: positions with sasa < 25.0 angstroms^2
  - for each alphabet in {RAA5, RAA7_FIXED}:
      - generate 2 reduced sequences: EFS-core preserved, SASA-core preserved
  - write 4 FASTA records per protein

output FASTAs (all in FINAL_RESULTS/efs_vs_sasa/):
  PF00069_efs_vs_sasa_RAA5.fasta   - 10 records (5 EFS + 5 SASA)
  PF00069_efs_vs_sasa_RAA7.fasta
  PF00089_efs_vs_sasa_RAA5.fasta
  PF00089_efs_vs_sasa_RAA7.fasta
  combined.fasta                    - 40 records, ready for ESMFold

also writes a per-protein core-stats CSV showing the EFS-SASA core overlap.
"""

import sys
from pathlib import Path
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

# import canonical alphabets
sys.path.insert(0, "/storage/kiran-stuff/RAA/CANONICAL_CODE")
from CANONICAL_ALPHABETS import RAA5, RAA7_FIXED, reduce_at_positions

BASE = Path("/storage/kiran-stuff/RAA")
REPS_DIR = BASE / "FINAL_RESULTS" / "reps_5_per_family"
OUT_DIR = BASE / "FINAL_RESULTS" / "efs_vs_sasa"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# thresholds
EFS_THRESHOLD = 20   # core_score >= 20 (C20)
SASA_THRESHOLD = 25.0  # angstroms^2; standard burial cutoff

# alphabets to test
ALPHABETS = {
    "RAA5": RAA5,
    "RAA7": RAA7_FIXED,
}


def get_reps(family):
    """load 5 representative protein IDs for a family from the FASTA."""
    fa = REPS_DIR / f"{family}_reps.fasta"
    return [r.id for r in SeqIO.parse(fa, "fasta")]


def extract_protein_data(family, protein_id):
    """
    extract sequence + core info from family_core_scores.csv.
    returns: (native_sequence, efs_positions_0idx, sasa_positions_0idx, df)
    positions are 0-indexed into native_sequence.
    """
    csv_path = BASE / family / "family_core_scores.csv"
    df = pd.read_csv(csv_path)
    sub = df[df["protein_id"] == protein_id].sort_values("res_index").reset_index(drop=True)
    if len(sub) == 0:
        return None, None, None, None

    native_seq = "".join(sub["aa"].tolist())

    # 0-indexed positions
    efs_positions = set(sub.index[sub["core_score"] >= EFS_THRESHOLD].tolist())
    sasa_positions = set(sub.index[sub["sasa"] < SASA_THRESHOLD].tolist())

    return native_seq, efs_positions, sasa_positions, sub


def reduce_non_core(native_seq, core_positions, alphabet):
    """
    preserve native AA at core positions, reduce all OTHER positions
    using alphabet. returns reduced sequence string.
    """
    non_core = set(range(len(native_seq))) - core_positions
    return reduce_at_positions(native_seq, non_core, alphabet)


def main():
    families = ["PF00069", "PF00089"]
    core_stats = []
    all_records = []  # the mega combined FASTA

    for family in families:
        reps = get_reps(family)
        print(f"\n=== {family}: {len(reps)} reps ===")

        per_alpha_records = {alpha: [] for alpha in ALPHABETS}

        for pid in reps:
            seq, efs_pos, sasa_pos, _ = extract_protein_data(family, pid)
            if seq is None:
                print(f"  ! no data for {pid}")
                continue

            L = len(seq)
            n_efs = len(efs_pos)
            n_sasa = len(sasa_pos)
            overlap = len(efs_pos & sasa_pos)
            efs_unique = len(efs_pos - sasa_pos)
            sasa_unique = len(sasa_pos - efs_pos)
            jacc = overlap / len(efs_pos | sasa_pos) if (efs_pos | sasa_pos) else 0.0

            core_stats.append({
                "family": family, "protein_id": pid, "length": L,
                "n_efs_core": n_efs, "n_sasa_core": n_sasa,
                "overlap": overlap,
                "efs_unique": efs_unique, "sasa_unique": sasa_unique,
                "jaccard": round(jacc, 3),
                "pct_efs_core": round(100 * n_efs / L, 1),
                "pct_sasa_core": round(100 * n_sasa / L, 1),
            })

            for alpha_name, alphabet in ALPHABETS.items():
                # EFS-core preserved + non-core reduced
                efs_reduced = reduce_non_core(seq, efs_pos, alphabet)
                rec_efs = SeqRecord(
                    Seq(efs_reduced),
                    id=f"{family}_{pid}_EFScore_{alpha_name}",
                    description=f"EFS core ({n_efs}/{L}) preserved, non-core -> {alpha_name}",
                )
                # SASA-core preserved + non-core reduced
                sasa_reduced = reduce_non_core(seq, sasa_pos, alphabet)
                rec_sasa = SeqRecord(
                    Seq(sasa_reduced),
                    id=f"{family}_{pid}_SASAcore_{alpha_name}",
                    description=f"SASA core ({n_sasa}/{L}) preserved, non-core -> {alpha_name}",
                )
                per_alpha_records[alpha_name].append(rec_efs)
                per_alpha_records[alpha_name].append(rec_sasa)
                all_records.append(rec_efs)
                all_records.append(rec_sasa)

            print(f"  {pid}: L={L}, EFS={n_efs} ({n_efs/L*100:.1f}%), "
                  f"SASA={n_sasa} ({n_sasa/L*100:.1f}%), jacc={jacc:.3f}")

        # write per-alphabet FASTAs
        for alpha_name, recs in per_alpha_records.items():
            out = OUT_DIR / f"{family}_efs_vs_sasa_{alpha_name}.fasta"
            SeqIO.write(recs, out, "fasta")
            print(f"  wrote {len(recs)} -> {out.name}")

    # write combined FASTA for one-shot ESMFold batch
    SeqIO.write(all_records, OUT_DIR / "combined.fasta", "fasta")
    print(f"\ncombined: {len(all_records)} records -> combined.fasta")

    # write core stats CSV
    stats_df = pd.DataFrame(core_stats)
    stats_df.to_csv(OUT_DIR / "core_stats.csv", index=False)
    print(f"core stats -> core_stats.csv")
    print("\nCore size summary:")
    print(stats_df[["family", "protein_id", "length", "n_efs_core",
                    "n_sasa_core", "efs_unique", "sasa_unique", "jaccard"]].to_string(index=False))


if __name__ == "__main__":
    main()
