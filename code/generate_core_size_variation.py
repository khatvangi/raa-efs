#!/usr/bin/env python3
"""
generate_core_size_variation.py - Core size variation experiment

Tests whether EFS matters at different core sizes:
- C10, C20, C30, C40 (core_score thresholds)
- For each cutoff, generates:
  1. EFS_PRESERVED - core native, non-core reduced (RAA10)
  2. ANTI_EFS_INVERTED - core chemistry inverted, non-core native

Expected results:
- EFS_PRESERVED should work at all cutoffs (90%+ pass)
- ANTI_EFS_INVERTED should fail at all cutoffs (proves core matters)
- Smaller core (C40) + inversion = more catastrophic failure

This proves: EFS core matters, regardless of size. Non-obvious.

Created: 2026-01-11
"""

import pandas as pd
import numpy as np
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys

# import canonical alphabets
sys.path.insert(0, str(Path(__file__).parent))
from CANONICAL_ALPHABETS import RAA10, CHEMISTRY_INVERSION, reduce_at_positions

def invert_core_chemistry(sequence: str, core_positions: set) -> str:
    """invert chemistry at core positions (hydrophobic -> charged)"""
    seq_list = list(sequence)
    for i in core_positions:
        if i < len(seq_list):
            aa = seq_list[i]
            seq_list[i] = CHEMISTRY_INVERSION.get(aa, aa)
    return "".join(seq_list)


def load_per_protein_core_scores(base_dir: Path, family: str, protein_id: str) -> pd.DataFrame:
    """load core scores for a specific protein"""
    path = base_dir / family / "family_EFS" / "per_protein" / f"{protein_id}_core_scores.csv"
    if path.exists():
        return pd.read_csv(path)
    raise FileNotFoundError(f"Core scores not found: {path}")


def get_core_positions_by_score(df: pd.DataFrame, threshold: float) -> set:
    """get positions with core_score >= threshold (0-indexed)"""
    # res_index is 1-based, convert to 0-based
    core_mask = df["core_score"] >= threshold
    positions = set(df.loc[core_mask, "res_index"].values - 1)  # convert to 0-indexed
    return positions


def main():
    base_dir = Path("/storage/kiran-stuff/RAA")
    output_dir = base_dir / "core_size_variation"
    output_dir.mkdir(parents=True, exist_ok=True)

    # proteins from chemistry inversion experiment
    proteins = {
        "PF00069": ["P0CP66", "Q6NW76", "Q9N0X0", "O75716", "P36887",
                    "Q0J4I1", "A2QRF6", "Q8N165", "O76484", "Q4U925"],
        "PF01145": ["P9WPR9", "Q2HJ97", "P35232", "O51221", "Q21190",
                    "A3QK16", "Q5ZMN3", "O94550", "O35129", "P0AA56"],
        "PF07686": ["P33681", "Q9H6B4", "P06907", "Q6AYD4", "Q2YHT5",
                    "Q149L7", "Q8TDQ1", "Q07699", "A5A6L6", "P40199"]
    }

    # core_score ranges 0-32, use meaningful thresholds
    # C8: ~70% core (large), C12: ~30% (medium), C20: ~20% (small), C28: ~15% (tiny)
    cutoffs = [8, 12, 20, 28]

    print("=" * 70)
    print("CORE SIZE VARIATION EXPERIMENT")
    print("=" * 70)
    print(f"Proteins: {sum(len(v) for v in proteins.values())}")
    print(f"Cutoffs: C{cutoffs}")
    print(f"Conditions: EFS_PRESERVED, ANTI_EFS_INVERTED")
    print("=" * 70)

    all_records = []
    stats = []

    for family, protein_list in proteins.items():
        print(f"\n--- {family} ---")

        for protein_id in protein_list:
            try:
                # load core scores
                df = load_per_protein_core_scores(base_dir, family, protein_id)
                df = df.sort_values("res_index")

                # get sequence
                sequence = "".join(df["aa"].values)

                for cutoff in cutoffs:
                    # get core positions at this cutoff
                    core_positions = get_core_positions_by_score(df, cutoff)
                    non_core_positions = set(range(len(sequence))) - core_positions

                    n_core = len(core_positions)
                    core_pct = n_core / len(sequence) * 100

                    # 1. EFS_PRESERVED - core native, non-core reduced
                    efs_seq = reduce_at_positions(sequence, non_core_positions, RAA10)
                    header1 = f"{protein_id}_C{cutoff}_EFS_PRESERVED"
                    all_records.append(SeqRecord(Seq(efs_seq), id=header1, description=""))

                    # 2. ANTI_EFS_INVERTED - core inverted, non-core native
                    inv_seq = invert_core_chemistry(sequence, core_positions)
                    header2 = f"{protein_id}_C{cutoff}_ANTI_EFS_INVERTED"
                    all_records.append(SeqRecord(Seq(inv_seq), id=header2, description=""))

                    stats.append({
                        "family": family,
                        "protein_id": protein_id,
                        "cutoff": cutoff,
                        "length": len(sequence),
                        "n_core": n_core,
                        "core_pct": core_pct
                    })

                print(f"  {protein_id}: {len(sequence)} AA, cores at C10-40 generated")

            except Exception as e:
                print(f"  {protein_id}: ERROR - {e}")

    # also add NATIVE for each protein (for reference)
    for family, protein_list in proteins.items():
        for protein_id in protein_list:
            try:
                df = load_per_protein_core_scores(base_dir, family, protein_id)
                df = df.sort_values("res_index")
                sequence = "".join(df["aa"].values)
                header = f"{protein_id}_NATIVE"
                all_records.append(SeqRecord(Seq(sequence), id=header, description=""))
            except:
                pass

    # write FASTA
    fasta_path = output_dir / "core_size_variation.fasta"
    with open(fasta_path, "w") as f:
        SeqIO.write(all_records, f, "fasta")

    # write stats
    stats_df = pd.DataFrame(stats)
    stats_path = output_dir / "experiment_stats.csv"
    stats_df.to_csv(stats_path, index=False)

    # summary by cutoff
    print("\n" + "=" * 70)
    print("CORE SIZE BY CUTOFF (mean across all proteins)")
    print("=" * 70)
    for cutoff in cutoffs:
        sub = stats_df[stats_df["cutoff"] == cutoff]
        print(f"C{cutoff}: {sub['core_pct'].mean():.1f}% core ({sub['n_core'].mean():.0f} positions)")

    print(f"\nTotal sequences: {len(all_records)}")
    print(f"  - 30 proteins × 4 cutoffs × 2 conditions = 240")
    print(f"  - + 30 NATIVE references")
    print(f"  - = {len(all_records)} total")
    print(f"\nOutput: {fasta_path}")
    print(f"Stats: {stats_path}")
    print("\nNEXT STEP:")
    print(f"  trill --fold esmfold {fasta_path}")
    print("=" * 70)


if __name__ == "__main__":
    main()
