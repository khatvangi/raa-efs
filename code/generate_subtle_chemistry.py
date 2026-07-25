#!/usr/bin/env python3
"""
generate_subtle_chemistry.py - Subtle chemistry variation experiment

Tests 6 different subtle chemistry changes at core positions:
A: Size reduction (keep hydrophobicity) - L→A, I→A, V→A, F→A, W→A, M→A
B: Polarity shift (similar size) - L→N, I→Q, V→T, F→Y, A→S
C: Aromatic→Aliphatic - F→L, Y→L, W→L
D: Charge neutralization - K→Q, R→Q, D→N, E→Q
E: H-bond removal - S→A, T→V, N→L, Q→L
F: Flexibility change - P→A, G→A

This creates a hierarchy of which chemistry aspects matter most for core stability.

Created: 2026-01-12
"""

import pandas as pd
import numpy as np
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys

# subtle chemistry mappings - only specified AAs change, others stay same
SUBTLE_A = {  # size reduction (keep hydrophobic)
    'L': 'A', 'I': 'A', 'V': 'A', 'F': 'A', 'W': 'A', 'M': 'A'
}

SUBTLE_B = {  # polarity shift (similar size)
    'L': 'N', 'I': 'Q', 'V': 'T', 'F': 'Y', 'A': 'S'
}

SUBTLE_C = {  # aromatic → aliphatic
    'F': 'L', 'Y': 'L', 'W': 'L'
}

SUBTLE_D = {  # charge neutralization
    'K': 'Q', 'R': 'Q', 'D': 'N', 'E': 'Q'
}

SUBTLE_E = {  # H-bond removal
    'S': 'A', 'T': 'V', 'N': 'L', 'Q': 'L'
}

SUBTLE_F = {  # flexibility change
    'P': 'A', 'G': 'A'
}

SUBTLE_OPTIONS = {
    'A_SIZE': SUBTLE_A,
    'B_POLAR': SUBTLE_B,
    'C_AROM': SUBTLE_C,
    'D_CHARGE': SUBTLE_D,
    'E_HBOND': SUBTLE_E,
    'F_FLEX': SUBTLE_F
}


def apply_subtle_change(sequence: str, core_positions: set, mapping: dict) -> str:
    """apply subtle chemistry change at core positions only"""
    seq_list = list(sequence)
    for i in core_positions:
        if i < len(seq_list):
            aa = seq_list[i]
            # only change if AA is in the mapping, otherwise keep original
            seq_list[i] = mapping.get(aa, aa)
    return "".join(seq_list)


def load_per_protein_core_scores(base_dir: Path, family: str, protein_id: str) -> pd.DataFrame:
    """load core scores for a specific protein"""
    path = base_dir / family / "family_EFS" / "per_protein" / f"{protein_id}_core_scores.csv"
    if path.exists():
        return pd.read_csv(path)
    raise FileNotFoundError(f"Core scores not found: {path}")


def get_core_positions_by_score(df: pd.DataFrame, threshold: float) -> set:
    """get positions with core_score >= threshold (0-indexed)"""
    core_mask = df["core_score"] >= threshold
    positions = set(df.loc[core_mask, "res_index"].values - 1)  # convert to 0-indexed
    return positions


def main():
    base_dir = Path("/storage/kiran-stuff/RAA")
    output_dir = base_dir / "subtle_chemistry"
    output_dir.mkdir(parents=True, exist_ok=True)

    # same proteins as before
    proteins = {
        "PF00069": ["P0CP66", "Q6NW76", "Q9N0X0", "O75716", "P36887",
                    "Q0J4I1", "A2QRF6", "Q8N165", "O76484", "Q4U925"],
        "PF01145": ["P9WPR9", "Q2HJ97", "P35232", "O51221", "Q21190",
                    "A3QK16", "Q5ZMN3", "O94550", "O35129", "P0AA56"],
        "PF07686": ["P33681", "Q9H6B4", "P06907", "Q6AYD4", "Q2YHT5",
                    "Q149L7", "Q8TDQ1", "Q07699", "A5A6L6", "P40199"]
    }

    cutoffs = [8, 12, 20, 28]

    print("=" * 70)
    print("SUBTLE CHEMISTRY EXPERIMENT")
    print("=" * 70)
    print(f"Proteins: {sum(len(v) for v in proteins.values())}")
    print(f"Cutoffs: C{cutoffs}")
    print(f"Options: {list(SUBTLE_OPTIONS.keys())}")
    print("=" * 70)

    all_records = []
    stats = []

    for family, protein_list in proteins.items():
        print(f"\n--- {family} ---")

        for protein_id in protein_list:
            try:
                df = load_per_protein_core_scores(base_dir, family, protein_id)
                df = df.sort_values("res_index")
                sequence = "".join(df["aa"].values)

                for cutoff in cutoffs:
                    core_positions = get_core_positions_by_score(df, cutoff)
                    n_core = len(core_positions)

                    for option_name, mapping in SUBTLE_OPTIONS.items():
                        # apply subtle change at core positions
                        modified_seq = apply_subtle_change(sequence, core_positions, mapping)

                        # count how many positions actually changed
                        n_changed = sum(1 for i in core_positions if i < len(sequence) and sequence[i] in mapping)

                        header = f"{protein_id}_C{cutoff}_{option_name}"
                        all_records.append(SeqRecord(Seq(modified_seq), id=header, description=""))

                        stats.append({
                            "family": family,
                            "protein_id": protein_id,
                            "cutoff": cutoff,
                            "option": option_name,
                            "length": len(sequence),
                            "n_core": n_core,
                            "n_changed": n_changed,
                            "pct_changed": n_changed / len(sequence) * 100
                        })

                print(f"  {protein_id}: {len(sequence)} AA, 4 cutoffs × 6 options generated")

            except Exception as e:
                print(f"  {protein_id}: ERROR - {e}")

    # add NATIVE for reference
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
    fasta_path = output_dir / "subtle_chemistry.fasta"
    with open(fasta_path, "w") as f:
        SeqIO.write(all_records, f, "fasta")

    # write stats
    stats_df = pd.DataFrame(stats)
    stats_path = output_dir / "experiment_stats.csv"
    stats_df.to_csv(stats_path, index=False)

    # summary
    print("\n" + "=" * 70)
    print("SUMMARY BY OPTION (mean positions changed)")
    print("=" * 70)
    for option in SUBTLE_OPTIONS.keys():
        sub = stats_df[stats_df["option"] == option]
        print(f"{option}: {sub['n_changed'].mean():.1f} positions ({sub['pct_changed'].mean():.1f}%)")

    print(f"\nTotal sequences: {len(all_records)}")
    print(f"  - 30 proteins × 4 cutoffs × 6 options = 720")
    print(f"  - + 30 NATIVE references")
    print(f"  - = {len(all_records)} total")
    print(f"\nOutput: {fasta_path}")
    print("\nNEXT STEP:")
    print(f"  trill --fold esmfold {fasta_path}")
    print("=" * 70)


if __name__ == "__main__":
    main()
