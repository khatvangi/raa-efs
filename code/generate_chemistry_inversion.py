#!/usr/bin/env python3
"""
generate_chemistry_inversion.py - Chemistry inversion experiment

Generates 5 conditions per protein for hypothesis testing:
1. NATIVE - original sequence (positive control)
2. EFS_PRESERVED - core native, non-core reduced (RAA10)
3. ANTI_EFS_REDUCED - core reduced (RAA10), non-core native
4. ANTI_EFS_INVERTED - core chemistry-inverted, non-core native (KEY TEST)
5. BLIND_REDUCED - all positions reduced (RAA10)

If hypothesis is correct:
- NATIVE, EFS_PRESERVED should fold well (TM > 0.85)
- BLIND_REDUCED should partially work (TM 0.6-0.8)
- ANTI_EFS_REDUCED should struggle (TM 0.5-0.7)
- ANTI_EFS_INVERTED should collapse (TM < 0.5) <-- KEY RESULT

Usage:
    python generate_chemistry_inversion.py \
        --families PF00089 PF00069 PF07686 \
        --proteins_per_family 20 20 10 \
        --core_threshold 0.3 \
        --output_dir /storage/kiran-stuff/RAA/chemistry_inversion/

Created: 2026-01-10
"""

import argparse
import pandas as pd
import numpy as np
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys

# import canonical alphabets
sys.path.insert(0, str(Path(__file__).parent))
from CANONICAL_ALPHABETS import (
    RAA10, CHEMISTRY_INVERSION,
    reduce_sequence, reduce_at_positions, invert_core_chemistry
)


def load_core_scores(family_dir: Path) -> pd.DataFrame:
    """load family core scores CSV"""
    # try different possible locations
    possible_paths = [
        family_dir / "family_EFS" / "family_core_scores_with_familyEFS.csv",
        family_dir / "family_EFS" / "family_core_scores.csv",
        family_dir / "family_core_scores.csv",
    ]

    for path in possible_paths:
        if path.exists():
            return pd.read_csv(path)

    raise FileNotFoundError(f"No core scores found in {family_dir}")


def load_sequences(family_dir: Path) -> dict:
    """load original sequences from FASTA or PDB directory"""
    sequences = {}

    # try FASTA first
    fasta_files = list(family_dir.glob("*.fasta")) + list(family_dir.glob("*.fa"))
    for fasta in fasta_files:
        for record in SeqIO.parse(fasta, "fasta"):
            protein_id = record.id.split("_")[0].split("/")[0]
            sequences[protein_id] = str(record.seq)

    if sequences:
        return sequences

    # try to extract from core scores CSV (if has sequence column)
    try:
        df = load_core_scores(family_dir)
        for protein_id, group in df.groupby("protein_id"):
            # reconstruct sequence from aa column
            group = group.sort_values("res_index")
            seq = "".join(group["aa"].values)
            sequences[protein_id] = seq
    except:
        pass

    return sequences


def get_core_positions(df: pd.DataFrame, protein_id: str, threshold_percentile: float = 0.3) -> set:
    """
    get core positions for a protein based on core_score percentile

    args:
        df: core scores dataframe
        protein_id: protein to get cores for
        threshold_percentile: top X% are considered core (0.3 = top 30%)

    returns:
        set of 0-indexed core positions
    """
    protein_df = df[df["protein_id"] == protein_id].copy()
    if protein_df.empty:
        return set()

    # use core_score column
    if "core_score" not in protein_df.columns:
        raise ValueError(f"core_score column not found for {protein_id}")

    # get threshold value for top X%
    threshold = protein_df["core_score"].quantile(1 - threshold_percentile)

    # get positions above threshold
    core_mask = protein_df["core_score"] >= threshold
    core_positions = set(protein_df.loc[core_mask, "res_index"].values)

    return core_positions


def generate_conditions(sequence: str, core_positions: set, protein_id: str) -> list:
    """
    generate 5 experimental conditions for a protein

    returns list of (header, sequence) tuples
    """
    non_core_positions = set(range(len(sequence))) - core_positions

    conditions = []

    # 1. NATIVE - original sequence
    conditions.append((
        f"{protein_id}_NATIVE",
        sequence
    ))

    # 2. EFS_PRESERVED - core native, non-core reduced (RAA10)
    efs_seq = reduce_at_positions(sequence, non_core_positions, RAA10)
    conditions.append((
        f"{protein_id}_EFS_PRESERVED",
        efs_seq
    ))

    # 3. ANTI_EFS_REDUCED - core reduced (RAA10), non-core native
    anti_efs_reduced = reduce_at_positions(sequence, core_positions, RAA10)
    conditions.append((
        f"{protein_id}_ANTI_EFS_REDUCED",
        anti_efs_reduced
    ))

    # 4. ANTI_EFS_INVERTED - core chemistry-inverted, non-core native (KEY TEST)
    anti_efs_inverted = invert_core_chemistry(sequence, core_positions)
    conditions.append((
        f"{protein_id}_ANTI_EFS_INVERTED",
        anti_efs_inverted
    ))

    # 5. BLIND_REDUCED - all positions reduced (RAA10)
    blind_seq = reduce_sequence(sequence, RAA10)
    conditions.append((
        f"{protein_id}_BLIND_REDUCED",
        blind_seq
    ))

    return conditions


def select_proteins(df: pd.DataFrame, n_proteins: int, seed: int = 42) -> list:
    """select n proteins with good core coverage"""
    np.random.seed(seed)

    # get unique proteins
    proteins = df["protein_id"].unique()

    # calculate core stats per protein
    stats = []
    for pid in proteins:
        pdata = df[df["protein_id"] == pid]
        core_pct = (pdata["core_score"] > pdata["core_score"].median()).mean()
        length = len(pdata)
        stats.append({"protein_id": pid, "core_pct": core_pct, "length": length})

    stats_df = pd.DataFrame(stats)

    # filter for reasonable lengths (100-500 AA)
    stats_df = stats_df[(stats_df["length"] >= 100) & (stats_df["length"] <= 500)]

    # select diverse set
    if len(stats_df) <= n_proteins:
        selected = stats_df["protein_id"].tolist()
    else:
        selected = stats_df.sample(n=n_proteins, random_state=seed)["protein_id"].tolist()

    return selected


def main():
    parser = argparse.ArgumentParser(description="Generate chemistry inversion experiment")
    parser.add_argument("--families", nargs="+", required=True,
                        help="Family directories (e.g., PF00089 PF00069)")
    parser.add_argument("--proteins_per_family", nargs="+", type=int, required=True,
                        help="Number of proteins per family")
    parser.add_argument("--base_dir", default="/storage/kiran-stuff/RAA",
                        help="Base RAA directory")
    parser.add_argument("--core_threshold", type=float, default=0.3,
                        help="Core threshold (top X%%, default 0.3 = top 30%%)")
    parser.add_argument("--output_dir", required=True,
                        help="Output directory for FASTA files")
    parser.add_argument("--seed", type=int, default=42,
                        help="Random seed for protein selection")
    args = parser.parse_args()

    if len(args.families) != len(args.proteins_per_family):
        print("ERROR: families and proteins_per_family must have same length")
        sys.exit(1)

    base_dir = Path(args.base_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("CHEMISTRY INVERSION EXPERIMENT - SEQUENCE GENERATION")
    print("=" * 70)
    print(f"Families: {args.families}")
    print(f"Proteins per family: {args.proteins_per_family}")
    print(f"Core threshold: top {args.core_threshold*100:.0f}%")
    print(f"Output: {output_dir}")
    print("=" * 70)

    all_records = []
    stats = []

    for family, n_proteins in zip(args.families, args.proteins_per_family):
        print(f"\n--- Processing {family} ({n_proteins} proteins) ---")

        family_dir = base_dir / family
        if not family_dir.exists():
            print(f"  WARNING: {family_dir} not found, skipping")
            continue

        # load core scores
        try:
            core_df = load_core_scores(family_dir)
            print(f"  Loaded core scores: {len(core_df)} residues")
        except FileNotFoundError as e:
            print(f"  ERROR: {e}")
            continue

        # select proteins
        selected = select_proteins(core_df, n_proteins, args.seed)
        print(f"  Selected {len(selected)} proteins")

        # generate conditions for each protein
        for protein_id in selected:
            # get sequence from core scores
            pdata = core_df[core_df["protein_id"] == protein_id].sort_values("res_index")
            if pdata.empty:
                continue

            sequence = "".join(pdata["aa"].values)

            # get core positions
            core_positions = get_core_positions(core_df, protein_id, args.core_threshold)

            if len(core_positions) < 10:
                print(f"    {protein_id}: only {len(core_positions)} core positions, skipping")
                continue

            # generate 5 conditions
            conditions = generate_conditions(sequence, core_positions, protein_id)

            for header, seq in conditions:
                record = SeqRecord(Seq(seq), id=header, description="")
                all_records.append(record)

            stats.append({
                "family": family,
                "protein_id": protein_id,
                "length": len(sequence),
                "n_core": len(core_positions),
                "core_pct": len(core_positions) / len(sequence) * 100
            })

            print(f"    {protein_id}: {len(sequence)} AA, {len(core_positions)} core ({len(core_positions)/len(sequence)*100:.1f}%)")

    # write combined FASTA
    fasta_path = output_dir / "chemistry_inversion_all.fasta"
    with open(fasta_path, "w") as f:
        SeqIO.write(all_records, f, "fasta")
    print(f"\n  Wrote {len(all_records)} sequences to {fasta_path}")

    # write stats
    stats_df = pd.DataFrame(stats)
    stats_path = output_dir / "experiment_stats.csv"
    stats_df.to_csv(stats_path, index=False)
    print(f"  Wrote stats to {stats_path}")

    # summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"Total proteins: {len(stats)}")
    print(f"Total sequences: {len(all_records)} (5 conditions each)")
    print(f"Output FASTA: {fasta_path}")
    print()
    print("CONDITIONS:")
    print("  1. NATIVE - original sequence (positive control)")
    print("  2. EFS_PRESERVED - core native, non-core reduced")
    print("  3. ANTI_EFS_REDUCED - core reduced, non-core native")
    print("  4. ANTI_EFS_INVERTED - core chemistry-inverted (KEY TEST)")
    print("  5. BLIND_REDUCED - all reduced")
    print()
    print("NEXT STEP:")
    print(f"  trill --fold esmfold {fasta_path}")
    print("=" * 70)


if __name__ == "__main__":
    main()
