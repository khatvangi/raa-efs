#!/usr/bin/env python3
"""
compare_core_native_vs_raa10.py

1. Compute core scores from EFS_PRESERVED (RAA10) PDBs
2. Compare with NATIVE core scores
3. Quantify: % overlap, which positions shifted

Question: Does the core shift after RAA10 reduction?
"""

import os
import glob
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, Tuple, List

from Bio.PDB import PDBParser, NeighborSearch, SASA
from Bio.PDB.Polypeptide import is_aa
from Bio.Data import IUPACData

# 3-letter to 1-letter mapping
_THREE_TO_ONE = {k.upper(): v for k, v in IUPACData.protein_letters_3to1.items()}
_THREE_TO_ONE.update({"MSE": "M", "SEC": "C", "PYL": "K"})

def three_to_one(resname: str) -> str:
    return _THREE_TO_ONE.get(resname.strip().upper(), "X")

def compute_core_scores_from_pdb(pdb_path: str) -> pd.DataFrame:
    """compute per-residue core scores from a PDB file"""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_path)
    model = structure[0]

    # find first chain with amino acids
    chain = None
    for c in model:
        if any(is_aa(res, standard=True) for res in c):
            chain = c
            break
    if chain is None:
        return None

    # compute SASA
    sr = SASA.ShrakeRupley()
    sr.compute(structure, level="R")

    # compute contacts
    ca_atoms = []
    for res in chain:
        if is_aa(res, standard=True) and "CA" in res:
            ca_atoms.append((res.id, res["CA"]))

    ns = NeighborSearch([atom for _, atom in ca_atoms])
    contact_counts = {}
    for res_id, atom in ca_atoms:
        neighbors = ns.search(atom.coord, 8.0)
        count = sum(1 for n in neighbors
                   if n.get_parent().id != res_id and is_aa(n.get_parent(), standard=True))
        contact_counts[res_id] = count

    # collect data
    rows = []
    res_index = 0
    for res in chain:
        if not is_aa(res, standard=True):
            continue
        res_index += 1

        aa = three_to_one(res.get_resname())
        sasa = getattr(res, "sasa", 999.0)
        contacts = contact_counts.get(res.id, 0)

        # flags
        hydrophobic = aa in "AVILMFWYP"
        charged = aa in "RKDE"
        cysteine = aa == "C"
        aromatic = aa in "FYWH"

        # compute raw score (same logic as original)
        raw = 0
        if sasa < 20.0:
            raw += 3
        elif sasa < 40.0:
            raw += 1
        if contacts >= 8:
            raw += 2
        elif contacts >= 5:
            raw += 1
        if hydrophobic and sasa < 30.0:
            raw += 2
        if cysteine:
            raw += 2
        if aromatic and sasa < 40.0:
            raw += 1

        core_score = raw * 4.0

        rows.append({
            'res_index': res_index,
            'aa': aa,
            'sasa': sasa,
            'contacts': contacts,
            'hydrophobic': int(hydrophobic),
            'core_score': core_score,
            'raw_score': raw
        })

    return pd.DataFrame(rows)


def load_native_core_scores(base_dir: Path, family: str, protein_id: str) -> pd.DataFrame:
    """load pre-computed native core scores"""
    path = base_dir / family / "family_EFS" / "per_protein" / f"{protein_id}_core_scores.csv"
    if path.exists():
        return pd.read_csv(path)
    return None


def compare_cores(native_df: pd.DataFrame, raa10_df: pd.DataFrame, threshold: int = 12):
    """compare core positions between native and RAA10"""
    if native_df is None or raa10_df is None:
        return None

    # get core positions (1-indexed)
    native_core = set(native_df[native_df['core_score'] >= threshold]['res_index'])
    raa10_core = set(raa10_df[raa10_df['core_score'] >= threshold]['res_index'])

    # compute overlap
    intersection = native_core & raa10_core
    union = native_core | raa10_core

    # positions that changed
    lost_from_core = native_core - raa10_core  # was core in native, not in RAA10
    gained_to_core = raa10_core - native_core  # was not core in native, is core in RAA10

    jaccard = len(intersection) / len(union) if union else 1.0
    overlap_native = len(intersection) / len(native_core) if native_core else 1.0
    overlap_raa10 = len(intersection) / len(raa10_core) if raa10_core else 1.0

    return {
        'native_core_size': len(native_core),
        'raa10_core_size': len(raa10_core),
        'intersection': len(intersection),
        'jaccard': jaccard,
        'overlap_with_native': overlap_native,
        'overlap_with_raa10': overlap_raa10,
        'lost_from_core': len(lost_from_core),
        'gained_to_core': len(gained_to_core),
        'lost_positions': sorted(lost_from_core),
        'gained_positions': sorted(gained_to_core)
    }


def main():
    base_dir = Path("/storage/kiran-stuff/RAA")

    # find all EFS_PRESERVED PDBs from core_size_variation
    efs_pdbs = glob.glob(str(base_dir / "core_size_variation" / "*_EFS_PRESERVED.pdb"))

    # protein -> family mapping
    protein_families = {
        "PF00069": ["P0CP66", "Q6NW76", "Q9N0X0", "O75716", "P36887",
                    "Q0J4I1", "A2QRF6", "Q8N165", "O76484", "Q4U925"],
        "PF01145": ["P9WPR9", "Q2HJ97", "P35232", "O51221", "Q21190",
                    "A3QK16", "Q5ZMN3", "O94550", "O35129", "P0AA56"],
        "PF07686": ["P33681", "Q9H6B4", "P06907", "Q6AYD4", "Q2YHT5",
                    "Q149L7", "Q8TDQ1", "Q07699", "A5A6L6", "P40199"]
    }

    # reverse lookup
    protein_to_family = {}
    for fam, prots in protein_families.items():
        for p in prots:
            protein_to_family[p] = fam

    results = []
    detailed_results = []

    print("=" * 70)
    print("COMPARING CORE: NATIVE vs RAA10 (EFS_PRESERVED)")
    print("=" * 70)

    # process only C12 cutoff for main comparison (most used)
    for pdb_path in sorted(efs_pdbs):
        filename = os.path.basename(pdb_path)
        # parse: {protein}_C{cutoff}_EFS_PRESERVED.pdb
        parts = filename.replace("_EFS_PRESERVED.pdb", "").split("_C")
        if len(parts) != 2:
            continue

        protein_id = parts[0]
        cutoff = int(parts[1])

        # only process C12 for now
        if cutoff != 12:
            continue

        family = protein_to_family.get(protein_id)
        if not family:
            continue

        # compute core from RAA10 structure
        raa10_df = compute_core_scores_from_pdb(pdb_path)
        if raa10_df is None:
            print(f"  {protein_id}: Failed to compute RAA10 core")
            continue

        # load native core
        native_df = load_native_core_scores(base_dir, family, protein_id)
        if native_df is None:
            print(f"  {protein_id}: Native core not found")
            continue

        # compare at multiple thresholds
        for thresh in [8, 12, 20, 28]:
            comp = compare_cores(native_df, raa10_df, threshold=thresh)
            if comp:
                results.append({
                    'protein': protein_id,
                    'family': family,
                    'threshold': thresh,
                    **{k: v for k, v in comp.items() if k not in ['lost_positions', 'gained_positions']}
                })

                if thresh == 12:  # detailed for C12
                    detailed_results.append({
                        'protein': protein_id,
                        'family': family,
                        'lost_positions': comp['lost_positions'],
                        'gained_positions': comp['gained_positions']
                    })

        print(f"  {protein_id}: native={len(native_df)}, raa10={len(raa10_df)} residues")

    # create summary dataframe
    df = pd.DataFrame(results)

    # output directory
    output_dir = base_dir / "core_native_vs_raa10"
    output_dir.mkdir(exist_ok=True)

    df.to_csv(output_dir / "core_comparison_results.csv", index=False)

    # summary by threshold
    print("\n" + "=" * 70)
    print("SUMMARY BY THRESHOLD")
    print("=" * 70)

    for thresh in [8, 12, 20, 28]:
        sub = df[df['threshold'] == thresh]
        if len(sub) == 0:
            continue
        print(f"\n--- Threshold: core_score >= {thresh} ---")
        print(f"  Proteins analyzed: {len(sub)}")
        print(f"  Mean Jaccard similarity: {sub['jaccard'].mean():.3f}")
        print(f"  Mean overlap with native: {sub['overlap_with_native'].mean():.3f}")
        print(f"  Mean positions lost from core: {sub['lost_from_core'].mean():.1f}")
        print(f"  Mean positions gained to core: {sub['gained_to_core'].mean():.1f}")

    # detailed analysis at C12
    print("\n" + "=" * 70)
    print("DETAILED ANALYSIS AT C12")
    print("=" * 70)

    c12 = df[df['threshold'] == 12]

    # best and worst overlap
    if len(c12) > 0:
        best = c12.loc[c12['jaccard'].idxmax()]
        worst = c12.loc[c12['jaccard'].idxmin()]

        print(f"\nBest preservation: {best['protein']} (Jaccard={best['jaccard']:.3f})")
        print(f"Worst preservation: {worst['protein']} (Jaccard={worst['jaccard']:.3f})")

        print(f"\nOverall at C12:")
        print(f"  Mean native core size: {c12['native_core_size'].mean():.1f}")
        print(f"  Mean RAA10 core size: {c12['raa10_core_size'].mean():.1f}")
        print(f"  Mean Jaccard: {c12['jaccard'].mean():.3f}")

    # save detailed results
    pd.DataFrame(detailed_results).to_csv(output_dir / "core_shift_details.csv", index=False)

    print(f"\nResults saved to: {output_dir}")


if __name__ == "__main__":
    main()
