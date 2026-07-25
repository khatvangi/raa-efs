#!/usr/bin/env python3
"""
generate_core_definition_controls.py

Test which core DEFINITION matters for structural stability:
1. SASA_ONLY - burial alone (positions with low SASA)
2. HYDRO_ONLY - hydrophobic residues only
3. COMBINED - current method (SASA + contacts + hydrophobicity bonuses)

For each definition, apply EFS_PRESERVED (RAA10 at non-core positions).
"""

import pandas as pd
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# RAA10 alphabet - the reduction that works
RAA10_MAP = {
    'L': 'L', 'I': 'L', 'V': 'L', 'M': 'L',
    'F': 'F', 'Y': 'F', 'W': 'F',
    'A': 'A',
    'G': 'G',
    'P': 'P',
    'S': 'S', 'T': 'S',
    'C': 'C',
    'N': 'N', 'Q': 'N',
    'D': 'D', 'E': 'D',
    'K': 'K', 'R': 'K',
    'H': 'H',
}

# the 30 proteins from subtle_chemistry
PROTEINS = [
    'A2QRF6', 'A3QK16', 'A5A6L6', 'O35129', 'O51221', 'O75716', 'O76484',
    'O94550', 'P06907', 'P0AA56', 'P0CP66', 'P33681', 'P35232', 'P36887',
    'P40199', 'P9WPR9', 'Q07699', 'Q0J4I1', 'Q149L7', 'Q21190', 'Q2HJ97',
    'Q2YHT5', 'Q4U925', 'Q5ZMN3', 'Q6AYD4', 'Q6NW76', 'Q8N165', 'Q8TDQ1',
    'Q9H6B4', 'Q9N0X0'
]

def load_core_scores_per_protein(base_dir, family, protein_id):
    """load per-protein core scores from family_EFS directory"""
    path = f"{base_dir}/{family}/family_EFS/per_protein/{protein_id}_core_scores.csv"
    if os.path.exists(path):
        return pd.read_csv(path)
    return None

def get_sequence(df, protein_id):
    """extract sequence from core scores dataframe"""
    prot_df = df[df['protein_id'] == protein_id].sort_values('res_index')
    return ''.join(prot_df['aa'].tolist())

def define_core_sasa_only(prot_df, sasa_threshold):
    """core = positions with SASA < threshold (buried)"""
    core_positions = set(prot_df[prot_df['sasa'] < sasa_threshold]['res_index'])
    return core_positions

def define_core_hydro_only(prot_df):
    """core = hydrophobic positions only"""
    core_positions = set(prot_df[prot_df['hydrophobic'] == 1]['res_index'])
    return core_positions

def define_core_combined(prot_df, score_threshold):
    """core = positions with core_score >= threshold (current method)"""
    core_positions = set(prot_df[prot_df['core_score'] >= score_threshold]['res_index'])
    return core_positions

def apply_reduction(sequence, core_positions):
    """apply RAA10 to non-core positions, preserve core"""
    result = []
    for i, aa in enumerate(sequence):
        pos = i + 1  # 1-indexed
        if pos in core_positions:
            # preserve core
            result.append(aa)
        else:
            # reduce non-core
            result.append(RAA10_MAP.get(aa, aa))
    return ''.join(result)

def main():
    base_dir = '/storage/kiran-stuff/RAA'
    output_dir = '/storage/kiran-stuff/RAA/core_definition_controls'
    os.makedirs(output_dir, exist_ok=True)

    # protein -> family mapping
    protein_families = {
        "PF00069": ["P0CP66", "Q6NW76", "Q9N0X0", "O75716", "P36887",
                    "Q0J4I1", "A2QRF6", "Q8N165", "O76484", "Q4U925"],
        "PF01145": ["P9WPR9", "Q2HJ97", "P35232", "O51221", "Q21190",
                    "A3QK16", "Q5ZMN3", "O94550", "O35129", "P0AA56"],
        "PF07686": ["P33681", "Q9H6B4", "P06907", "Q6AYD4", "Q2YHT5",
                    "Q149L7", "Q8TDQ1", "Q07699", "A5A6L6", "P40199"]
    }

    all_records = []
    stats = []

    # core definitions to test
    definitions = [
        ('SASA_25', 'sasa', 25),      # SASA < 25: ~30% core
        ('SASA_40', 'sasa', 40),      # SASA < 40: ~47% core
        ('HYDRO', 'hydro', None),     # hydrophobic only: ~46% core
        ('COMBINED_20', 'combined', 20),  # score >= 20: ~37% core
        ('COMBINED_12', 'combined', 12),  # score >= 12: ~46% core
    ]

    for family, protein_list in protein_families.items():
        print(f"\n--- {family} ---")
        for protein_id in protein_list:
            prot_df = load_core_scores_per_protein(base_dir, family, protein_id)
            if prot_df is None:
                print(f"  WARNING: {protein_id} not found")
                continue

            prot_df = prot_df.sort_values('res_index')
            sequence = ''.join(prot_df['aa'].tolist())
            seq_len = len(sequence)

            # add native
            rec = SeqRecord(Seq(sequence), id=f"{protein_id}_NATIVE", description="")
            all_records.append(rec)

            for def_name, def_type, threshold in definitions:
                if def_type == 'sasa':
                    core_pos = define_core_sasa_only(prot_df, threshold)
                elif def_type == 'hydro':
                    core_pos = define_core_hydro_only(prot_df)
                else:  # combined
                    core_pos = define_core_combined(prot_df, threshold)

                reduced_seq = apply_reduction(sequence, core_pos)

                # count changes
                n_changed = sum(1 for a, b in zip(sequence, reduced_seq) if a != b)
                core_pct = len(core_pos) / seq_len * 100

                rec = SeqRecord(Seq(reduced_seq),
                              id=f"{protein_id}_{def_name}",
                              description="")
                all_records.append(rec)

                stats.append({
                    'protein': protein_id,
                    'family': family,
                    'definition': def_name,
                    'core_size': len(core_pos),
                    'core_pct': core_pct,
                    'positions_changed': n_changed,
                    'change_pct': n_changed / seq_len * 100
                })

            print(f"  {protein_id}: {seq_len} AA, 5 definitions generated")

    # write fasta
    fasta_path = os.path.join(output_dir, 'core_definition_controls.fasta')
    SeqIO.write(all_records, fasta_path, 'fasta')
    print(f"Wrote {len(all_records)} sequences to {fasta_path}")

    # write stats
    stats_df = pd.DataFrame(stats)
    stats_path = os.path.join(output_dir, 'core_definition_stats.csv')
    stats_df.to_csv(stats_path, index=False)

    # summary
    print("\n=== SUMMARY BY DEFINITION ===")
    summary = stats_df.groupby('definition').agg({
        'core_pct': 'mean',
        'change_pct': 'mean'
    }).round(1)
    print(summary)

if __name__ == '__main__':
    main()
