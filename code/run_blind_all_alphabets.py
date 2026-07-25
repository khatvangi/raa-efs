#!/usr/bin/env python3
"""
run_blind_all_alphabets.py - Blind (uniform) alphabet reduction runner

Applies any combination of {RAA5, RAA7, RAA8, RAA9, RAA10, RAA13} to every
residue in every sequence of a FASTA file, writing one reduced FASTA per
alphabet. All alphabets are imported from CANONICAL_ALPHABETS.py to
guarantee a single source of truth.

This replaces the older run_blind_13-10-7-5_fasta.py, which defined
alphabets inline (some of them buggy). Use only this script going forward.

Typical usage (from the PF00069 directory, for example):

    python run_blind_all_alphabets.py \\
        --input PF00069_natives.fasta \\
        --output_prefix PF00069 \\
        --alphabets 13 10 9 8 7 5

Output filenames follow the existing convention:
    <prefix>_blind_<N>aa.fasta    for each N in --alphabets

so that downstream pipelines (TRILL/ESMFold folding, USalign comparison)
pick them up without modification.

Created: 2026-04-18
"""

import argparse
import os
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# import canonical alphabets. adjust sys.path if running from a different
# directory; the standard project layout places CANONICAL_ALPHABETS.py at
# /storage/kiran-stuff/RAA/CANONICAL_CODE/.
_CANONICAL_DIR = "/storage/kiran-stuff/RAA/CANONICAL_CODE"
if _CANONICAL_DIR not in sys.path:
    sys.path.insert(0, _CANONICAL_DIR)

from CANONICAL_ALPHABETS import (
    RAA5, RAA7_FIXED, RAA8, RAA9, RAA10, RAA13,
    reduce_sequence, validate_alphabet,
)

# registry maps the CLI integer argument to (name, mapping dict).
# note that the dict-internal group counts differ slightly from the CLI
# names for RAA10 (11 actual output letters) and RAA13 (17 actual), but
# the CLI names match project convention and output filenames.
ALPHABET_REGISTRY = {
    5:  ("RAA5",       RAA5),
    7:  ("RAA7_FIXED", RAA7_FIXED),
    8:  ("RAA8",       RAA8),
    9:  ("RAA9",       RAA9),
    10: ("RAA10",      RAA10),
    13: ("RAA13",      RAA13),
}


def setup_arguments():
    p = argparse.ArgumentParser(
        description="Apply canonical blind alphabet reductions to FASTA sequences",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    p.add_argument("--input", "-i", required=True,
                   help="Input FASTA file (native sequences)")
    p.add_argument("--output_prefix", "-o", required=True,
                   help="Prefix for output files, e.g. 'PF00069' -> "
                        "'PF00069_blind_8aa.fasta'")
    p.add_argument("--alphabets", nargs="+", type=int,
                   choices=sorted(ALPHABET_REGISTRY.keys()),
                   default=[13, 10, 9, 8, 7, 5],
                   help="Alphabet sizes to generate "
                        "(default: 13 10 9 8 7 5)")
    p.add_argument("--suffix", default="_blind",
                   help="Suffix added to each FASTA record ID "
                        "(default: _blind)")
    p.add_argument("--skip_non_standard", action="store_true",
                   help="Skip sequences containing non-standard residues "
                        "(e.g. X, B, Z). Default is to fail-fast so "
                        "upstream FASTA cleaning errors are visible.")
    return p.parse_args()


def apply_reduction(input_fasta, output_file, alphabet_size, alphabet_name,
                    alphabet_map, header_suffix, skip_non_standard):
    """reduce every residue in every sequence; write to output_file."""
    print(f"\n--- {alphabet_name} ({alphabet_size}-letter) ---")
    print(f"  Output: {output_file}")

    reduced_records = []
    n_processed = 0
    n_errors = 0
    first_errors = []

    for record in SeqIO.parse(input_fasta, "fasta"):
        original_seq = str(record.seq).upper()
        try:
            reduced_seq = reduce_sequence(original_seq, alphabet_map)
            if "X" in reduced_seq and not skip_non_standard:
                bad = {aa for aa in original_seq if aa not in alphabet_map}
                raise ValueError(f"non-standard residues {bad} in "
                                 f"{record.id}")
            if "X" in reduced_seq and skip_non_standard:
                n_errors += 1
                if len(first_errors) < 5:
                    first_errors.append(
                        f"  skipped {record.id} (non-standard residues)"
                    )
                continue
        except Exception as exc:
            n_errors += 1
            if len(first_errors) < 5:
                first_errors.append(f"  error on {record.id}: {exc}")
            continue

        actual_alphabet_size = len(set(reduced_seq))
        new_id = (f"{record.id}{header_suffix}_{alphabet_size}aa"
                  f"_alphab{actual_alphabet_size}")
        new_record = SeqRecord(
            Seq(reduced_seq),
            id=new_id,
            description=(f"{alphabet_name} blind reduction "
                         f"(len={len(original_seq)}, "
                         f"unique_letters={actual_alphabet_size})"),
        )
        reduced_records.append(new_record)
        n_processed += 1

        if n_processed % 200 == 0:
            print(f"    processed {n_processed} sequences")

    # write FASTA
    out_dir = os.path.dirname(output_file)
    if out_dir and not os.path.isdir(out_dir):
        os.makedirs(out_dir, exist_ok=True)
    SeqIO.write(reduced_records, output_file, "fasta")

    print(f"  wrote {n_processed} sequences, skipped/errored {n_errors}")
    for line in first_errors:
        print(line)

    return n_processed, n_errors


def print_alphabet_summary(alphabet_name, alphabet_map):
    """show the group structure of an alphabet in reverse-map form."""
    reverse = {}
    for orig, reduced in alphabet_map.items():
        reverse.setdefault(reduced, []).append(orig)
    print(f"\n  {alphabet_name} groups:")
    for reduced_label in sorted(reverse):
        members = ",".join(sorted(reverse[reduced_label]))
        print(f"    {reduced_label} <- {members}")


def main():
    args = setup_arguments()

    if not os.path.isfile(args.input):
        print(f"Error: input file not found: {args.input}", file=sys.stderr)
        sys.exit(1)

    # validate all requested alphabets upfront
    for size in args.alphabets:
        name, mapping = ALPHABET_REGISTRY[size]
        if not validate_alphabet(mapping):
            print(f"Error: {name} does not map all 20 standard AAs",
                  file=sys.stderr)
            sys.exit(2)

    print("=== BLIND ALPHABET REDUCTION (canonical) ===")
    print(f"Input:            {args.input}")
    print(f"Output prefix:    {args.output_prefix}")
    print(f"Alphabet sizes:   {sorted(args.alphabets, reverse=True)}")
    print(f"Header suffix:    {args.suffix}")
    print(f"Non-standard:     "
          f"{'skip' if args.skip_non_standard else 'error'}")
    print("=" * 50)

    results = {}
    for size in sorted(args.alphabets, reverse=True):
        name, mapping = ALPHABET_REGISTRY[size]
        print_alphabet_summary(name, mapping)
        out_file = f"{args.output_prefix}_blind_{size}aa.fasta"
        n_proc, n_err = apply_reduction(
            args.input, out_file,
            alphabet_size=size,
            alphabet_name=name,
            alphabet_map=mapping,
            header_suffix=args.suffix,
            skip_non_standard=args.skip_non_standard,
        )
        results[size] = (n_proc, n_err, out_file)

    # final summary
    print("\n" + "=" * 60)
    print("BLIND REDUCTION COMPLETE")
    print("=" * 60)
    for size in sorted(results.keys(), reverse=True):
        n_proc, n_err, out_file = results[size]
        total = n_proc + n_err
        rate = (100.0 * n_proc / total) if total else 0.0
        print(f"  RAA{size:<2d}: {n_proc:5d} sequences "
              f"({rate:5.1f}% of input) -> {out_file}")

    print("\nNext steps:")
    print("  1. Fold each FASTA with TRILL / ESMFold:")
    print("       trill ... esmfold --fasta <file>_blind_<N>aa.fasta")
    print("  2. Score predicted vs native structures with US-align.")
    print("  3. Compare mean TM-score across alphabet sizes to locate "
          "the RAA-N -> fold collapse cliff on this dataset.")


if __name__ == "__main__":
    main()
