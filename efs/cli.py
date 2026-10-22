"""Command-line interface for the Essential Folding Score.

    efs score      STRUCTURE...              per-residue scores as TSV or CSV
    efs reduce     STRUCTURE...              reduced sequences as FASTA
    efs control    STRUCTURE...              size-matched random control FASTA
    efs alphabets                            list alphabets and true sizes
    efs identify   REDUCED.fasta NATIVE.fasta   which alphabet produced a file
"""

from __future__ import annotations

import argparse
import csv
import sys
from dataclasses import asdict
from typing import Dict, Iterable, List, Optional, Sequence, TextIO, Tuple

from . import __version__
from .alphabets import ALIASES, CANONICAL, LEGACY, output_letters, size as alphabet_size
from .reduce import (
    identify_alphabet_from_reduction,
    random_control_matching,
    reduce_by_efs,
    reduce_uniform,
)
from .score import (
    DEFAULT_CORE_THRESHOLD,
    PROFILES,
    PUBLISHED,
    ScoredChain,
    StructureError,
    score_structure,
)


def _record_id(path: str) -> str:
    import os

    return os.path.splitext(os.path.basename(path))[0]


def _read_fasta(handle: TextIO) -> List[Tuple[str, str]]:
    records: List[Tuple[str, str]] = []
    name: Optional[str] = None
    chunks: List[str] = []
    for line in handle:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                records.append((name, "".join(chunks).upper()))
            name = line[1:].split()[0] if len(line) > 1 else ""
            chunks = []
        else:
            chunks.append(line)
    if name is not None:
        records.append((name, "".join(chunks).upper()))
    return records


def _score_many(
    paths: Sequence[str],
    chain: Optional[str],
    profile_name: str,
    use_ss: bool,
    stderr: TextIO,
) -> List[Tuple[str, ScoredChain]]:
    profile = PROFILES[profile_name]
    out: List[Tuple[str, ScoredChain]] = []
    for path in paths:
        try:
            out.append(
                (
                    _record_id(path),
                    score_structure(
                        path,
                        chain_id=chain,
                        profile=profile,
                        use_secondary_structure=use_ss,
                    ),
                )
            )
        except (StructureError, OSError, ValueError) as exc:
            print(f"{path}: skipped ({exc})", file=stderr)
    return out


def cmd_score(args: argparse.Namespace, stdout: TextIO, stderr: TextIO) -> int:
    scored = _score_many(args.structure, args.chain, args.profile, args.dssp, stderr)
    if not scored:
        return 1
    delimiter = "," if args.format == "csv" else "\t"
    writer = csv.writer(stdout, delimiter=delimiter, lineterminator="\n")
    writer.writerow(
        [
            "protein_id",
            "chain_id",
            "res_index",
            "pdb_resnum",
            "aa",
            "sasa",
            "contacts",
            "secondary_structure",
            "raw_score",
            "efs",
            "is_core",
        ]
    )
    for record_id, chain in scored:
        for residue in chain.residues:
            row = asdict(residue)
            writer.writerow(
                [
                    record_id,
                    row["chain_id"],
                    row["index"],
                    row["pdb_resnum"],
                    row["aa"],
                    f"{row['sasa']:.2f}",
                    row["contacts"],
                    row["secondary_structure"] or "",
                    row["raw_score"],
                    row["score"],
                    int(residue.is_core(args.threshold)),
                ]
            )
    if args.summary:
        for record_id, chain in scored:
            n = len(chain.core_positions(args.threshold))
            print(
                f"{record_id}\t{len(chain)} residues\t{n} core "
                f"({100 * chain.core_fraction(args.threshold):.1f}%) "
                f"at C{args.threshold}",
                file=stderr,
            )
    return 0


def cmd_reduce(args: argparse.Namespace, stdout: TextIO, stderr: TextIO) -> int:
    scored = _score_many(args.structure, args.chain, args.profile, args.dssp, stderr)
    if not scored:
        return 1
    for record_id, chain in scored:
        if args.blind:
            result = reduce_uniform(chain.sequence, args.alphabet)
        else:
            result = reduce_by_efs(chain, args.alphabet, threshold=args.threshold)
        stdout.write(result.fasta(record_id))
        if args.stats:
            print(
                f"{record_id}\t{result.length}\t{len(result.preserved)}\t"
                f"{result.preserved_fraction:.4f}\t"
                f"{result.percent_identity_to_native:.2f}",
                file=stderr,
            )
    return 0


def cmd_control(args: argparse.Namespace, stdout: TextIO, stderr: TextIO) -> int:
    scored = _score_many(args.structure, args.chain, args.profile, args.dssp, stderr)
    if not scored:
        return 1
    for seed in args.seed:
        for record_id, chain in scored:
            result = random_control_matching(
                chain, args.alphabet, seed=seed, threshold=args.threshold
            )
            stdout.write(result.fasta(record_id))
            if args.stats:
                print(
                    f"{record_id}\t{seed}\t{result.length}\t"
                    f"{len(result.preserved)}\t"
                    f"{result.percent_identity_to_native:.2f}",
                    file=stderr,
                )
    return 0


def cmd_alphabets(args: argparse.Namespace, stdout: TextIO, stderr: TextIO) -> int:
    def describe(name: str) -> str:
        letters = "".join(sorted(output_letters(name)))
        return f"{name:<20} {alphabet_size(name):>3} groups   emits {letters}"

    print("Canonical", file=stdout)
    for name in CANONICAL:
        print("  " + describe(name), file=stdout)
    print("\nHistorical names", file=stdout)
    for old, new in ALIASES.items():
        print(f"  {old:<20} -> {new}", file=stdout)
    print("\nLegacy (deprecated; for reproducing older results only)", file=stdout)
    for name in LEGACY:
        print("  " + describe(name), file=stdout)
    return 0


def cmd_identify(args: argparse.Namespace, stdout: TextIO, stderr: TextIO) -> int:
    with open(args.reduced) as handle:
        reduced = _read_fasta(handle)
    with open(args.native) as handle:
        native = dict(_read_fasta(handle))
    if not reduced:
        print(f"{args.reduced}: no records", file=stderr)
        return 1

    unresolved = 0
    for name, seq in reduced:
        match = native.get(name)
        if match is None:
            # Reduced records are often named <id>_<condition>; fall back to
            # the longest native id that prefixes this one.
            prefixes = [k for k in native if name.startswith(k)]
            match = native[max(prefixes, key=len)] if prefixes else None
        if match is None or len(match) != len(seq):
            print(f"{name}\tno matching native sequence", file=stdout)
            unresolved += 1
            continue
        candidates = identify_alphabet_from_reduction(seq, match)
        print(f"{name}\t{','.join(candidates) if candidates else 'none'}", file=stdout)
    return 1 if unresolved == len(reduced) else 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="efs",
        description=(
            "Essential Folding Score: per-residue structural criticality, and "
            "position-specific reduced-alphabet sequences derived from it."
        ),
    )
    parser.add_argument("--version", action="version", version=f"efs {__version__}")
    sub = parser.add_subparsers(dest="command", required=True)

    def add_scoring_options(p: argparse.ArgumentParser) -> None:
        p.add_argument("structure", nargs="+", help="PDB or mmCIF files")
        p.add_argument("--chain", default=None, help="chain id (default: first protein chain)")
        p.add_argument(
            "--profile",
            default=PUBLISHED.name,
            choices=sorted(PROFILES),
            help=(
                "published reproduces the deposited dataset; "
                "methods_as_written follows the paper's Methods text, which "
                "differs at proline and histidine"
            ),
        )
        p.add_argument(
            "--dssp",
            action="store_true",
            help="add a secondary-structure point (needs the DSSP binary; off in the published scores)",
        )
        p.add_argument(
            "--threshold",
            type=int,
            default=DEFAULT_CORE_THRESHOLD,
            help=f"core score cutoff (default {DEFAULT_CORE_THRESHOLD})",
        )

    p_score = sub.add_parser("score", help="write per-residue scores")
    add_scoring_options(p_score)
    p_score.add_argument("--format", choices=("tsv", "csv"), default="tsv")
    p_score.add_argument("--summary", action="store_true", help="per-protein core counts on stderr")
    p_score.set_defaults(func=cmd_score)

    p_reduce = sub.add_parser("reduce", help="write reduced sequences as FASTA")
    add_scoring_options(p_reduce)
    p_reduce.add_argument("--alphabet", default="RAA5", help="alphabet name (see `efs alphabets`)")
    p_reduce.add_argument(
        "--blind",
        action="store_true",
        help="reduce every position uniformly instead of sparing the core",
    )
    p_reduce.add_argument("--stats", action="store_true", help="per-protein statistics on stderr")
    p_reduce.set_defaults(func=cmd_reduce)

    p_control = sub.add_parser(
        "control",
        help="random-position control matched per protein to the EFS core size",
    )
    add_scoring_options(p_control)
    p_control.add_argument("--alphabet", default="RAA5")
    p_control.add_argument(
        "--seed", type=int, nargs="+", default=[1], help="one output set per seed"
    )
    p_control.add_argument("--stats", action="store_true")
    p_control.set_defaults(func=cmd_control)

    p_alpha = sub.add_parser("alphabets", help="list alphabets with their true group counts")
    p_alpha.set_defaults(func=cmd_alphabets)

    p_ident = sub.add_parser(
        "identify",
        help="report which alphabet could have produced a reduced FASTA",
    )
    p_ident.add_argument("reduced", help="FASTA of reduced sequences")
    p_ident.add_argument("native", help="FASTA of the corresponding native sequences")
    p_ident.set_defaults(func=cmd_identify)

    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return args.func(args, sys.stdout, sys.stderr)


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
