#!/usr/bin/env python3
"""Check that this package reproduces the published per-residue scores.

The deposited dataset stores, for every residue, the SASA and contact count
that went into the score as well as the score itself. That makes the scoring
function checkable without recomputing any structure: read the stored
features, re-run the arithmetic, and compare.

Run it against the per-protein score files, for example::

    python tools/validate_against_published.py \\
        /storage/kiran-stuff/RAA/PF00069/family_EFS/per_protein/*_core_scores.csv

Exit status is 0 when every row under at least one scoring profile matches,
1 otherwise. Where the two profiles disagree, the residues responsible are
listed: those are the proline and histidine positions where the original
implementation and the paper's Methods text part company.
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter
from typing import Dict, Iterable, List, Sequence, Tuple

sys.path.insert(0, __file__.rsplit("/tools/", 1)[0])

from efs.score import PROFILES, raw_score  # noqa: E402

REQUIRED_COLUMNS = {"aa", "sasa", "contacts", "raw_score"}


def check_file(path: str, use_ss: bool) -> Dict[str, Counter]:
    """Return, per profile, a counter of matched/mismatched rows."""
    results = {name: Counter() for name in PROFILES}
    mismatched_residues = {name: Counter() for name in PROFILES}
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle)
        missing = REQUIRED_COLUMNS - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path}: missing columns {sorted(missing)}")
        for row in reader:
            aa = row["aa"].strip().upper()
            sasa = float(row["sasa"])
            contacts = int(float(row["contacts"]))
            ss = row.get("ss") if use_ss else None
            want = int(float(row["raw_score"]))
            for name, profile in PROFILES.items():
                got = raw_score(aa, sasa, contacts, ss, profile)
                if got == want:
                    results[name]["match"] += 1
                else:
                    results[name]["mismatch"] += 1
                    mismatched_residues[name][aa] += 1
    for name in results:
        results[name]["residues"] = mismatched_residues[name]
    return results


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv", nargs="+", help="per-protein *_core_scores.csv files")
    parser.add_argument(
        "--use-ss",
        action="store_true",
        help="include the secondary-structure point (the deposited scores do not)",
    )
    args = parser.parse_args(argv)

    totals = {name: Counter() for name in PROFILES}
    offenders = {name: Counter() for name in PROFILES}
    failed_files: List[Tuple[str, str]] = []

    for path in args.csv:
        try:
            per_file = check_file(path, args.use_ss)
        except (OSError, ValueError) as exc:
            failed_files.append((path, str(exc)))
            continue
        for name, counts in per_file.items():
            totals[name]["match"] += counts["match"]
            totals[name]["mismatch"] += counts["mismatch"]
            offenders[name].update(counts["residues"])

    for path, reason in failed_files:
        print(f"unreadable: {reason}", file=sys.stderr)

    ok = False
    for name in sorted(totals):
        counts = totals[name]
        n = counts["match"] + counts["mismatch"]
        if not n:
            continue
        print(
            f"{name:<20} {counts['match']}/{n} rows reproduce "
            f"({100 * counts['match'] / n:.3f}%)"
        )
        if counts["mismatch"]:
            worst = ", ".join(
                f"{aa}:{k}" for aa, k in offenders[name].most_common(6)
            )
            print(f"{'':<20} mismatches by residue: {worst}")
        else:
            ok = True

    if not totals or all(
        (c["match"] + c["mismatch"]) == 0 for c in totals.values()
    ):
        print("no rows were read", file=sys.stderr)
        return 1
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
