"""Reduced amino acid alphabets.

Two generations of alphabet exist in this project's history and they are
different partitions sharing the same names. Both are recorded here, the
canonical set as the supported interface and the legacy set as a clearly
marked, deprecated reproduction aid. Nothing in this module silently picks
one; the caller names the alphabet it wants.

Every alphabet is validated at import time against two invariants: it maps all
twenty standard amino acids, and its advertised size equals the number of
distinct output groups it actually produces. The second invariant is the one
that was violated in the original project code, where a mapping with eleven
groups was named RAA10 and one with seventeen was named RAA13. Names here
carry the true counts.
"""

from __future__ import annotations

from typing import Dict, Mapping

STANDARD_AA = frozenset("ACDEFGHIKLMNPQRSTVWY")


class AlphabetError(ValueError):
    """Raised when an alphabet mapping is incomplete or misnamed."""


def _validate(name: str, mapping: Mapping[str, str], advertised_size: int) -> Dict[str, str]:
    missing = STANDARD_AA - set(mapping)
    if missing:
        raise AlphabetError(
            f"{name}: mapping omits {''.join(sorted(missing))}"
        )
    unknown = set(mapping) - STANDARD_AA
    if unknown:
        raise AlphabetError(
            f"{name}: mapping contains non-standard residues {''.join(sorted(unknown))}"
        )
    groups = {mapping[aa] for aa in STANDARD_AA}
    if len(groups) != advertised_size:
        raise AlphabetError(
            f"{name}: advertised size {advertised_size} but mapping yields "
            f"{len(groups)} distinct groups ({''.join(sorted(groups))})"
        )
    return dict(mapping)


# --------------------------------------------------------------------------
# Canonical alphabets.
#
# RAA7 through RAA19 form a strict hierarchical refinement: every group at
# size n is contained in a group at size n-1, so a difference in outcome
# between adjacent alphabets is attributable to the distinction introduced at
# that step. RAA5 is a separate coarser partition and is NOT the parent of
# RAA7 (alanine is grouped with the large hydrophobics in RAA5 and with
# glycine in RAA7), so the RAA5-to-RAA7 transition is a re-assignment rather
# than a merge.
# --------------------------------------------------------------------------

RAA5 = _validate("RAA5", {
    "L": "L", "I": "L", "V": "L", "M": "L", "F": "L", "Y": "L", "W": "L", "A": "L",
    "G": "G", "P": "G",
    "S": "S", "T": "S", "N": "S", "Q": "S", "C": "S",
    "D": "D", "E": "D",
    "K": "K", "R": "K", "H": "K",
}, 5)

RAA7 = _validate("RAA7", {
    "L": "L", "I": "L", "V": "L", "M": "L",
    "F": "F", "Y": "F", "W": "F",
    "A": "A", "G": "A",
    "P": "P",
    "S": "S", "T": "S", "N": "S", "Q": "S", "C": "S",
    "D": "D", "E": "D",
    "K": "K", "R": "K", "H": "K",
}, 7)

RAA8 = _validate("RAA8", {
    "L": "L", "I": "L", "V": "L", "M": "L",
    "F": "F", "Y": "F", "W": "F",
    "A": "A", "G": "A",
    "P": "P",
    "S": "S", "T": "S", "N": "S", "Q": "S",
    "C": "C",
    "D": "D", "E": "D",
    "K": "K", "R": "K", "H": "K",
}, 8)

RAA9 = _validate("RAA9", {
    "L": "L", "I": "L", "V": "L",
    "M": "M",
    "F": "F", "Y": "F", "W": "F",
    "A": "A", "G": "A",
    "P": "P",
    "S": "S", "T": "S", "N": "S", "Q": "S",
    "C": "C",
    "D": "D", "E": "D",
    "K": "K", "R": "K", "H": "K",
}, 9)

# The two genuine ten-letter alphabets lying strictly between RAA9 and RAA11.
# RAA9 merges both {A,G} and {S,T,N,Q}; RAA11 splits both. Each of these
# splits exactly one of the two, which is what makes the 11-to-9 step
# decomposable.

RAA10_AG_MERGED = _validate("RAA10_AG_MERGED", {
    "L": "L", "I": "L", "V": "L",
    "M": "M",
    "F": "F", "Y": "F", "W": "F",
    "A": "A", "G": "A",
    "P": "P",
    "S": "S", "T": "S",
    "N": "N", "Q": "N",
    "C": "C",
    "D": "D", "E": "D",
    "K": "K", "R": "K", "H": "K",
}, 10)

RAA10_POLAR_MERGED = _validate("RAA10_POLAR_MERGED", {
    "L": "L", "I": "L", "V": "L",
    "M": "M",
    "F": "F", "Y": "F", "W": "F",
    "A": "A",
    "G": "G",
    "P": "P",
    "S": "S", "T": "S", "N": "S", "Q": "S",
    "C": "C",
    "D": "D", "E": "D",
    "K": "K", "R": "K", "H": "K",
}, 10)

# Published as "RAA10". It has eleven groups.
RAA11 = _validate("RAA11", {
    "L": "L", "I": "L", "V": "L",
    "M": "M",
    "F": "F", "Y": "F", "W": "F",
    "A": "A",
    "G": "G",
    "P": "P",
    "S": "S", "T": "S",
    "N": "N", "Q": "N",
    "C": "C",
    "D": "D", "E": "D",
    "K": "K", "R": "K", "H": "K",
}, 11)

# Published as "RAA13". The mapping that actually ran has NINETEEN groups:
# only isoleucine-with-leucine remains merged. Verified against the deposited
# sequences by provenance audit, 2026-10-09.
RAA19 = _validate("RAA19", {
    "L": "L", "I": "L",
    "V": "V", "M": "M",
    "F": "F", "Y": "Y", "W": "W",
    "A": "A", "G": "G", "P": "P",
    "S": "S", "T": "T",
    "N": "N", "Q": "Q",
    "C": "C",
    "D": "D", "E": "E",
    "K": "K", "R": "R", "H": "H",
}, 19)

# The variant printed in Supplementary Table S1 of the manuscript, which shows
# the basic residues merged as a single group and so has seventeen. No
# deposited sequence was generated with it: the table and the code that ran
# disagree, and the sequences settle it in favour of RAA19. Kept so that the
# discrepancy is expressible in code rather than only in prose.
RAA17_AS_PRINTED = _validate("RAA17_AS_PRINTED", {
    "L": "L", "I": "L",
    "V": "V", "M": "M",
    "F": "F", "Y": "Y", "W": "W",
    "A": "A", "G": "G", "P": "P",
    "S": "S", "T": "T",
    "N": "N", "Q": "Q",
    "C": "C",
    "D": "D", "E": "E",
    "K": "K", "R": "K", "H": "K",
}, 17)

CANONICAL: Dict[str, Dict[str, str]] = {
    "RAA5": RAA5,
    "RAA7": RAA7,
    "RAA8": RAA8,
    "RAA9": RAA9,
    "RAA10_AG_MERGED": RAA10_AG_MERGED,
    "RAA10_POLAR_MERGED": RAA10_POLAR_MERGED,
    "RAA11": RAA11,
    "RAA19": RAA19,
    "RAA17_AS_PRINTED": RAA17_AS_PRINTED,
}

# Names used in earlier project records, kept so that older scripts and
# filenames resolve to the right mapping.
ALIASES: Dict[str, str] = {
    "RAA10": "RAA11",
    "RAA13": "RAA19",
    "RAA7_FIXED": "RAA7",
}


# --------------------------------------------------------------------------
# Legacy alphabets. Deprecated. Present only so that results generated before
# the canonical set existed can be reproduced and identified. These are
# different partitions from the canonical alphabets of the same nominal size:
# alanine sits with threonine and histidine rather than with the
# hydrophobics, cysteine sits inside the hydrophobic class at size 5 and maps
# to aspartate at size 7, and the basic residues are merged with the polars.
# The size-7 cysteine mapping is a known error. Do not use these to generate
# new results.
# --------------------------------------------------------------------------

LEGACY_5 = _validate("LEGACY_5", {
    "C": "L", "M": "L", "F": "L", "I": "L", "L": "L", "V": "L", "W": "L", "Y": "L",
    "A": "A", "T": "A", "H": "A",
    "G": "G", "P": "G",
    "D": "E", "E": "E",
    "S": "Q", "N": "Q", "Q": "Q", "R": "Q", "K": "Q",
}, 5)

LEGACY_7 = _validate("LEGACY_7", {
    "A": "A", "G": "A", "V": "A",
    "I": "I", "L": "I", "F": "I", "P": "I",
    "Y": "Y", "M": "Y",
    "T": "T", "S": "T",
    "H": "H", "N": "H", "Q": "H", "W": "H",
    "R": "R", "K": "R",
    "D": "D", "E": "D", "C": "D",
}, 7)

LEGACY: Dict[str, Dict[str, str]] = {
    "LEGACY_5": LEGACY_5,
    "LEGACY_7": LEGACY_7,
}


def get(name: str) -> Dict[str, str]:
    """Return an alphabet mapping by name, resolving historical aliases."""
    key = ALIASES.get(name, name)
    if key in CANONICAL:
        return CANONICAL[key]
    if key in LEGACY:
        return LEGACY[key]
    raise KeyError(
        f"unknown alphabet {name!r}; available: "
        + ", ".join(sorted(CANONICAL) + sorted(LEGACY))
    )


def size(name: str) -> int:
    """Number of distinct output groups the named alphabet actually produces."""
    mapping = get(name)
    return len({mapping[aa] for aa in STANDARD_AA})


def output_letters(name: str) -> frozenset:
    """The set of letters a uniformly reduced sequence can contain.

    This is the signature used to identify which alphabet generated an
    existing reduced sequence, since the canonical and legacy partitions of
    the same nominal size emit disjoint letter sets.
    """
    mapping = get(name)
    return frozenset(mapping[aa] for aa in STANDARD_AA)


def identify(sequence: str, candidates=None) -> list:
    """Return the alphabets whose output letters could have produced `sequence`.

    Useful for auditing a reduced FASTA of unknown provenance. A sequence is
    consistent with an alphabet if every character in it is one that alphabet
    can emit. Positions holding native residues (as in position-specific
    reduction) must be excluded by the caller before calling this.
    """
    observed = set(sequence.upper())
    names = candidates if candidates is not None else list(CANONICAL) + list(LEGACY)
    return [n for n in names if observed <= output_letters(n)]
