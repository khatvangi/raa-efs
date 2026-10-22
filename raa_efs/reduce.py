"""Position-specific and uniform reduction of protein sequences.

Uniform ("blind") reduction maps every position through a reduced alphabet.
Position-specific reduction retains native identity at a chosen set of
positions and reduces the rest. The set may come from the Essential Folding
Score, from any externally supplied list of positions, or — for controls —
from a seeded random draw matched in size to another set.

The random control exists here as a first-class function rather than as a
one-off script because it is the baseline that makes a position-specific
result interpretable: preserving a quarter of positions necessarily moves the
variant closer to its native sequence, so a rescue is only attributable to
*which* positions were chosen if a size-matched random choice does worse.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Sequence, Set

from .alphabets import STANDARD_AA, get as get_alphabet
from .score import DEFAULT_CORE_THRESHOLD, ScoredChain


@dataclass(frozen=True)
class ReductionResult:
    """A reduced sequence together with what was done to produce it."""

    sequence: str
    native: str
    alphabet: str
    preserved: frozenset
    label: str = ""

    @property
    def length(self) -> int:
        return len(self.native)

    @property
    def preserved_fraction(self) -> float:
        return len(self.preserved) / self.length if self.length else 0.0

    @property
    def positions_changed(self) -> int:
        return sum(1 for a, b in zip(self.native, self.sequence) if a != b)

    @property
    def percent_identity_to_native(self) -> float:
        """Percentage of positions whose residue is unchanged.

        Higher than `preserved_fraction`, because a reduced position whose
        native residue is its own group representative is unchanged too. This
        is the quantity that separates "the variant is closer to native" from
        "the right positions were protected", and it should be reported
        alongside any position-specific result.
        """
        if not self.length:
            return 0.0
        same = sum(1 for a, b in zip(self.native, self.sequence) if a == b)
        return 100.0 * same / self.length

    def fasta(self, record_id: str, width: int = 60) -> str:
        name = record_id if not self.label else f"{record_id}_{self.label}"
        lines = [f">{name}"]
        lines += [
            self.sequence[i : i + width] for i in range(0, len(self.sequence), width)
        ]
        return "\n".join(lines) + "\n"


def _check_sequence(sequence: str) -> str:
    seq = sequence.strip().upper()
    if not seq:
        raise ValueError("empty sequence")
    unknown = set(seq) - STANDARD_AA
    if unknown:
        raise ValueError(
            "sequence contains non-standard residues "
            f"{''.join(sorted(unknown))}; remove or substitute them first"
        )
    return seq


def _check_positions(positions: Iterable[int], length: int) -> frozenset:
    pos = frozenset(int(p) for p in positions)
    out_of_range = {p for p in pos if p < 1 or p > length}
    if out_of_range:
        raise ValueError(
            f"positions out of range for a sequence of length {length}: "
            f"{sorted(out_of_range)[:5]}"
        )
    return pos


def reduce_uniform(sequence: str, alphabet: str) -> ReductionResult:
    """Map every position through `alphabet`."""
    seq = _check_sequence(sequence)
    mapping = get_alphabet(alphabet)
    return ReductionResult(
        sequence="".join(mapping[aa] for aa in seq),
        native=seq,
        alphabet=alphabet,
        preserved=frozenset(),
        label=f"blind_{alphabet}",
    )


def reduce_position_specific(
    sequence: str,
    alphabet: str,
    preserved_positions: Iterable[int],
    label: str = "",
) -> ReductionResult:
    """Keep native residues at `preserved_positions` (1-based), reduce the rest."""
    seq = _check_sequence(sequence)
    preserved = _check_positions(preserved_positions, len(seq))
    mapping = get_alphabet(alphabet)
    reduced = "".join(
        aa if (i + 1) in preserved else mapping[aa] for i, aa in enumerate(seq)
    )
    return ReductionResult(
        sequence=reduced,
        native=seq,
        alphabet=alphabet,
        preserved=preserved,
        label=label or f"preserved_{alphabet}",
    )


def reduce_by_efs(
    scored: ScoredChain,
    alphabet: str,
    threshold: int = DEFAULT_CORE_THRESHOLD,
    label: str = "",
) -> ReductionResult:
    """Reduce every position except those whose EFS meets `threshold`."""
    return reduce_position_specific(
        scored.sequence,
        alphabet,
        scored.core_positions(threshold),
        label=label or f"EFS_C{threshold}_{alphabet}",
    )


def random_matched_control(
    sequence: str,
    alphabet: str,
    n_preserved: int,
    seed: int,
    label: str = "",
) -> ReductionResult:
    """Preserve `n_preserved` uniformly random positions and reduce the rest.

    The seed is required, not optional, so that every control sequence is
    reproducible from its record alone.
    """
    seq = _check_sequence(sequence)
    if not 0 <= n_preserved <= len(seq):
        raise ValueError(
            f"n_preserved={n_preserved} outside [0, {len(seq)}]"
        )
    try:
        import numpy as np
    except ImportError as exc:  # pragma: no cover
        raise ImportError(
            "random_matched_control requires numpy; install it with "
            "`pip install numpy`"
        ) from exc
    rng = np.random.default_rng(seed)
    chosen = rng.choice(len(seq), size=n_preserved, replace=False)
    positions = {int(i) + 1 for i in chosen}
    return reduce_position_specific(
        seq, alphabet, positions, label=label or f"RANDOM{seed}_{alphabet}"
    )


def random_control_matching(
    scored: ScoredChain,
    alphabet: str,
    seed: int,
    threshold: int = DEFAULT_CORE_THRESHOLD,
    label: str = "",
) -> ReductionResult:
    """Random control sized to this chain's own EFS core, position by position.

    Matching the count per protein rather than using a global percentage is
    what makes the comparison fair: core fraction varies from under 15% to
    over 75% across families, and a fixed global fraction would compare
    different amounts of preserved chemistry.
    """
    n = len(scored.core_positions(threshold))
    return random_matched_control(
        scored.sequence,
        alphabet,
        n_preserved=n,
        seed=seed,
        label=label or f"RANDOM{seed}_C{threshold}_{alphabet}",
    )


def identify_alphabet_from_reduction(
    result_sequence: str,
    native: str,
    preserved_positions: Optional[Iterable[int]] = None,
    candidates: Optional[Sequence[str]] = None,
) -> List[str]:
    """Which alphabets could have produced this reduced sequence.

    Positions that retained native identity are excluded before testing, so
    this works on position-specific reductions as well as uniform ones. With
    `preserved_positions` unknown, every position where the sequence matches
    the native is excluded, which is conservative: it can only widen the
    candidate set, never produce a false identification.
    """
    from .alphabets import identify

    if len(result_sequence) != len(native):
        raise ValueError("reduced sequence and native differ in length")
    if preserved_positions is None:
        tested = "".join(
            b for a, b in zip(native, result_sequence) if a != b
        )
    else:
        preserved = _check_positions(preserved_positions, len(native))
        tested = "".join(
            c for i, c in enumerate(result_sequence, start=1) if i not in preserved
        )
    return identify(tested, candidates)


def summarise(results: Dict[str, ReductionResult]) -> List[dict]:
    """Flatten a set of reductions into rows for a table."""
    rows = []
    for name, r in results.items():
        rows.append(
            {
                "id": name,
                "label": r.label,
                "alphabet": r.alphabet,
                "length": r.length,
                "n_preserved": len(r.preserved),
                "preserved_fraction": round(r.preserved_fraction, 4),
                "positions_changed": r.positions_changed,
                "percent_identity_to_native": round(r.percent_identity_to_native, 2),
            }
        )
    return rows
