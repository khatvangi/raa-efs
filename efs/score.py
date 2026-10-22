"""Essential Folding Score: per-residue structural criticality from a structure.

The score is a small integer sum over burial, packing and residue chemistry,
multiplied by a scale factor. It is deliberately discrete and
parameter-light: it takes a structure and nothing else — no alignment, no
language model, no training — which is what allows it to be computed before
any sequence is modified and to transfer across families.

Two scoring profiles are provided. `PUBLISHED` reproduces the scores used in
the Ribbeck & Kiran dataset exactly, including the two places where the
original implementation differs from the description in that paper's Methods:
proline is counted in the hydrophobic set, and histidine in the aromatic set.
`METHODS_AS_WRITTEN` implements the Methods text literally. They disagree only
at proline and histidine positions. The default is `PUBLISHED`, so that
numbers computed with this package match the deposited data; pass the other
profile to see what the paper's prose describes.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from typing import Dict, FrozenSet, List, Optional, Sequence, Tuple

try:  # pragma: no cover - exercised only when Biopython is absent
    from Bio.PDB import MMCIFParser, NeighborSearch, PDBParser
    from Bio.PDB.SASA import ShrakeRupley
    from Bio.PDB.Polypeptide import is_aa
    from Bio.Data import IUPACData
except ImportError as exc:  # pragma: no cover
    raise ImportError(
        "efs requires Biopython; install it with `pip install biopython`"
    ) from exc


_THREE_TO_ONE = {k.upper(): v for k, v in IUPACData.protein_letters_3to1.items()}
# Modified residues that occupy a standard position and should be scored as
# their parent amino acid. Selenomethionine in particular is common in
# crystal structures and dropping it would punch holes in the sequence.
_THREE_TO_ONE.update({"MSE": "M", "SEC": "C", "PYL": "K"})

HYDROPHOBIC_PUBLISHED = frozenset("AVILMFWYP")
HYDROPHOBIC_METHODS = frozenset("LIVFMWYA")
AROMATIC_PUBLISHED = frozenset("FYWH")
AROMATIC_METHODS = frozenset("FYW")


@dataclass(frozen=True)
class ScoringProfile:
    """Parameters of the Essential Folding Score.

    All thresholds are exposed rather than hard-coded so the score can be
    recalibrated for a different structural regime (membrane proteins, say,
    where burial means something else) without editing this module.
    """

    name: str
    deep_burial_sasa: float = 20.0
    deep_burial_points: int = 3
    moderate_burial_sasa: float = 40.0
    moderate_burial_points: int = 1
    high_packing_contacts: int = 8
    high_packing_points: int = 2
    moderate_packing_contacts: int = 5
    moderate_packing_points: int = 1
    hydrophobic_set: FrozenSet[str] = HYDROPHOBIC_PUBLISHED
    hydrophobic_sasa: float = 30.0
    hydrophobic_points: int = 2
    cysteine_points: int = 2
    aromatic_set: FrozenSet[str] = AROMATIC_PUBLISHED
    aromatic_sasa: float = 30.0
    aromatic_points: int = 1
    secondary_structure_points: int = 1
    scale: int = 4
    contact_cutoff: float = 8.0
    sasa_probe_radius: float = 1.4

    def max_raw(self) -> int:
        """Largest raw score any residue can attain under this profile.

        Computed rather than asserted, because the bonuses are not all
        mutually compatible: cysteine is neither in the hydrophobic nor the
        aromatic set under either profile, so the maximum is the better of
        the cysteine route and the hydrophobic-aromatic route.
        """
        base = self.deep_burial_points + self.high_packing_points
        chem_routes = [self.cysteine_points]
        both_sets = self.hydrophobic_set & self.aromatic_set
        if both_sets:
            chem_routes.append(self.hydrophobic_points + self.aromatic_points)
        else:
            chem_routes.append(max(self.hydrophobic_points, self.aromatic_points))
        return base + max(chem_routes) + self.secondary_structure_points

    def max_score(self) -> int:
        return self.max_raw() * self.scale


PUBLISHED = ScoringProfile(name="published")
METHODS_AS_WRITTEN = ScoringProfile(
    name="methods_as_written",
    hydrophobic_set=HYDROPHOBIC_METHODS,
    aromatic_set=AROMATIC_METHODS,
)

PROFILES: Dict[str, ScoringProfile] = {
    PUBLISHED.name: PUBLISHED,
    METHODS_AS_WRITTEN.name: METHODS_AS_WRITTEN,
}

DEFAULT_CORE_THRESHOLD = 20


@dataclass(frozen=True)
class ResidueScore:
    """One residue's score and the quantities it was computed from."""

    index: int           # 1-based position along the scored chain
    pdb_resnum: int
    chain_id: str
    aa: str
    sasa: float
    contacts: int
    secondary_structure: Optional[str]
    raw_score: int
    score: int

    def is_core(self, threshold: int = DEFAULT_CORE_THRESHOLD) -> bool:
        return self.score >= threshold


@dataclass(frozen=True)
class ScoredChain:
    """Scores for every standard residue of one chain, in sequence order."""

    source: str
    chain_id: str
    profile: ScoringProfile
    residues: Tuple[ResidueScore, ...] = field(default_factory=tuple)

    @property
    def sequence(self) -> str:
        return "".join(r.aa for r in self.residues)

    def __len__(self) -> int:
        return len(self.residues)

    def scores(self) -> List[int]:
        return [r.score for r in self.residues]

    def core_positions(self, threshold: int = DEFAULT_CORE_THRESHOLD) -> FrozenSet[int]:
        """1-based indices whose score meets the threshold.

        A threshold above the maximum attainable score raises rather than
        returning an empty set. Selecting nothing is indistinguishable, in the
        output, from selecting nothing on purpose: a position-specific
        reduction with an empty core is byte-identical to a uniform one, so an
        unreachable threshold silently relabels a blind arm as a
        core-preserving one. That is not hypothetical — the published C40
        condition is exactly this, the score being bounded at 32.
        """
        ceiling = self.profile.max_score()
        if threshold > ceiling:
            raise ValueError(
                f"threshold {threshold} exceeds the maximum attainable score "
                f"{ceiling} for profile {self.profile.name!r}; no position can "
                "ever qualify, so this condition would be a uniform reduction "
                "wearing a core-preserving label"
            )
        return frozenset(r.index for r in self.residues if r.score >= threshold)

    def core_fraction(self, threshold: int = DEFAULT_CORE_THRESHOLD) -> float:
        if not self.residues:
            return 0.0
        return len(self.core_positions(threshold)) / len(self.residues)


class StructureError(ValueError):
    """Raised when a structure cannot be scored."""


def _parse(path: str, quiet: bool = True):
    lowered = path.lower()
    if lowered.endswith((".cif", ".mmcif")):
        parser = MMCIFParser(QUIET=quiet)
    else:
        parser = PDBParser(QUIET=quiet)
    return parser.get_structure("structure", path)


def _three_to_one(resname: str) -> str:
    return _THREE_TO_ONE.get(resname.strip().upper(), "X")


def _select_chain(model, chain_id: Optional[str]):
    if chain_id is not None:
        if chain_id not in model:
            available = ", ".join(c.id for c in model)
            raise StructureError(
                f"chain {chain_id!r} not present; chains available: {available}"
            )
        return model[chain_id]
    for chain in model:
        if any(is_aa(res, standard=True) for res in chain):
            return chain
    raise StructureError("no chain containing standard amino acids was found")


def _compute_sasa(structure, probe_radius: float) -> Dict[Tuple[str, tuple], float]:
    """Per-residue SASA in the context of the whole structure.

    SASA is computed on the full structure, not the isolated chain, so that
    an interface residue buried by a partner chain is recognised as buried.
    """
    ShrakeRupley(probe_radius=probe_radius).compute(structure, level="R")
    out: Dict[Tuple[str, tuple], float] = {}
    for model in structure:
        for chain in model:
            for res in chain:
                if is_aa(res, standard=True) or res.get_resname().upper() in _THREE_TO_ONE:
                    out[(chain.id, res.id)] = float(getattr(res, "sasa", 0.0))
        break  # first model only
    return out


def _compute_contacts(residues, cutoff: float) -> Dict[Tuple[str, tuple], int]:
    """Count of other CA atoms within `cutoff` angstroms of each CA."""
    keyed = [
        ((res.get_parent().id, res.id), res["CA"])
        for res in residues
        if "CA" in res
    ]
    counts = {key: 0 for key, _ in keyed}
    if not keyed:
        return counts
    search = NeighborSearch([atom for _, atom in keyed])
    for key, atom in keyed:
        neighbours = search.search(atom.coord, cutoff)
        counts[key] = sum(
            1
            for n in neighbours
            if (n.get_parent().get_parent().id, n.get_parent().id) != key
        )
    return counts


def _dssp_map(model, path: str) -> Dict[Tuple[str, tuple], str]:
    """Secondary structure per residue, or an empty map if DSSP is unusable.

    DSSP is an external binary and frequently absent. Its contribution to the
    score is optional and small; failing to find it is not an error.
    """
    try:
        from Bio.PDB import DSSP
    except ImportError:
        return {}
    try:
        dssp = DSSP(model, path)
    except Exception:
        return {}
    return {
        key: (row[2] if row[2] else "-")
        for key, row in dssp.property_dict.items()
    }


def raw_score(
    aa: str,
    sasa: float,
    contacts: int,
    secondary_structure: Optional[str] = None,
    profile: ScoringProfile = PUBLISHED,
) -> int:
    """Unscaled Essential Folding Score for a single residue.

    Pure function of the four per-residue quantities, so it can be tested and
    reused without a structure in hand.
    """
    raw = 0
    if sasa < profile.deep_burial_sasa:
        raw += profile.deep_burial_points
    elif sasa < profile.moderate_burial_sasa:
        raw += profile.moderate_burial_points

    if contacts >= profile.high_packing_contacts:
        raw += profile.high_packing_points
    elif contacts >= profile.moderate_packing_contacts:
        raw += profile.moderate_packing_points

    if secondary_structure in ("H", "E"):
        raw += profile.secondary_structure_points

    if aa in profile.hydrophobic_set and sasa < profile.hydrophobic_sasa:
        raw += profile.hydrophobic_points

    if aa == "C":
        raw += profile.cysteine_points

    if aa in profile.aromatic_set and sasa < profile.aromatic_sasa:
        raw += profile.aromatic_points

    return raw


def score_structure(
    path: str,
    chain_id: Optional[str] = None,
    profile: ScoringProfile = PUBLISHED,
    use_secondary_structure: bool = False,
) -> ScoredChain:
    """Compute the Essential Folding Score for every residue of one chain.

    `use_secondary_structure` is off by default and requires the DSSP binary;
    the deposited dataset was computed without it, and turning it on shifts
    every helix and strand residue up by one raw point.
    """
    structure = _parse(path)
    try:
        model = next(iter(structure))
    except StopIteration:
        raise StructureError(f"{path}: structure contains no models")

    chain = _select_chain(model, chain_id)
    residues = [
        res
        for res in chain
        if is_aa(res, standard=True) or res.get_resname().upper() in _THREE_TO_ONE
    ]
    if not residues:
        raise StructureError(f"{path}: chain {chain.id} has no scorable residues")

    sasa = _compute_sasa(structure, profile.sasa_probe_radius)
    contacts = _compute_contacts(residues, profile.contact_cutoff)
    ss = _dssp_map(model, path) if use_secondary_structure else {}
    effective = profile if use_secondary_structure else replace(
        profile, secondary_structure_points=0
    )

    scored: List[ResidueScore] = []
    for i, res in enumerate(residues, start=1):
        key = (chain.id, res.id)
        aa = _three_to_one(res.get_resname())
        if aa == "X":
            continue
        res_sasa = sasa.get(key, 0.0)
        res_contacts = contacts.get(key, 0)
        res_ss = ss.get(key)
        raw = raw_score(aa, res_sasa, res_contacts, res_ss, effective)
        scored.append(
            ResidueScore(
                index=i,
                pdb_resnum=res.id[1],
                chain_id=chain.id,
                aa=aa,
                sasa=res_sasa,
                contacts=res_contacts,
                secondary_structure=res_ss,
                raw_score=raw,
                score=raw * effective.scale,
            )
        )

    # Re-index so that indices are contiguous after any unknown residues were
    # dropped; downstream reduction maps positions onto the sequence this
    # object reports, and a gap would silently misalign them.
    scored = [replace(r, index=i) for i, r in enumerate(scored, start=1)]

    return ScoredChain(
        source=path,
        chain_id=chain.id,
        profile=effective,
        residues=tuple(scored),
    )


def score_sequence_features(
    sequence: str,
    sasa: Sequence[float],
    contacts: Sequence[int],
    secondary_structure: Optional[Sequence[Optional[str]]] = None,
    profile: ScoringProfile = PUBLISHED,
) -> List[int]:
    """Score from precomputed per-residue features rather than a structure file.

    For callers who already have SASA and contact counts from another source
    and want the same arithmetic applied.
    """
    n = len(sequence)
    if len(sasa) != n or len(contacts) != n:
        raise ValueError(
            f"feature length mismatch: sequence {n}, sasa {len(sasa)}, "
            f"contacts {len(contacts)}"
        )
    if secondary_structure is not None and len(secondary_structure) != n:
        raise ValueError("secondary_structure length does not match sequence")
    ss = secondary_structure or [None] * n
    return [
        raw_score(aa, sasa[i], contacts[i], ss[i], profile) * profile.scale
        for i, aa in enumerate(sequence)
    ]
