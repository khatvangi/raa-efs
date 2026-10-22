# efs — Essential Folding Score

A per-residue score, computed from a structure alone, that identifies the
positions where a protein sequence cannot be written in a reduced amino acid
alphabet without losing the fold. It needs no alignment, no language model
and no training, so it can be computed before any sequence is modified and
applied to a family it has never seen.

Reduced-alphabet work has, for thirty years, applied one mapping uniformly to
every position and asked how few letters a protein can be written in. The
answer turns out to depend less on how far the alphabet is compressed than on
where the compression falls: five letters suffice almost everywhere, and fail
at roughly a quarter of positions. This package computes that quarter and
performs the reduction.

## Install

    pip install -e .   # installs raa-efs, import raa_efs

Requires Python 3.9 or later, Biopython and NumPy. DSSP is optional and off
by default.

## Use

    from raa_efs import score_structure, reduce_by_efs, random_control_matching

    scored = score_structure("kinase.pdb")
    scored.core_fraction(threshold=20)        # ~0.26 for a typical kinase

    variant = reduce_by_efs(scored, "RAA5", threshold=20)
    variant.sequence
    variant.percent_identity_to_native

    control = random_control_matching(scored, "RAA5", seed=1)

`random_control_matching` preserves the same *number* of positions as the EFS
core for that protein but chooses them at random. Report it alongside any
position-specific result. Preserving a quarter of positions necessarily moves
a variant closer to its native sequence, which on its own can improve a
structure prediction; the rescue is attributable to *which* positions were
protected only if the size-matched random draw does worse. `percent_identity_
to_native` quantifies the same confound directly, and is always at least the
preserved fraction because a reduced position whose native residue is its own
group representative is unchanged too.

From the command line:

    raa-efs score     structure.pdb --summary          # per-residue TSV
    raa-efs reduce    structure.pdb --alphabet RAA5    # EFS-guided FASTA
    raa-efs reduce    structure.pdb --blind            # uniform reduction
    raa-efs control   structure.pdb --seed 1 2 3       # matched random controls
    raa-efs alphabets                                  # alphabets and true sizes
    raa-efs identify  reduced.fasta native.fasta       # which alphabet made this

## The score

Each residue accumulates points for burial (SASA below 20 Å² scores 3, below
40 Å² scores 1), packing (8 or more Cα neighbours within 8 Å scores 2, 5 to 7
scores 1), a buried hydrophobic bonus of 2, a cysteine bonus of 2 awarded
regardless of burial because an exposed cysteine may still hold a disulfide,
and a buried aromatic bonus of 1. The sum is multiplied by 4. Positions
scoring 20 or above are the core by default.

Every threshold is a field on `ScoringProfile` rather than a constant in the
code, so the score can be recalibrated — for membrane proteins, where burial
means something else — without editing the module.

The maximum attainable score is computed from the profile rather than
asserted, and it is 32 under the published parameters: the cysteine bonus and
the hydrophobic-plus-aromatic route are mutually exclusive, so the ceiling is
lower than adding every bonus suggests. A core threshold above that ceiling
raises instead of returning an empty set, because an empty core produces a
sequence byte-identical to a uniform reduction, and a uniform reduction
labelled as core-preserving is invisible in the output. The published C40
condition is exactly that: a cutoff of 40 against a ceiling of 32, selecting
nothing in all 627 proteins.

Two profiles ship. `PUBLISHED` is the default and reproduces the deposited
dataset exactly. `METHODS_AS_WRITTEN` implements the Methods text of the
accompanying manuscript literally. They differ in two places: the original
implementation counts proline in the hydrophobic set and histidine in the
aromatic set, and the Methods text does not. The profiles therefore disagree
only at buried proline and buried histidine positions. `tools/validate_
against_published.py` re-runs the arithmetic against the stored SASA and
contact columns of the deposited per-residue files and reports which profile
reproduces them.

## Alphabets

`raa-efs alphabets` prints every alphabet with the number of distinct groups it
actually produces, which is checked at import time against the number its
name claims. That check exists because this project got it wrong: the
mapping recorded as "RAA10" has eleven groups and the one recorded as "RAA13"
has nineteen. They are available here under their true names, `RAA11` and
`RAA19`, with the old names kept as aliases so older scripts and filenames
still resolve.

The largest alphabet is worse than misnamed: the code that generated the
deposited sequences keeps lysine, arginine and histidine in separate groups,
giving nineteen, while Supplementary Table S1 prints them merged, giving
seventeen. The sequences settle it in favour of nineteen. The printed variant
is kept as `RAA17_AS_PRINTED` so the discrepancy is expressible in code; no
deposited sequence was generated with it.

`RAA7` through `RAA19` form a strict hierarchical refinement — every group at
size *n* is contained in a group at size *n*−1 — which is what allows a
difference in folding outcome between adjacent alphabets to be attributed to
the single distinction introduced at that step. The test suite checks this
property rather than trusting the comment. `RAA5` is a separate coarser
partition and is *not* the parent of `RAA7`: alanine sits with the large
hydrophobics in one and with glycine in the other.

`RAA10_AG_MERGED` and `RAA10_POLAR_MERGED` are the two genuine ten-letter
alphabets lying strictly between `RAA9` and `RAA11`. The step from eleven
letters to nine merges two things at once, alanine with glycine and the
hydroxyl class with the amide class; these two alphabets each merge one, so
the loss can be apportioned between them.

A second generation of alphabets exists in this project's history under the
same names and with different partitions — alanine with threonine and
histidine, cysteine inside the hydrophobic class, the basic residues merged
with the polars, and, at size 7, cysteine mapped to aspartate. They are kept
as `LEGACY_5` and `LEGACY_7`, marked deprecated, so that results generated
before the canonical set existed can be reproduced and, more importantly,
identified. The canonical and legacy five-letter alphabets emit disjoint
letter sets (`LGSDK` against `LAGEQ`), so a reduced sequence of unknown
provenance can be traced to the mapping that produced it:

    raa-efs identify old_results.fasta natives.fasta

Do not generate new results with the legacy mappings.

## Tests

    pip install -e ".[test]"
    pytest

The suite checks the scoring arithmetic band by band, the hierarchical
refinement property, the per-protein size matching of the random control, and
that a misnamed alphabet is rejected at construction rather than propagating
into a figure.

## Citing

Ribbeck, H. & Kiran, B. Position-specific alphabet reduction: a
structure-based score identifies where the protein code cannot be compressed.
Data at https://doi.org/10.5281/zenodo.21560977.

## Earlier contents of this repository

Until October 2026 this repository held the v5-era analysis scripts, result
tables and a manuscript draft that was withdrawn. That state is preserved at
the tag `v5-analysis` for provenance. Its README claimed that
`code/CANONICAL_ALPHABETS.py` was the single source of truth for all six
alphabets and that every published number reproduced from the deposited
tables; neither was true, because the legacy mappings that generated several
of those values were never deposited. Those mappings are available here as
`LEGACY_5` and `LEGACY_7`, and `raa-efs identify` traces a sequence of unknown
provenance to the mapping that produced it.
