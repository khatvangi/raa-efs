"""Tests for the efs package.

The scoring tests are written against the arithmetic rather than against
stored outputs, so they check that the implementation does what it claims
rather than that it reproduces a particular run. The alphabet tests exist
because the error they guard against — an alphabet whose name disagrees with
its group count — is one this project actually made and carried into print.
"""

import io
import re

import pytest

from efs import alphabets, reduce as reduce_mod, score as score_mod
from efs.alphabets import AlphabetError, STANDARD_AA


# --------------------------------------------------------------------------
# Alphabets
# --------------------------------------------------------------------------


@pytest.mark.parametrize("name", sorted(alphabets.CANONICAL))
def test_canonical_alphabet_covers_all_standard_residues(name):
    mapping = alphabets.get(name)
    assert set(mapping) == set(STANDARD_AA)


@pytest.mark.parametrize("name", sorted(alphabets.CANONICAL))
def test_canonical_alphabet_name_matches_group_count(name):
    # RAA5 -> 5, RAA11 -> 11, RAA10_AG_MERGED -> 10, and so on. This is the
    # check the project lacked: it is what makes an eleven-group mapping
    # called RAA10 impossible to load.
    digits = "".join(c for c in name.split("_")[0] if c.isdigit())
    assert digits, f"{name} carries no size in its name"
    assert alphabets.size(name) == int(digits)


def test_the_as_printed_and_as_run_largest_alphabets_differ_in_the_basic_group():
    # Supplementary Table S1 prints K, R and H as one group; the code that
    # generated the deposited sequences keeps them separate. The sequences
    # settle it in favour of the nineteen-group version.
    assert alphabets.RAA19["R"] != alphabets.RAA19["K"]
    assert alphabets.RAA17_AS_PRINTED["R"] == alphabets.RAA17_AS_PRINTED["K"]
    assert alphabets.size("RAA19") == 19
    assert alphabets.size("RAA17_AS_PRINTED") == 17


@pytest.mark.parametrize("name", sorted(alphabets.LEGACY))
def test_legacy_alphabet_name_matches_group_count(name):
    digits = "".join(c for c in name if c.isdigit())
    assert alphabets.size(name) == int(digits)


def test_a_misnamed_alphabet_is_rejected_at_construction():
    eleven_groups = dict(alphabets.RAA11)
    with pytest.raises(AlphabetError, match="advertised size 10"):
        alphabets._validate("pretend_RAA10", eleven_groups, 10)


def test_an_incomplete_alphabet_is_rejected():
    incomplete = {aa: aa for aa in "ACDEFGHIKLMNPQRSTVW"}  # no tyrosine
    with pytest.raises(AlphabetError, match="omits Y"):
        alphabets._validate("incomplete", incomplete, 19)


def test_historical_names_resolve_to_their_true_sized_alphabets():
    assert alphabets.get("RAA10") is alphabets.RAA11
    assert alphabets.get("RAA13") is alphabets.RAA19
    assert alphabets.size("RAA10") == 11
    assert alphabets.size("RAA13") == 19


def test_unknown_alphabet_names_raise():
    with pytest.raises(KeyError):
        alphabets.get("RAA42")


def test_hierarchical_refinement_holds_from_raa7_up():
    # Every pair of residues sharing a group in the finer alphabet must share
    # a group in the coarser one. This is the property that lets a change in
    # folding outcome between adjacent alphabets be attributed to the single
    # distinction introduced at that step.
    chain = ["RAA7", "RAA8", "RAA9", "RAA11", "RAA19"]
    for coarse_name, fine_name in zip(chain, chain[1:]):
        coarse = alphabets.get(coarse_name)
        fine = alphabets.get(fine_name)
        for a in sorted(STANDARD_AA):
            for b in sorted(STANDARD_AA):
                if fine[a] == fine[b]:
                    assert coarse[a] == coarse[b], (
                        f"{a} and {b} share a group in {fine_name} "
                        f"but not in {coarse_name}"
                    )


def test_raa5_is_not_a_parent_of_raa7():
    # Documented re-assignment: alanine sits with the large hydrophobics in
    # RAA5 and with glycine in RAA7, so the refinement chain does not extend
    # down to RAA5. If this ever starts passing as a refinement, the claim in
    # the Methods about RAA5 being a separate partition has silently changed.
    raa5, raa7 = alphabets.RAA5, alphabets.RAA7
    violations = [
        (a, b)
        for a in STANDARD_AA
        for b in STANDARD_AA
        if raa7[a] == raa7[b] and raa5[a] != raa5[b]
    ]
    assert violations


def test_the_two_ten_letter_alphabets_sit_between_raa9_and_raa11():
    nine, eleven = alphabets.RAA9, alphabets.RAA11
    for name in ("RAA10_AG_MERGED", "RAA10_POLAR_MERGED"):
        ten = alphabets.get(name)
        for a in STANDARD_AA:
            for b in STANDARD_AA:
                if eleven[a] == eleven[b]:
                    assert ten[a] == ten[b], f"{name} is not coarser than RAA11"
                if ten[a] == ten[b]:
                    assert nine[a] == nine[b], f"{name} is not finer than RAA9"


def test_the_two_ten_letter_alphabets_differ_from_each_other():
    ag = alphabets.RAA10_AG_MERGED
    polar = alphabets.RAA10_POLAR_MERGED
    assert ag["A"] == ag["G"] and polar["A"] != polar["G"]
    assert polar["S"] == polar["N"] and ag["S"] != ag["N"]


def test_canonical_and_legacy_five_letter_alphabets_emit_disjoint_letters():
    # This is what makes provenance identifiable after the fact.
    canonical = alphabets.output_letters("RAA5")
    legacy = alphabets.output_letters("LEGACY_5")
    assert canonical == set("LGSDK")
    assert legacy == set("LAGEQ")
    assert canonical & legacy == {"L", "G"}


def test_legacy_seven_letter_alphabet_maps_cysteine_to_aspartate():
    # The known bug, pinned so nobody "fixes" the legacy map and thereby
    # makes old results irreproducible.
    assert alphabets.LEGACY_7["C"] == "D"
    assert alphabets.RAA7["C"] != "D"


def test_identify_distinguishes_the_two_generations():
    native = "ACDEFGHIKLMNPQRSTVWY"
    canonical = "".join(alphabets.RAA5[aa] for aa in native)
    legacy = "".join(alphabets.LEGACY_5[aa] for aa in native)
    assert "RAA5" in alphabets.identify(canonical)
    assert "LEGACY_5" not in alphabets.identify(canonical)
    assert "LEGACY_5" in alphabets.identify(legacy)
    assert "RAA5" not in alphabets.identify(legacy)


# --------------------------------------------------------------------------
# Scoring arithmetic
# --------------------------------------------------------------------------


def test_a_buried_packed_leucine_scores_at_the_top_of_the_hydrophobic_route():
    raw = score_mod.raw_score("L", sasa=5.0, contacts=12)
    assert raw == 3 + 2 + 2  # deep burial, high packing, buried hydrophobic
    assert raw * score_mod.PUBLISHED.scale == 28


def test_an_exposed_loop_residue_scores_zero():
    assert score_mod.raw_score("S", sasa=120.0, contacts=2) == 0


def test_burial_bands_are_exclusive():
    assert score_mod.raw_score("S", sasa=19.9, contacts=0) == 3
    assert score_mod.raw_score("S", sasa=20.0, contacts=0) == 1
    assert score_mod.raw_score("S", sasa=39.9, contacts=0) == 1
    assert score_mod.raw_score("S", sasa=40.0, contacts=0) == 0


def test_packing_bands_are_exclusive():
    assert score_mod.raw_score("S", sasa=100.0, contacts=8) == 2
    assert score_mod.raw_score("S", sasa=100.0, contacts=7) == 1
    assert score_mod.raw_score("S", sasa=100.0, contacts=5) == 1
    assert score_mod.raw_score("S", sasa=100.0, contacts=4) == 0


def test_cysteine_scores_even_when_fully_exposed():
    # Deliberate: a solvent-exposed cysteine can still be holding a disulfide.
    assert score_mod.raw_score("C", sasa=150.0, contacts=0) == 2


def test_buried_aromatic_takes_both_the_hydrophobic_and_aromatic_bonus():
    raw = score_mod.raw_score("F", sasa=10.0, contacts=10)
    assert raw == 3 + 2 + 2 + 1


def test_profiles_differ_only_at_proline_and_histidine():
    disagreements = {
        aa
        for aa in STANDARD_AA
        if score_mod.raw_score("".join(aa), 5.0, 10, None, score_mod.PUBLISHED)
        != score_mod.raw_score("".join(aa), 5.0, 10, None, score_mod.METHODS_AS_WRITTEN)
    }
    assert disagreements == {"P", "H"}


def test_secondary_structure_contributes_only_for_helix_and_strand():
    base = score_mod.raw_score("S", sasa=100.0, contacts=0)
    assert score_mod.raw_score("S", 100.0, 0, "H") == base + 1
    assert score_mod.raw_score("S", 100.0, 0, "E") == base + 1
    assert score_mod.raw_score("S", 100.0, 0, "T") == base
    assert score_mod.raw_score("S", 100.0, 0, None) == base


def test_max_raw_is_computed_not_asserted():
    profile = score_mod.PUBLISHED
    observed = max(
        score_mod.raw_score(aa, 0.0, 99, "H", profile) for aa in STANDARD_AA
    )
    assert observed == profile.max_raw()


def test_score_sequence_features_matches_per_residue_scoring():
    seq = "LLSC"
    sasa = [5.0, 100.0, 10.0, 150.0]
    contacts = [12, 1, 6, 0]
    got = score_mod.score_sequence_features(seq, sasa, contacts)
    want = [
        score_mod.raw_score(aa, sasa[i], contacts[i]) * score_mod.PUBLISHED.scale
        for i, aa in enumerate(seq)
    ]
    assert got == want


def test_score_sequence_features_rejects_mismatched_feature_lengths():
    with pytest.raises(ValueError, match="length mismatch"):
        score_mod.score_sequence_features("LLL", [1.0, 2.0], [1, 2, 3])


# --------------------------------------------------------------------------
# Reduction
# --------------------------------------------------------------------------


NATIVE = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ"


def test_uniform_reduction_emits_only_that_alphabets_letters():
    result = reduce_mod.reduce_uniform(NATIVE, "RAA5")
    assert set(result.sequence) <= alphabets.output_letters("RAA5")
    assert result.preserved == frozenset()
    assert len(result.sequence) == len(NATIVE)


def test_position_specific_reduction_keeps_exactly_the_named_positions():
    preserved = {1, 5, 9}
    result = reduce_mod.reduce_position_specific(NATIVE, "RAA5", preserved)
    for i, (native_aa, reduced_aa) in enumerate(zip(NATIVE, result.sequence), start=1):
        if i in preserved:
            assert reduced_aa == native_aa
        else:
            assert reduced_aa == alphabets.RAA5[native_aa]


def test_preserving_every_position_returns_the_native_sequence():
    result = reduce_mod.reduce_position_specific(
        NATIVE, "RAA5", range(1, len(NATIVE) + 1)
    )
    assert result.sequence == NATIVE
    assert result.percent_identity_to_native == 100.0


def test_percent_identity_exceeds_preserved_fraction():
    # Reduced positions whose native residue is its own group representative
    # are unchanged too, so identity is always at least the preserved
    # fraction. Conflating the two is how a proximity effect gets mistaken
    # for a placement effect.
    result = reduce_mod.reduce_position_specific(NATIVE, "RAA5", {1, 2, 3})
    assert result.percent_identity_to_native >= 100.0 * result.preserved_fraction


def test_out_of_range_positions_are_rejected():
    with pytest.raises(ValueError, match="out of range"):
        reduce_mod.reduce_position_specific(NATIVE, "RAA5", {0})
    with pytest.raises(ValueError, match="out of range"):
        reduce_mod.reduce_position_specific(NATIVE, "RAA5", {len(NATIVE) + 1})


def test_non_standard_residues_are_rejected_rather_than_passed_through():
    with pytest.raises(ValueError, match="non-standard"):
        reduce_mod.reduce_uniform("MKTXYI", "RAA5")


def test_random_control_is_reproducible_and_seed_dependent():
    a = reduce_mod.random_matched_control(NATIVE, "RAA5", n_preserved=8, seed=1)
    b = reduce_mod.random_matched_control(NATIVE, "RAA5", n_preserved=8, seed=1)
    c = reduce_mod.random_matched_control(NATIVE, "RAA5", n_preserved=8, seed=2)
    assert a.sequence == b.sequence
    assert a.sequence != c.sequence


def test_random_control_preserves_exactly_the_requested_count():
    for n in (0, 1, 7, len(NATIVE)):
        result = reduce_mod.random_matched_control(NATIVE, "RAA5", n, seed=3)
        assert len(result.preserved) == n


def test_random_control_rejects_an_impossible_count():
    with pytest.raises(ValueError):
        reduce_mod.random_matched_control(NATIVE, "RAA5", len(NATIVE) + 1, seed=1)


def test_identify_alphabet_from_a_position_specific_reduction():
    preserved = {2, 4, 6, 8}
    result = reduce_mod.reduce_position_specific(NATIVE, "RAA5", preserved)
    candidates = reduce_mod.identify_alphabet_from_reduction(
        result.sequence, NATIVE, preserved
    )
    assert "RAA5" in candidates
    assert "LEGACY_5" not in candidates


def test_identify_without_knowing_the_preserved_positions():
    result = reduce_mod.reduce_position_specific(NATIVE, "RAA5", {2, 4})
    candidates = reduce_mod.identify_alphabet_from_reduction(result.sequence, NATIVE)
    assert "RAA5" in candidates


def test_fasta_output_wraps_and_labels():
    result = reduce_mod.reduce_uniform(NATIVE, "RAA5")
    text = result.fasta("P12345", width=10)
    lines = text.strip().split("\n")
    assert lines[0] == ">P12345_blind_RAA5"
    assert all(len(line) <= 10 for line in lines[1:])
    assert "".join(lines[1:]) == result.sequence


# --------------------------------------------------------------------------
# Structure-level behaviour, exercised on a synthetic PDB
# --------------------------------------------------------------------------


def _synthetic_pdb(residues, spacing=3.8):
    """A straight CA-only chain; enough for parsing, SASA and contact code."""
    lines = []
    three = {
        "A": "ALA", "C": "CYS", "D": "ASP", "E": "GLU", "F": "PHE",
        "G": "GLY", "H": "HIS", "I": "ILE", "K": "LYS", "L": "LEU",
        "M": "MET", "N": "ASN", "P": "PRO", "Q": "GLN", "R": "ARG",
        "S": "SER", "T": "THR", "V": "VAL", "W": "TRP", "Y": "TYR",
    }
    for i, aa in enumerate(residues, start=1):
        x = i * spacing
        lines.append(
            f"ATOM  {i:>5}  CA  {three[aa]} A{i:>4}    "
            f"{x:>8.3f}{0.0:>8.3f}{0.0:>8.3f}  1.00  0.00           C"
        )
    lines.append("END")
    return "\n".join(lines) + "\n"


@pytest.fixture
def straight_chain(tmp_path):
    path = tmp_path / "chain.pdb"
    path.write_text(_synthetic_pdb("LLCFSSGGKKDD"))
    return str(path)


def test_scoring_a_structure_returns_one_entry_per_residue(straight_chain):
    scored = score_mod.score_structure(straight_chain)
    assert len(scored) == 12
    assert scored.sequence == "LLCFSSGGKKDD"
    assert [r.index for r in scored.residues] == list(range(1, 13))


def test_scored_indices_are_contiguous_and_one_based(straight_chain):
    scored = score_mod.score_structure(straight_chain)
    assert scored.residues[0].index == 1
    assert scored.residues[-1].index == len(scored)


def test_core_positions_respect_the_threshold(straight_chain):
    scored = score_mod.score_structure(straight_chain)
    assert len(scored.core_positions(threshold=0)) == len(scored)
    assert scored.core_positions(threshold=scored.profile.max_score()) is not None


def test_an_unreachable_threshold_raises_rather_than_selecting_nothing(straight_chain):
    # The published C40 condition: the score is bounded at 32, so a cutoff of
    # 40 selects no position and the resulting "core-preserving" sequences are
    # byte-identical to the uniform reduction. An empty return would hide it.
    scored = score_mod.score_structure(straight_chain)
    ceiling = scored.profile.max_score()
    assert ceiling == 32
    with pytest.raises(ValueError, match="exceeds the maximum attainable score"):
        scored.core_positions(threshold=40)
    with pytest.raises(ValueError):
        scored.core_fraction(threshold=ceiling + 1)


def test_reduce_by_efs_refuses_an_unreachable_threshold(straight_chain):
    scored = score_mod.score_structure(straight_chain)
    with pytest.raises(ValueError, match="exceeds the maximum attainable score"):
        reduce_mod.reduce_by_efs(scored, "RAA5", threshold=40)


def test_requesting_a_missing_chain_names_the_chains_that_exist(straight_chain):
    with pytest.raises(score_mod.StructureError, match="chains available: A"):
        score_mod.score_structure(straight_chain, chain_id="Z")


def test_reduce_by_efs_preserves_exactly_the_core(straight_chain):
    scored = score_mod.score_structure(straight_chain)
    threshold = 8
    result = reduce_mod.reduce_by_efs(scored, "RAA5", threshold=threshold)
    assert result.preserved == scored.core_positions(threshold)
    assert result.label == "EFS_C8_RAA5"


def test_random_control_matching_is_size_matched_to_the_core(straight_chain):
    scored = score_mod.score_structure(straight_chain)
    threshold = 8
    efs = reduce_mod.reduce_by_efs(scored, "RAA5", threshold=threshold)
    control = reduce_mod.random_control_matching(
        scored, "RAA5", seed=1, threshold=threshold
    )
    assert len(control.preserved) == len(efs.preserved)
    assert control.preserved_fraction == efs.preserved_fraction


def test_summarise_reports_one_row_per_reduction(straight_chain):
    scored = score_mod.score_structure(straight_chain)
    rows = reduce_mod.summarise(
        {
            "blind": reduce_mod.reduce_uniform(scored.sequence, "RAA5"),
            "efs": reduce_mod.reduce_by_efs(scored, "RAA5", threshold=8),
        }
    )
    assert {row["id"] for row in rows} == {"blind", "efs"}
    assert all(row["length"] == len(scored) for row in rows)


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------


def _run_cli(argv):
    from efs.cli import build_parser

    out, err = io.StringIO(), io.StringIO()
    args = build_parser().parse_args(argv)
    code = args.func(args, out, err)
    return code, out.getvalue(), err.getvalue()


def test_cli_alphabets_reports_true_group_counts():
    code, out, _ = _run_cli(["alphabets"])
    assert code == 0
    collapsed = re.sub(r"[ \t]+", " ", out)
    assert "RAA11 11 groups" in collapsed
    assert "RAA19 19 groups" in collapsed
    assert "RAA17_AS_PRINTED 17 groups" in collapsed
    assert "RAA10 -> RAA11" in collapsed
    assert "RAA13 -> RAA19" in collapsed


def test_cli_score_emits_a_header_and_one_row_per_residue(straight_chain):
    code, out, _ = _run_cli(["score", straight_chain])
    assert code == 0
    lines = out.strip().split("\n")
    assert lines[0].split("\t")[0] == "protein_id"
    assert len(lines) == 13


def test_cli_reduce_writes_fasta(straight_chain):
    code, out, _ = _run_cli(["reduce", straight_chain, "--alphabet", "RAA5"])
    assert code == 0
    assert out.startswith(">chain_EFS_C20_RAA5")


def test_cli_control_writes_one_record_per_seed(straight_chain):
    code, out, _ = _run_cli(["control", straight_chain, "--seed", "1", "2", "3"])
    assert code == 0
    assert out.count(">") == 3


def test_cli_skips_unreadable_structures_and_reports_them(tmp_path):
    bad = tmp_path / "empty.pdb"
    bad.write_text("")
    code, _, err = _run_cli(["score", str(bad)])
    assert code == 1
    assert "skipped" in err
