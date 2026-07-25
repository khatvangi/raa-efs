#!/usr/bin/env python3
"""
CANONICAL_ALPHABETS.py - Single source of truth for RAA alphabets

All reduction scripts MUST import from this file.
Do NOT define alphabets elsewhere.

Created: 2026-01-10
Updated: 2026-04-18 - Added RAA8 and RAA9 as hierarchical interpolants
                     between RAA7 and RAA10.
"""

# standard amino acids
STANDARD_AAS = set("ACDEFGHIKLMNPQRSTVWY")


# =============================================================================
# DESIGN NOTE: The progression RAA7 -> RAA8 -> RAA9 -> RAA10 -> RAA13 forms a
# STRICT hierarchical refinement. Every group at size n is a subset of its
# parent group at size n-1. Any difference in folding outcome between RAAn
# and RAA(n+1) in that chain is attributable to exactly the split(s)
# introduced at that step.
#
# RAA5 is a DIFFERENT PARTITION -- not a strict parent of RAA7. Specifically,
# A moves from {LIVMFYWA} (hydrophobic in RAA5) into a merged {A,G} small
# group in RAA7, while G moves out of {G,P}. The RAA5<->RAA7 transition is a
# re-assignment, not a pure split/merge. Interpret RAA5 as a separate very
# coarse baseline rather than as the parent of the hierarchical chain.
#
# Splits introduced at each step of the strict chain:
#   RAA7 -> RAA8:  C separated from polar group {STNQC} -> {STNQ} + {C}.
#                  Rationale: cysteine's disulfide chemistry is the single
#                  most distinctive feature after aromatic/aliphatic.
#   RAA8 -> RAA9:  M separated from large hydrophobic {LIVM} -> {LIV} + {M}.
#                  Rationale: methionine sulfur contributes specific
#                  methionine-aromatic interactions and distinctive packing.
#   RAA9 -> RAA10: {AG} splits into {A} and {G}; {STNQ} splits into {ST}
#                  and {NQ} (hydroxyl vs amide).
#   RAA10 -> RAA13: acidic {DE} splits to D/E; aromatic {FYW} splits to
#                   F/Y/W; various further refinements.
# =============================================================================


# =============================================================================
# RAA5 - 5 letter alphabet
# =============================================================================
# Groups: L (hydrophobic), G (conformational), S (polar), D (acidic), K (basic)

RAA5 = {
    # hydrophobic (large) -> L
    'L': 'L', 'I': 'L', 'V': 'L', 'M': 'L', 'F': 'L',
    'Y': 'L', 'W': 'L', 'A': 'L',
    # conformationally special -> G (includes P)
    'G': 'G', 'P': 'G',
    # polar (including C for simplicity at this level) -> S
    'S': 'S', 'T': 'S', 'N': 'S', 'Q': 'S', 'C': 'S',
    # acidic -> D
    'D': 'D', 'E': 'D',
    # basic -> K
    'K': 'K', 'R': 'K', 'H': 'K',
}

RAA5_GROUPS = {
    'L': ['L', 'I', 'V', 'M', 'F', 'Y', 'W', 'A'],  # hydrophobic
    'G': ['G', 'P'],                                  # conformational
    'S': ['S', 'T', 'N', 'Q', 'C'],                  # polar
    'D': ['D', 'E'],                                  # acidic
    'K': ['K', 'R', 'H'],                            # basic
}


# =============================================================================
# RAA7_FIXED - 7 letter alphabet (CORRECTED - chemically sound)
# =============================================================================
# Groups: L (large hydrophobic), F (aromatic), A (small), P (proline),
#         S (polar incl. C), D (acidic), K (basic)
# NOTE: C kept with polar at this level (separates at RAA8)
# NOTE: V stays with hydrophobic (not small like buggy version)
# NOTE: P kept separate (conformationally unique)

RAA7_FIXED = {
    # large hydrophobic -> L
    'L': 'L', 'I': 'L', 'V': 'L', 'M': 'L',
    # aromatic -> F
    'F': 'F', 'Y': 'F', 'W': 'F',
    # small -> A
    'A': 'A', 'G': 'A',
    # proline separate (conformationally unique)
    'P': 'P',
    # polar (including C) -> S
    'S': 'S', 'T': 'S', 'N': 'S', 'Q': 'S', 'C': 'S',
    # acidic -> D
    'D': 'D', 'E': 'D',
    # basic -> K
    'K': 'K', 'R': 'K', 'H': 'K',
}

RAA7_FIXED_GROUPS = {
    'L': ['L', 'I', 'V', 'M'],      # large hydrophobic
    'F': ['F', 'Y', 'W'],           # aromatic
    'A': ['A', 'G'],                # small
    'P': ['P'],                     # proline (unique)
    'S': ['S', 'T', 'N', 'Q', 'C'], # polar
    'D': ['D', 'E'],                # acidic
    'K': ['K', 'R', 'H'],           # basic
}


# =============================================================================
# RAA8 - 8 letter alphabet (NEW)
# =============================================================================
# RAA7 with cysteine separated from the polar group.
# Rationale: cysteine's disulfide bonding capability is chemically the most
# distinctive feature in the polar class. The historical RAA7 bug (C -> D)
# and the high fold-collapse rate when C was mismapped both demonstrate
# that C deserves its own group before any other polar distinction is made.
# ThermoMPNN data: C mutations mean delta-delta-G = +0.43 kcal/mol (lowest
# among all AAs), but C -> any other AA at a disulfide-bonding position
# incurs additional +0.5-1.0 kcal/mol penalty from disulfide loss.
#
# Groups: L (large hydrophobic), F (aromatic), A (small), P (proline),
#         S (polar, no C), C (cysteine), D (acidic), K (basic)

RAA8 = {
    # large hydrophobic -> L
    'L': 'L', 'I': 'L', 'V': 'L', 'M': 'L',
    # aromatic -> F
    'F': 'F', 'Y': 'F', 'W': 'F',
    # small -> A
    'A': 'A', 'G': 'A',
    # proline separate (conformationally unique)
    'P': 'P',
    # polar (C removed) -> S
    'S': 'S', 'T': 'S', 'N': 'S', 'Q': 'S',
    # cysteine separate (disulfide bonding)
    'C': 'C',
    # acidic -> D
    'D': 'D', 'E': 'D',
    # basic -> K
    'K': 'K', 'R': 'K', 'H': 'K',
}

RAA8_GROUPS = {
    'L': ['L', 'I', 'V', 'M'],   # large hydrophobic
    'F': ['F', 'Y', 'W'],        # aromatic
    'A': ['A', 'G'],             # small
    'P': ['P'],                  # proline (unique)
    'S': ['S', 'T', 'N', 'Q'],   # polar (no C)
    'C': ['C'],                  # cysteine (disulfide)
    'D': ['D', 'E'],             # acidic
    'K': ['K', 'R', 'H'],        # basic
}


# =============================================================================
# RAA9 - 9 letter alphabet (NEW)
# =============================================================================
# RAA8 with methionine separated from the aliphatic hydrophobic group.
# Rationale: methionine's sulfur atom enables specific methionine-aromatic
# interactions (Valley et al. 2012) and distinctive packing behavior that
# is not captured by L/I/V equivalence. Pulling M out at RAA9 leaves the
# RAA9 -> RAA10 refinement as exactly two remaining splits: {AG} -> {A}+{G}
# and {STNQ} -> {ST}+{NQ}.
#
# Groups: L (aliphatic), M (methionine), F (aromatic), A (small), P (proline),
#         S (polar, no C), C (cysteine), D (acidic), K (basic)

RAA9 = {
    # aliphatic hydrophobic (M removed) -> L
    'L': 'L', 'I': 'L', 'V': 'L',
    # methionine separate (sulfur, distinctive packing)
    'M': 'M',
    # aromatic -> F
    'F': 'F', 'Y': 'F', 'W': 'F',
    # small -> A
    'A': 'A', 'G': 'A',
    # proline separate
    'P': 'P',
    # polar (no C) -> S
    'S': 'S', 'T': 'S', 'N': 'S', 'Q': 'S',
    # cysteine separate
    'C': 'C',
    # acidic -> D
    'D': 'D', 'E': 'D',
    # basic -> K
    'K': 'K', 'R': 'K', 'H': 'K',
}

RAA9_GROUPS = {
    'L': ['L', 'I', 'V'],        # aliphatic (no M)
    'M': ['M'],                  # methionine (sulfur)
    'F': ['F', 'Y', 'W'],        # aromatic
    'A': ['A', 'G'],             # small
    'P': ['P'],                  # proline
    'S': ['S', 'T', 'N', 'Q'],   # polar
    'C': ['C'],                  # cysteine
    'D': ['D', 'E'],             # acidic
    'K': ['K', 'R', 'H'],        # basic
}


# =============================================================================
# RAA10 - 10 letter alphabet
# =============================================================================
# More granular groupings, preserves more chemical distinction.
# (Historical note: has 11 actual output groups; name retained for
# consistency with project records.)

RAA10 = {
    # aliphatic hydrophobic -> L
    'L': 'L', 'I': 'L', 'V': 'L',
    # methionine separate (sulfur)
    'M': 'M',
    # aromatic -> F
    'F': 'F', 'Y': 'F', 'W': 'F',
    # small -> separate
    'A': 'A', 'G': 'G',
    # proline separate
    'P': 'P',
    # hydroxyl -> S
    'S': 'S', 'T': 'S',
    # amide -> N
    'N': 'N', 'Q': 'N',
    # cysteine separate (disulfide)
    'C': 'C',
    # acidic -> D
    'D': 'D', 'E': 'D',
    # basic -> K
    'K': 'K', 'R': 'K', 'H': 'K',
}

RAA10_GROUPS = {
    'L': ['L', 'I', 'V'],     # aliphatic hydrophobic
    'M': ['M'],               # methionine
    'F': ['F', 'Y', 'W'],     # aromatic
    'A': ['A'],               # alanine
    'G': ['G'],               # glycine
    'P': ['P'],               # proline
    'S': ['S', 'T'],          # hydroxyl
    'N': ['N', 'Q'],          # amide
    'C': ['C'],               # cysteine
    'D': ['D', 'E'],          # acidic
    'K': ['K', 'R', 'H'],     # basic
}


# =============================================================================
# RAA13 - 13 letter alphabet
# =============================================================================
# Near-full, only minimal groupings

RAA13 = {
    # L and I together (very similar)
    'L': 'L', 'I': 'L',
    # everything else mostly separate
    'V': 'V',
    'M': 'M',
    'F': 'F', 'Y': 'Y', 'W': 'W',
    'A': 'A', 'G': 'G',
    'P': 'P',
    'S': 'S', 'T': 'T',
    'N': 'N', 'Q': 'Q',
    'C': 'C',
    'D': 'D', 'E': 'E',
    'K': 'K', 'R': 'R', 'H': 'H',
}

RAA13_GROUPS = {
    'L': ['L', 'I'],
    'V': ['V'],
    'M': ['M'],
    'F': ['F'],
    'Y': ['Y'],
    'W': ['W'],
    'A': ['A'],
    'G': ['G'],
    'P': ['P'],
    'S': ['S'],
    'T': ['T'],
    'N': ['N'],
    'Q': ['Q'],
    'C': ['C'],
    'D': ['D'],
    'E': ['E'],
    'K': ['K'],
    'R': ['R'],
    'H': ['H'],
}


# =============================================================================
# CHEMISTRY_INVERSION - for hypothesis testing
# =============================================================================
# Inverts core chemistry to test if fold collapses
# Hydrophobic -> Charged (destroys core packing)
# Cysteine -> Positive (breaks disulfides)
# Proline -> Glycine (loses rigidity)

CHEMISTRY_INVERSION = {
    # hydrophobic -> charged (destroys core packing)
    'L': 'D', 'I': 'E', 'V': 'K', 'M': 'R',
    'F': 'D', 'W': 'E', 'Y': 'K',
    'A': 'D',
    # cysteine -> positive (breaks disulfides)
    'C': 'K',
    # proline -> glycine (loses rigidity)
    'P': 'G',
    # polar/charged -> keep as-is (already exposed)
    'S': 'S', 'T': 'T', 'N': 'N', 'Q': 'Q',
    'D': 'D', 'E': 'E', 'K': 'K', 'R': 'R', 'H': 'H',
    'G': 'G',
}


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def reduce_sequence(sequence: str, alphabet: dict) -> str:
    """reduce a sequence using the given alphabet mapping"""
    return ''.join(alphabet.get(aa, 'X') for aa in sequence.upper())


def get_alphabet(name: str) -> dict:
    """get alphabet by name"""
    alphabets = {
        'RAA5': RAA5,
        'RAA7': RAA7_FIXED,  # always use fixed version
        'RAA7_FIXED': RAA7_FIXED,
        'RAA8': RAA8,
        'RAA9': RAA9,
        'RAA10': RAA10,
        'RAA13': RAA13,
        'INVERT': CHEMISTRY_INVERSION,
        'CHEMISTRY_INVERSION': CHEMISTRY_INVERSION,
    }
    return alphabets.get(name.upper(), None)


def reduce_at_positions(sequence: str, positions: set, alphabet: dict) -> str:
    """
    reduce only at specified positions, keep rest native

    args:
        sequence: original amino acid sequence
        positions: set of 0-indexed positions to reduce
        alphabet: reduction mapping to use

    returns:
        modified sequence
    """
    result = list(sequence.upper())
    for pos in positions:
        if 0 <= pos < len(result):
            aa = result[pos]
            result[pos] = alphabet.get(aa, aa)
    return ''.join(result)


def invert_core_chemistry(sequence: str, core_positions: set) -> str:
    """
    invert chemistry at core positions (for hypothesis testing)

    hydrophobic core -> charged = should collapse structure
    """
    return reduce_at_positions(sequence, core_positions, CHEMISTRY_INVERSION)


def validate_alphabet(alphabet: dict) -> bool:
    """check that alphabet maps all 20 standard AAs"""
    return all(aa in alphabet for aa in STANDARD_AAS)


def is_refinement(coarse_groups: dict, fine_groups: dict) -> bool:
    """
    check that fine_groups is a strict refinement of coarse_groups.
    returns True if every pair of AAs that share a group in fine_groups
    also share a group in coarse_groups.
    """
    def aa_to_label(groups):
        m = {}
        for label, members in groups.items():
            for aa in members:
                m[aa] = label
        return m

    coarse_map = aa_to_label(coarse_groups)
    fine_map = aa_to_label(fine_groups)
    aas = sorted(STANDARD_AAS)
    for i, aa1 in enumerate(aas):
        for aa2 in aas[i+1:]:
            if fine_map.get(aa1) == fine_map.get(aa2):
                if coarse_map.get(aa1) != coarse_map.get(aa2):
                    return False
    return True


# =============================================================================
# DEPRECATED - DO NOT USE
# =============================================================================
# these are the buggy mappings for reference only

RAA7_BUGGY = {
    # WARNING: chemically incorrect - DO NOT USE
    'A': 'A', 'G': 'A', 'V': 'A',      # V->A WRONG (V is hydrophobic)
    'I': 'I', 'L': 'I', 'F': 'I', 'P': 'I',  # P->I WRONG (P is structural)
    'Y': 'Y', 'M': 'Y',
    'T': 'T', 'S': 'T',
    'H': 'H', 'N': 'H', 'Q': 'H', 'W': 'H',
    'R': 'R', 'K': 'R',
    'D': 'D', 'E': 'D', 'C': 'D',      # C->D WRONG (C forms disulfides)
}


# =============================================================================
# MAIN - self-test
# =============================================================================

if __name__ == "__main__":
    print("=== CANONICAL ALPHABETS SELF-TEST ===\n")

    # test sequence
    test_seq = "ACDEFGHIKLMNPQRSTVWY"

    for name, alpha in [('RAA5', RAA5),
                        ('RAA7_FIXED', RAA7_FIXED),
                        ('RAA8', RAA8),
                        ('RAA9', RAA9),
                        ('RAA10', RAA10),
                        ('RAA13', RAA13)]:
        reduced = reduce_sequence(test_seq, alpha)
        unique = sorted(set(reduced))
        print(f"{name}: {len(unique)} letters -> {reduced}")
        print(f"       unique: {''.join(unique)}")
        assert validate_alphabet(alpha), f"{name} missing AAs!"
        print()

    # test chemistry inversion
    print("CHEMISTRY_INVERSION:")
    inverted = reduce_sequence(test_seq, CHEMISTRY_INVERSION)
    print(f"  Original: {test_seq}")
    print(f"  Inverted: {inverted}")
    print()

    # verify hierarchical refinement property (RAA7 onwards only)
    # NOTE: RAA5 is NOT a strict parent of RAA7 because A moves between
    # groups (see design note above). We assert RAA5<->RAA7 is non-refinement.
    print("Hierarchical refinement checks:")
    assert not is_refinement(RAA5_GROUPS, RAA7_FIXED_GROUPS), \
        "RAA5 unexpectedly appears to be a strict parent of RAA7"
    assert is_refinement(RAA7_FIXED_GROUPS, RAA8_GROUPS), "RAA8 not refinement of RAA7"
    assert is_refinement(RAA8_GROUPS, RAA9_GROUPS), "RAA9 not refinement of RAA8"
    assert is_refinement(RAA9_GROUPS, RAA10_GROUPS), "RAA10 not refinement of RAA9"
    assert is_refinement(RAA10_GROUPS, RAA13_GROUPS), "RAA13 not refinement of RAA10"
    print("  RAA5 is a separate partition (A/G re-assigned) - OK")
    print("  RAA7 -> RAA8 -> RAA9 -> RAA10 -> RAA13: all strict refinements OK")
    print()

    # test position-specific reduction
    print("Position-specific test:")
    seq = "LLLLLDDDDD"
    core = {0, 1, 2, 3, 4}  # first 5 positions
    result = reduce_at_positions(seq, core, RAA5)
    print(f"  Seq: {seq}, Core positions: {core}")
    print(f"  RAA5 at core only: {result}")

    print("\n=== ALL TESTS PASSED ===")
