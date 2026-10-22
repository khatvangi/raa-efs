"""Essential Folding Score: where a protein sequence can and cannot be compressed.

Typical use::

    from raa_efs import score_structure, reduce_by_efs

    scored = score_structure("1ake.pdb")
    variant = reduce_by_efs(scored, "RAA5", threshold=20)
    print(variant.sequence, variant.preserved_fraction)

and, for the control that makes that result interpretable::

    from raa_efs import random_control_matching

    control = random_control_matching(scored, "RAA5", seed=1)
"""

from .alphabets import (
    ALIASES,
    CANONICAL,
    LEGACY,
    STANDARD_AA,
    AlphabetError,
    get as get_alphabet,
    identify as identify_alphabet,
    output_letters,
    size as alphabet_size,
)
from .reduce import (
    ReductionResult,
    identify_alphabet_from_reduction,
    random_control_matching,
    random_matched_control,
    reduce_by_efs,
    reduce_position_specific,
    reduce_uniform,
    summarise,
)
from .score import (
    DEFAULT_CORE_THRESHOLD,
    METHODS_AS_WRITTEN,
    PROFILES,
    PUBLISHED,
    ResidueScore,
    ScoredChain,
    ScoringProfile,
    StructureError,
    raw_score,
    score_sequence_features,
    score_structure,
)

__version__ = "2.0.0"

__all__ = [
    "ALIASES",
    "CANONICAL",
    "DEFAULT_CORE_THRESHOLD",
    "LEGACY",
    "METHODS_AS_WRITTEN",
    "PROFILES",
    "PUBLISHED",
    "STANDARD_AA",
    "AlphabetError",
    "ReductionResult",
    "ResidueScore",
    "ScoredChain",
    "ScoringProfile",
    "StructureError",
    "alphabet_size",
    "get_alphabet",
    "identify_alphabet",
    "identify_alphabet_from_reduction",
    "output_letters",
    "random_control_matching",
    "random_matched_control",
    "raw_score",
    "reduce_by_efs",
    "reduce_position_specific",
    "reduce_uniform",
    "score_sequence_features",
    "score_structure",
    "summarise",
    "__version__",
]
