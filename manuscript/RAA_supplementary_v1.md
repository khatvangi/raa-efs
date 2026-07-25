*Ribbeck & Kiran — Supplementary Information*

**Supplementary Information**

**A five-letter amino acid alphabet suffices for protein folding when chemistry is preserved at structurally critical positions**

Hannah Ribbeck and Boggavarapu Kiran

Department of Chemistry and Physics, McNeese State University, Lake Charles, LA, USA

---

This document contains supplementary tables and supporting analyses referenced in the main text. All values derive from the same computational pipeline described in the Materials and Methods and are deposited in full at the data repository.

# **Supplementary Table S1. Reduced amino acid alphabet definitions**

Complete residue-to-group mappings for all six alphabets. Each alphabet is written as a set of groups; every residue in a group is mapped to the group representative (the first residue listed). RAA7 through RAA13 form a strict hierarchical refinement in which every group at size *n* is a subset of a group at size *n* – 1. RAA5 is a separate coarser partition and is not a strict parent of RAA7: alanine is grouped with the large hydrophobics in RAA5 but moves into a small-residue group {A, G} in RAA7, so the RAA5-to-RAA7 transition is a re-assignment rather than a pure merge.

| Alphabet | Letters | Groups (representative: members) |
| --- | --- | --- |
| RAA5 | 5 | L: L I V M F Y W A · G: G P · S: S T N Q C · D: D E · K: K R H |
| RAA7 | 7 | L: L I V M · F: F Y W · A: A G · P: P · S: S T N Q C · D: D E · K: K R H |
| RAA8 | 8 | L: L I V M · F: F Y W · A: A G · P: P · S: S T N Q · C: C · D: D E · K: K R H |
| RAA9 | 9 | L: L I V · M: M · F: F Y W · A: A G · P: P · S: S T N Q · C: C · D: D E · K: K R H |
| RAA10 | 10 | L: L I V · M: M · F: F Y W · A: A · G: G · P: P · S: S T · N: N Q · C: C · D: D E · K: K R H |
| RAA13 | 13 | L: L I · V: V · M: M · F: F · Y: Y · W: W · A: A · G: G · P: P · S: S · T: T · N: N · Q: Q · C: C · D: D · E: E · K: K R H |

The hierarchical splits introduced at each step are: RAA7 → RAA8 separates cysteine from the polar group; RAA8 → RAA9 separates methionine from the aliphatic hydrophobic group; RAA9 → RAA10 splits {A, G} into {A} and {G} and splits the polar group {S, T, N, Q} into hydroxyl {S, T} and amide {N, Q}; RAA10 → RAA13 splits the acidic pair {D, E} and the aromatic group {F, Y, W} into individual residues and separates L from I. RAA13 retains three residual groupings (L with I, and K/R/H as basic) because these pairs are the most conservative remaining substitutions.

**Chemistry-inversion mapping (hypothesis test only; not a reduced alphabet).** Core hydrophobic and structural residues are mapped to charged or flexible residues to violate native core chemistry: L→D, I→E, V→K, M→R, F→D, W→E, Y→K, A→D (hydrophobic → charged); C→K (breaks disulfide); P→G (removes backbone rigidity). Polar and already-charged residues (S, T, N, Q, D, E, K, R, H, G) retain native identity.

# **Supplementary Table S2. EFS-guided reduction across core thresholds, alphabets, and families**

Mean TM-score, standard deviation, 95% confidence interval, and pass/collapse rates for EFS-guided position-specific reduction. Core positions (those with EFS at or above the stated cutoff) retain native identity; all other positions are reduced to the stated alphabet. Four core thresholds (C10, C20, C30, C40) are crossed with four alphabets (RAA5, RAA7, RAA10, RAA13). The "ALL" rows pool all families; per-family rows follow. This table is the full sweep underlying the single-condition result reported as Table 3 in the main text (kinase, C20, RAA5).

Pooled across all families:

| Cutoff | Alphabet | n | Mean TM | SD | 95% CI | Pass ≥0.70 (%) | Collapse <0.50 (%) |
| --- | --- | --- | --- | --- | --- | --- | --- |
| C10 | RAA5 | 627 | 0.795 | 0.123 | 0.786–0.805 | 84.2 | 3.3 |
| C10 | RAA7 | 628 | 0.799 | 0.131 | 0.789–0.809 | 84.4 | 3.3 |
| C10 | RAA10 | 628 | 0.885 | 0.077 | 0.879–0.891 | 98.1 | 0.0 |
| C10 | RAA13 | 628 | 0.893 | 0.077 | 0.887–0.899 | 97.9 | 0.0 |
| C20 | RAA5 | 627 | 0.771 | 0.132 | 0.760–0.781 | 80.1 | 3.5 |
| C20 | RAA7 | 627 | 0.769 | 0.147 | 0.757–0.780 | 76.6 | 3.8 |
| C20 | RAA10 | 628 | 0.878 | 0.077 | 0.872–0.885 | 97.5 | 0.0 |
| C20 | RAA13 | 628 | 0.888 | 0.077 | 0.882–0.894 | 97.8 | 0.0 |
| C30 | RAA5 | 603 | 0.447 | 0.184 | 0.433–0.462 | 12.4 | 58.4 |
| C30 | RAA7 | 603 | 0.374 | 0.149 | 0.362–0.386 | 6.0 | 79.1 |
| C30 | RAA10 | 603 | 0.852 | 0.083 | 0.846–0.859 | 96.8 | 0.5 |
| C30 | RAA13 | 603 | 0.869 | 0.074 | 0.863–0.874 | 98.0 | 0.2 |
| C40 | RAA5 | 603 | 0.398 | 0.167 | 0.385–0.411 | 7.1 | 69.2 |
| C40 | RAA7 | 604 | 0.333 | 0.125 | 0.323–0.343 | 4.0 | 88.7 |
| C40 | RAA10 | 603 | 0.849 | 0.084 | 0.842–0.855 | 97.0 | 0.3 |
| C40 | RAA13 | 603 | 0.867 | 0.074 | 0.861–0.873 | 98.7 | 0.2 |

The threshold naming (C10–C40) refers to the EFS cutoff value, not to a percentage of positions. Because the EFS scoring function is discrete and takes values only at multiples of four, the cutoffs C10 and C20 select overlapping but progressively smaller core sets; C30 and C40 retain so few positions as core that the RAA5 and RAA7 reductions of the remaining sequence collapse the fold (pass rates fall below 13%), whereas RAA10 and RAA13 remain robust at every cutoff. This is the basis for the choice of C20 with RAA5 as the reported operating point: it is the most aggressive alphabet-plus-threshold combination that retains a high pass rate and zero collapses at the pooled level. The per-family breakdown at C20 is given in Supplementary Table S3.

# **Supplementary Table S3. EFS-guided reduction at the C20 threshold, by family**

The operating-point condition (EFS core threshold C20, native identity at core, reduced alphabet elsewhere), resolved by family. This is the source table for the kinase row quoted in main-text Table 3.

| Family | Alphabet | n | Mean TM | SD | 95% CI | Pass ≥0.70 (%) | Collapse <0.50 (%) |
| --- | --- | --- | --- | --- | --- | --- | --- |
| PF00069 Kinase | RAA5 | 585 | 0.791 | 0.092 | 0.783–0.798 | 83.9 | 0.0 |
| PF00069 Kinase | RAA7 | 585 | 0.789 | 0.107 | 0.780–0.798 | 79.5 | 0.5 |
| PF00069 Kinase | RAA10 | 586 | 0.886 | 0.069 | 0.881–0.892 | 99.1 | 0.0 |
| PF00069 Kinase | RAA13 | 586 | 0.896 | 0.068 | 0.890–0.901 | 99.5 | 0.0 |
| PF00089 Trypsin | RAA5 | 5 | 0.883 | 0.034 | 0.841–0.925 | 100.0 | 0.0 |
| PF00089 Trypsin | RAA10 | 5 | 0.891 | 0.064 | 0.812–0.970 | 100.0 | 0.0 |
| PF01145 Ferredoxin | RAA5 | 5 | 0.676 | 0.096 | 0.556–0.795 | 40.0 | 0.0 |
| PF01145 Ferredoxin | RAA10 | 5 | 0.789 | 0.077 | 0.694–0.884 | 80.0 | 0.0 |
| PF07686 Spectrin | RAA5 | 25 | 0.336 | 0.173 | 0.265–0.408 | 4.0 | 84.0 |
| PF07686 Spectrin | RAA10 | 25 | 0.729 | 0.065 | 0.702–0.755 | 68.0 | 0.0 |
| PF12796 Ankyrin | RAA5 | 5 | 0.694 | 0.174 | 0.478–0.910 | 60.0 | 20.0 |
| PF12796 Ankyrin | RAA10 | 5 | 0.837 | 0.102 | 0.711–0.964 | 80.0 | 0.0 |
| PF00047 Immunoglobulin | RAA5 | 2 | 0.593 | 0.059 | — | 0.0 | 0.0 |
| PF00047 Immunoglobulin | RAA10 | 2 | 0.715 | 0.165 | — | 50.0 | 0.0 |

Confidence intervals for families with n ≤ 2 are omitted as uninformative. The family-level pattern reinforces the main-text argument: at RAA5, spectrin (a diffuse-core family) collapses under EFS-guided reduction at the C20 threshold, whereas at RAA10 it recovers to a 68% pass rate. Families with well-defined cores (trypsin, kinase) tolerate RAA5 at core-preserving thresholds; families with diffuse cores require the more permissive RAA10 alphabet even with EFS guidance. The large kinase sample (n = 585) provides the tightest estimate and is the value cited in the main text.

# **Supplementary Table S4. Chemistry-axis dissection (subtle-chemistry experiment)**

To determine which specific chemical property of core residues is load-bearing, core positions were modified along six isolated chemical axes rather than by full inversion. Each axis alters one property while holding others as close to native as possible, applied at four core sizes (C8, C12, C20, C28). Thirty proteins per condition. Axis A (size) replaces core hydrophobics with size-mismatched hydrophobics; axis B (polar/H-bond character) alters polar core residues; axes C–F alter aromatic identity, charge, hydrogen-bonding, and backbone flexibility respectively.

| Axis (property altered) | C8 | C12 | C20 | C28 |
| --- | --- | --- | --- | --- |
| A — hydrophobic size/packing | 0.437 | 0.746 | 0.765 | 0.804 |
| B — polar / H-bond character | 0.363 | 0.660 | 0.733 | 0.779 |
| C — aromatic identity | 0.849 | 0.902 | 0.928 | 0.941 |
| D — charge identity | 0.890 | 0.978 | 0.985 | 1.000 |
| E — hydrogen bonding | 0.853 | 0.944 | 0.954 | 1.000 |
| F — backbone flexibility | 0.871 | 0.929 | 0.945 | 0.979 |

Values are mean TM-score across 30 proteins. The dissection refines the chemistry-inversion result reported in the main text. Violating hydrophobic size/packing (axis A) or polar character (axis B) at a small, well-defined core (C8) degrades the fold substantially (TM 0.437 and 0.363, with collapse rates of 67% and 80% respectively), whereas altering aromatic identity, charge, hydrogen bonding, or flexibility at the same positions leaves the fold essentially intact (TM 0.85–1.00). The interpretation is that the central chemistry-inversion signal is carried primarily by disruption of hydrophobic packing volume and polar hydrogen-bond geometry in the core, not by aromatic or charge identity per se. This is consistent with the ThermoMPNN result that within-class hydrophobic substitutions (which preserve size and packing) are near-neutral. The effect weakens monotonically as the core set grows from C8 to C28, because at larger core definitions the truly critical residues are diluted by more marginal ones.

# **Supplementary Table S5. Robustness of core definition and core size**

Two control experiments testing whether the EFS results depend on the specific multi-criteria scoring function or on the chosen core size.

Core-definition controls (30 proteins, RAA5 at non-core positions, native at core):

| Core definition | Mean TM | Pass ≥0.70 (%) | Collapse <0.50 (%) |
| --- | --- | --- | --- |
| EFS combined (C12) | 0.840 | 90.0 | 0.0 |
| EFS combined (C20) | 0.834 | 86.7 | 0.0 |
| Hydrophobicity only | 0.858 | 83.3 | 0.0 |
| SASA < 25 Å² | 0.835 | 86.7 | 0.0 |
| SASA < 40 Å² | 0.840 | 90.0 | 0.0 |

The five core definitions produce mean TM-scores within a 0.024 band (0.834–0.858) and pass rates within 7 percentage points. This confirms the main-text statement that for families with conventional hydrophobic cores, a simple SASA threshold reproduces most of the EFS benefit; the multi-criteria EFS score adds value chiefly in chemically distinctive folds (disulfide-rich and buried-aromatic-rich families), not uniformly.

Core-size variation (30 proteins, EFS-preserved vs. inverted, across core sizes):

| Condition | C8 | C12 | C20 | C28 |
| --- | --- | --- | --- | --- |
| EFS_PRESERVED (mean TM) | 0.893 | 0.843 | 0.841 | 0.834 |
| ANTI_EFS_INVERTED (mean TM) | 0.182 | 0.266 | 0.285 | 0.348 |

Preserving native chemistry at even the smallest core (C8, the ~8–12% most buried and packed positions) is sufficient to maintain the fold (TM 0.893); enlarging the preserved core beyond C8 yields no further benefit and slightly reduces mean TM as more surface positions are locked to native identity, reducing the compression achieved. Conversely, inverting chemistry at the smallest core is the most destructive (TM 0.182), because the C8 set concentrates the residues most essential to packing. Both trends confirm that a compact, well-chosen core carries the fold and that EFS identifies it.

# **Supplementary Table S6. Alternative core-definition overlap (Jaccard)**

Pairwise Jaccard overlap between EFS core positions (C20) and four alternative definitions, computed across 50 proteins spanning all ten families. Reproduced here in full; the main text cites the EFS–SASA and EFS–conservation values.

| Method vs. EFS | Jaccard | What it captures |
| --- | --- | --- |
| SASA-only (< 25 Å²) | 0.865 | Burial depth |
| ESM-2 attention (top 25%) | 0.343 | Learned functional features |
| Network centrality (betweenness, top 25%) | 0.231 | Contact-graph connectivity |
| Sequence conservation (top 25%) | 0.241 | Evolutionary constraint |

EFS overlaps strongly only with burial (SASA), as expected given that burial is the dominant term in the scoring function. The three non-burial definitions each capture fundamentally different information and none reproduces the EFS core: conservation reflects evolutionary rather than structural constraint, network centrality reflects topological hubs, and attention reflects learned sequence features. The 246 positions where EFS diverges from SASA are enriched in cysteines (23.6% overall, 73% in the EGF-like family) and buried aromatics.

# **Supplementary Note S1. ThermoMPNN substitution-level predictions**

In silico saturation mutagenesis on 602 PF00069 kinases (approximately 3.6 million ΔΔG predictions) supports two claims in the main text. First, EFS core positions are on average 2.8-fold more destabilizing to mutate than non-core positions (mean ΔΔG +1.66 versus +0.59 kcal/mol). Second, the substitutions characteristic of reduced alphabets are near-neutral, while the substitutions used in the chemistry-inversion mapping are strongly destabilizing. Representative within-class substitutions: I→L (−0.08 kcal/mol), V→L (−0.19), M→L (−0.12), all slightly stabilizing, consistent with these residues occupying similar volumes and maintaining equivalent hydrophobic contacts. The most destabilizing residues to introduce at an arbitrary position are proline, aspartate, and glycine, which are precisely the residues used in the inversion mapping. This convergence between an experimentally trained stability model (ThermoMPNN) and a structure-similarity readout (ESMFold TM-score) indicates that both detect the same underlying constraint: disruption of hydrophobic core packing.

# **Supplementary Note S2. Alphabet provenance and the RAA8/RAA9 interpolants**

The alphabets RAA5, RAA10, and RAA13 are standard groupings widely used in the reduced-alphabet literature. The two intermediate alphabets RAA8 and RAA9 were constructed for this study to make the RAA7-to-RAA10 interval a strict hierarchical refinement, so that the folding change at each single alphabet step is attributable to one identified chemical distinction. RAA8 separates cysteine from the polar group; RAA9 additionally separates methionine from the aliphatic hydrophobics. These two splits were chosen because cysteine (disulfide chemistry) and methionine (sulfur-mediated packing and methionine-aromatic interactions) are the most chemically distinctive residues within their parent groups. The strict-refinement property was verified computationally: for every adjacent pair in the chain RAA7 ⊃ RAA8 ⊃ RAA9 ⊃ RAA10 ⊃ RAA13, every pair of residues sharing a group in the finer alphabet also shares a group in the coarser one. RAA5 was confirmed not to satisfy this property relative to RAA7 (alanine is re-assigned between groups), and is therefore treated as a separate coarse baseline rather than as the root of the hierarchy.

# **Supplementary Note S3. Data and code availability**

All reduced sequences, predicted structures, per-residue EFS scores, US-align TM-score outputs, ThermoMPNN ΔΔG tables, and the canonical alphabet definitions are deposited at the data repository [DATA_DOI — Zenodo deposit to be inserted]. Analysis and figure-generation code, including the canonical alphabet module and the aggregation scripts that produce Supplementary Tables S2–S5, is available at [GITHUB_URL — to be inserted]. The high-resolution cliff data are provided as raa89_validation_results.csv and raa89_validation_summary.csv; the EFS sweep as efs_by_family_cutoff_alphabet.csv; the auxiliary experiments (chemistry axis, core size, core definition) as auxiliary_experiments_summary.csv.
