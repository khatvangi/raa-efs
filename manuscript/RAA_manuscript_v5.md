*Ribbeck & Kiran — Position-specific alphabet reduction*

**A five-letter amino acid alphabet suffices for protein folding when chemistry is preserved at structurally critical positions**

Hannah Ribbeck and Boggavarapu Kiran*

*Department of Chemistry and Physics, McNeese State University, Lake Charles, LA, USA*

*Hannah Ribbeck ORCID: 0009-0005-3084-5607. Corresponding author: Boggavarapu Kiran (kiran@mcneese.edu), ORCID: 0000-0003-0751-6459.*

# **Abstract**

Reduced amino acid alphabets compress the twenty-letter protein code into fewer chemical classes, but three decades of work have applied this compression uniformly and have not asked where, between thirteen and five letters, fold integrity actually fails. We mapped the failure at single-alphabet resolution by folding 25 matched proteins across five Pfam families under six alphabets (RAA13 through RAA5). The collapse is sudden and family-specific: mean TM-score drops 0.135 units in a single step from RAA10 to RAA9, almost entirely driven by trypsin-like serine proteases (ΔTM = −0.409 in PF00089), where the merger of S/T with N/Q and of A with G eliminates hydrogen bonds essential to the β-protease topology. Position-specific reduction breaks this floor. Preserving native residues at positions identified by an Essential Folding Score (EFS)—a composite of burial, contact density, and chemistry—while reducing the remaining ~74% of positions to RAA5 restored folding on 602 kinases to a mean TM-score of 0.79 (84% pass rate, n = 585), approaching blind RAA13 performance with a five-letter alphabet at non-core sites. A chemistry-inversion experiment isolated the mechanism: substituting hydrophobic core residues with charged amino acids (L → D, I → E, V → K) collapsed 60% of folds, while conservative substitutions at the same positions preserved them (mean TM 0.85 vs. 0.38). The 0.47 TM-score gap between identical position sets under opposite chemical regimes, supported by ThermoMPNN saturation mutagenesis on 3.6 million substitutions showing a 2.8-fold ΔΔG penalty at EFS core positions, establishes that core chemical identity—not positional count or evolutionary conservation—determines fold survival under compression. The location of compression, not the alphabet size, sets the lower bound of protein foldability.

**Keywords:** *reduced amino acid alphabets, protein folding, Essential Folding Score, protein engineering, structural core, minimum viable alphabet*

# **Introduction**

Proteins encode their three-dimensional structures using a twenty-letter amino acid alphabet, but this alphabet is not uniformly necessary at every sequence position. Buried residues in protein cores tolerate conservative substitutions—leucine for isoleucine, valine for leucine—with minimal structural consequence, while surface residues evolve rapidly under weak structural constraint. This asymmetry is well established in evolutionary biology and protein engineering, yet its implications for a fundamental question in protein science remain underexplored: how far can the amino acid alphabet be compressed before protein folds collapse, and does the answer depend on where in the structure the compression is applied?

Experimental work over the past three decades has established rough bounds on the first half of this question. Riddle et al. (1997) demonstrated that the SH3 domain folds to near-native structure using only five amino acid types, with hydrophobic patterning at core positions preserved by design. Longo et al. (2013) showed that a β-trefoil domain composed predominantly of the ten prebiotically plausible amino acids adopts a stable, halophilic fold. Computational analyses by Murphy et al. (2000) and Solis (2019) converged on 10–12 letters as the practical optimum for preserving fold-recognition information. Throughout this literature, one design choice is essentially universal: a single reduced alphabet is applied to every position in the protein sequence. What has been missing is a systematic high-resolution map of where between the uniformly tolerated regime (13 letters) and the catastrophic regime (5 letters) the collapse actually occurs, and whether that collapse point is family-specific.

This gap matters because the uniform approach ignores a well-characterized physical reality. Protein cores are enriched in hydrophobic residues whose packing geometry and chemical character determine thermodynamic stability (Dill, 1990). Mutations at core positions are on average far more destabilizing than equivalent changes at the surface, a pattern confirmed both by experimental ΔΔG measurements and by machine-learning stability predictors such as ThermoMPNN (Dieckhaus et al., 2024). A reduction strategy that distinguishes core from surface should therefore outperform global reduction, and the advantage should grow exactly in the regime where uniform reduction begins to destroy core packing.

We developed such a strategy, built around the Essential Folding Score (EFS)—a composite metric that identifies structurally critical residues on the basis of solvent-accessible surface area, contact density, and amino acid chemistry, including explicit bonuses for cysteines (which form disulfide bonds regardless of burial) and for buried aromatic residues. EFS differs from conservation-based core definitions, which conflate functional selection pressures with structural ones. Conservation identifies positions under evolutionary constraint; EFS identifies positions under physical constraint.

To test this framework, we first mapped the blind-reduction cliff at fine resolution using six alphabets (RAA13, RAA10, RAA9, RAA8, RAA7, RAA5) applied to 25 proteins across five Pfam families, then asked whether EFS-guided reduction rescues fold integrity below the cliff. The two new alphabets introduced here (RAA8 and RAA9) were designed as strict hierarchical refinements that interpolate cleanly between RAA7 and RAA10, isolating cysteine (at RAA8) and methionine (at RAA9) as the first chemically distinctive residues to separate. This high-resolution blind scan identifies the exact compression ratio at which uniform reduction fails, and localizes the failure to specific protein families.

# **Results**

## ***Uniform alphabet reduction fails abruptly between RAA10 and RAA9***

To determine the precise compression ratio at which uniform alphabet reduction stops preserving protein folds, we applied six reduced alphabets—RAA13 (13 letters), RAA10 (10 letters), RAA9 (9 letters), RAA8 (8 letters), RAA7 (7 letters), and RAA5 (5 letters)—to a matched panel of 25 proteins (5 randomly selected per family from five Pfam families: kinase PF00069, trypsin-like serine protease PF00089, ferredoxin-like PF01145, spectrin/VWA PF07686, ankyrin repeat PF12796). The same 25 native structures were used at every alphabet size, removing any family-composition confound from the comparison. Each reduced sequence was folded with ESMFold (Lin et al., 2023) and compared to its native structure by TM-score using US-align (Zhang et al., 2022).

The alphabets from RAA7 through RAA13 form a strict hierarchical refinement: each group at size n is a subset of its parent group at size n – 1, so differences in folding outcome between adjacent alphabets in this range can be attributed to the specific chemical distinction introduced at that step. RAA5 is a coarser baseline that is not a strict parent of RAA7 (alanine regroups between the two), and the RAA7→RAA5 step is therefore treated as a re-assignment rather than a single merge. RAA7 → RAA8 separates cysteine from the polar group; RAA8 → RAA9 separates methionine from the aliphatic hydrophobic group; RAA9 → RAA10 splits {A, G} into two groups and divides polar residues into hydroxyl (S, T) and amide (N, Q) classes; RAA10 → RAA13 distinguishes the individual aromatic residues and the acidic pair. Full alphabet mappings appear in Supplementary Table S1.

The mean TM-score collapses abruptly at a specific step. At RAA13 the mean TM-score is 0.954, falling to 0.817 at RAA10—a modest loss of 0.137 TM-score units distributed over three chemical splits (aromatic partitioning, D/E separation, and several minor refinements). The next step, RAA10 → RAA9, drops the mean TM-score by 0.135 in a single move, from 0.817 to 0.682. This step is formally a single refinement: the merger of alanine with glycine and of S/T with N/Q. Continuing the reduction to RAA8 costs an additional 0.034 TM-score units (mean = 0.648), RAA7 a further 0.074 (mean = 0.574), and RAA5 a further 0.261 (mean = 0.313). The entire cliff between moderately tolerated compression and catastrophic collapse is therefore concentrated in the single step from RAA10 to RAA9, with the remainder of the decline distributed more gradually (Table 1, Figure 1).

**Table 1.** High-resolution blind reduction cliff curve. Mean TM-score across 25 proteins (5 per family, 5 families) for each of six reduced alphabets applied uniformly to every residue. ΔTM is the change from the preceding (less aggressive) alphabet. The step from RAA10 to RAA9 accounts for 21% of the total TM-score loss between RAA13 and RAA5 despite being a single chemical refinement.

| **Alphabet** | **Letters** | **Mean TM** | **SD** | **ΔTM vs. prev.** |
| --- | --- | --- | --- | --- |
| RAA13 | 13 | 0.954 | 0.048 | — |
| RAA10 | 10 | 0.817 | 0.149 | −0.137 |
| RAA9 | 9 | 0.682 | 0.220 | −0.135 |
| RAA8 | 8 | 0.648 | 0.210 | −0.034 |
| RAA7 | 7 | 0.574 | 0.242 | −0.074 |
| RAA5 | 5 | 0.313 | 0.134 | −0.261 |

The location of the cliff suggests that the chemical distinction encoded by the S/T versus N/Q split and the A versus G split carries substantially more structural information than any of the hydrophobic or aromatic splits above it. This seemed inconsistent with the conventional view that the hydrophobic core dominates fold specification. We therefore asked whether the cliff appears in every family or is driven by a subset of them.

## ***The RAA10 → RAA9 cliff is concentrated in trypsin-like proteases***

Decomposing the cliff curve by family reveals that the single-step collapse at RAA10 → RAA9 is overwhelmingly driven by trypsin-like serine proteases (PF00089), which fall from a mean TM-score of 0.876 at RAA10 to 0.467 at RAA9, a drop of 0.409 TM-score units in one alphabet step (Table 2). Ferredoxin-like domains and ankyrin repeats lose approximately 0.10 TM-score units at the same step, spectrin domains lose 0.05, and kinases lose only 0.026. If kinases alone had been used as the test set, as in our earlier work, the cliff would have appeared between RAA5 and RAA7 rather than between RAA9 and RAA10. Family choice determines where the cliff appears.

The trypsin-specific collapse identifies a chemical distinction that β-protease folds cannot afford to lose. The S/T group splits off from N/Q at RAA10: S and T carry a hydroxyl group and serve as H-bond donors, whereas N and Q carry an amide group and act as both donors and acceptors. Serine protease folds are densely cross-linked by backbone and side-chain hydrogen bonds, many of which involve serine and threonine hydroxyls in specific geometric configurations; merging these residues with the amide class destroys a large fraction of those contacts simultaneously. The A/G split that separates at the same step plays a secondary role: glycine's distinctive φ/ψ flexibility permits tight turns that alanine cannot adopt, and trypsin folds contain multiple such turns in their antiparallel β-sheet topology. The combination of these two losses at one alphabet step explains why the collapse is sudden rather than gradual for this family.

**Table 2.** Per-family decomposition of the blind reduction cliff at RAA10 → RAA9. Mean TM-scores are computed across five proteins per family. The drop is concentrated almost entirely in PF00089 trypsin-like proteases.

| **Family** | **TM at RAA10** | **TM at RAA9** | **ΔTM** |
| --- | --- | --- | --- |
| PF00089 Trypsin | 0.876 | 0.467 | −0.409 |
| PF01145 Ferredoxin | 0.844 | 0.745 | −0.099 |
| PF12796 Ankyrin | 0.731 | 0.635 | −0.096 |
| PF07686 Spectrin | 0.735 | 0.689 | −0.046 |
| PF00069 Kinase | 0.899 | 0.873 | −0.026 |

Two broader conclusions follow. First, RAA10 is the minimum alphabet size at which uniform reduction preserves fold integrity across all five families tested (cross-family mean TM-score = 0.817, no family mean below 0.73). Below RAA10, protein family becomes an important determinant of folding success, and no single uniform alphabet below ten letters is safe for all families. Second, the failure below RAA10 is not a generic "not enough chemical distinctions" effect; it is the loss of specific chemical features that specific folds depend on. This suggests that the correct remedy is not to find an even finer global alphabet but to apply the available alphabet selectively, preserving native identity at structurally critical positions while reducing elsewhere.

## ***EFS-guided position-specific reduction rescues fold integrity at RAA5***

We hypothesized that the RAA5 collapse reflects indiscriminate modification of structurally critical core residues rather than an intrinsic insufficiency of five chemical classes. To test this, we computed the Essential Folding Score for each residue in a larger panel of 602 PF00069 kinase structures and generated position-specific reduced variants: residues with EFS ≥ 20 (the C20 core threshold, encompassing approximately 26% of positions on average) retained their native amino acid identity, while all remaining positions were mapped to RAA5.

EFS-guided RAA5 variants achieved a mean TM-score of 0.791 with an 83.9% pass rate (TM ≥ 0.70) and zero collapses (n = 585 kinases with complete core scoring). The blind RAA5 score on the same protein family was 0.399 with a 7.1% pass rate and 68.9% collapse rate. Protecting roughly one-quarter of positions from reduction doubled the mean TM-score and increased the pass rate twelvefold (Table 3). Because both conditions apply the same five-letter alphabet and differ only in where it is applied, this comparison isolates placement from alphabet size: identical chemical resolution yields a collapsed fold when imposed uniformly and a native-like one when core positions are exempted (Figure 3). The C20 core threshold is one point on a threshold–alphabet surface: across the pooled dataset, EFS-guided reduction retains a pass rate above 80% and near-zero collapse for RAA5 and RAA7 at thresholds C10 and C20, while RAA10 and RAA13 remain robust at every threshold tested (C10 through C40); the full sweep with confidence intervals appears in Supplementary Tables S2 and S3.

**Table 3.** EFS-guided position-specific reduction versus blind uniform reduction on the PF00069 kinase family. EFS-guided RAA5 preserves native residues at positions with EFS ≥ 20 (approximately 26% of positions) and reduces the remainder; blind reductions apply the alphabet to every position.

| **Strategy** | **Mean TM** | **Pass rate (≥ 0.70)** | **Collapse rate (< 0.50)** |
| --- | --- | --- | --- |
| Blind RAA13 | 0.870 | 99.7% | 0.0% |
| Blind RAA10 | 0.879 | 98.8% | 0.0% |
| Blind RAA5 | 0.399 | 7.1% | 68.9% |
| EFS-guided RAA5 | 0.791 | 83.9% | 0.0% |

This result confirms that the RAA5 failure is not an intrinsic limit of five-letter chemical resolution. Five classes of amino acid are sufficient to maintain fold topology at surface and loop positions; the fold collapses only when core positions—those responsible for hydrophobic packing, disulfide bonding, and aromatic stacking—are subjected to the same coarse-grained mapping. The next questions concern what positions EFS actually identifies, whether the advantage generalizes beyond kinases, and what the mechanism of the core's contribution to fold stability is.

## ***EFS captures chemically important positions beyond simple burial***

To assess whether multi-criteria scoring adds value over simpler alternatives, we computed five core definitions for 50 proteins spanning 10 Pfam families: EFS (multi-criteria), SASA-only (burial below 25 Å²), network centrality (top 25% by betweenness in a Cα contact graph), ESM-2 attention (top 25% by transformer attention), and sequence conservation (top 25% most conserved columns in a family multiple sequence alignment). Pairwise overlap was quantified by the Jaccard similarity coefficient.

EFS and SASA-only overlapped substantially (Jaccard = 0.865), as expected: burial depth is the dominant contributor to EFS. The methods diverged on 246 positions across the dataset. Of these EFS-unique positions, 23.6% were cysteines—residues that are frequently solvent-exposed yet structurally essential through disulfide bonding (Figure 2). In the EGF-like domain family (PF00008), 73% of EFS-unique positions were cysteines, consistent with the known dependence of EGF folds on multiple disulfide cross-links. The remaining EFS-unique positions were enriched in aromatic residues with moderate burial and high contact density: residues contributing to core packing via ring stacking but falling below a strict SASA threshold.

Conservation showed the lowest overlap with EFS of any pairwise comparison (Jaccard = 0.241), reflecting the well-established distinction between evolutionary and structural constraint. A conserved active-site residue may be fully solvent-exposed; a buried hydrophobic residue may vary freely among L, I, V, and M because all four maintain packing. Network centrality (Jaccard = 0.231) and ESM-2 attention (0.343) captured fundamentally different structural information—connectivity hubs and learned functional features, respectively—confirming that no single alternative method recapitulates the EFS core definition (Supplementary Table S6).

## ***EFS-guided reduction outperforms conservation across diverse families***

To test EFS against its most natural competitor at scale, we conducted head-to-head comparisons on 20 proteins spanning all 10 families. For each protein, an EFS-reduced variant (native at EFS core, RAA5 elsewhere) and a conservation-reduced variant (native at top 25% conserved positions, RAA5 elsewhere) were folded and scored. EFS produced higher TM-scores in 13 of 20 cases (65%), with a mean advantage of +0.026 TM-score units (Table 4).

The pattern of wins was informative. EFS showed its largest advantages in families where conservation is least reliable as a proxy for structural importance. In spectrin/VWA domains (PF07686)—a family with diffuse hydrophobic content and low overall sequence conservation—EFS-reduced variants scored 0.601 versus 0.427 for conservation, a gap of +0.174 TM-score units, the largest single-family difference in the dataset. In ankyrin repeats (PF12796), where repeat architecture distributes structural constraints across many degenerate copies, EFS outperformed by +0.050. Conservation performed comparably in families with high sequence conservation (globins ΔTM = −0.011, trypsin ΔTM = −0.019), where evolutionary and structural constraints are already closely aligned. Three families—EGF-like domains, immunoglobulins, and HMG-CoA reductase—produced TM-scores below 0.50 with both methods, reflecting small disulfide-dependent or intrinsically disordered architectures that lack the stable hydrophobic cores EFS is designed to identify.

**Table 4.** Cross-family performance of EFS-guided versus conservation-guided reduction (RAA5 at non-core positions, n = 2 per family). ΔTM = TM(EFS) – TM(Conservation). Classes: I, TM ≥ 0.80; II, TM 0.60–0.80 with EFS advantage; III, both methods effective; IV, TM < 0.50.

| **Family** | **Pfam ID** | **Class** | **Core %** | **TM (EFS)** | **TM (Cons)** | **ΔTM** |
| --- | --- | --- | --- | --- | --- | --- |
| TIM barrel | PF00121 | I | 44.0 | 0.980 | 0.968 | +0.011 |
| Trypsin | PF00089 | I | 38.8 | 0.820 | 0.838 | −0.019 |
| Ferredoxin | PF01145 | I | 14.9 | 0.819 | 0.794 | +0.025 |
| Kinase | PF00069 | II | 25.7 | 0.745 | 0.729 | +0.015 |
| Spectrin/VWA | PF07686 | II | 17.6 | 0.601 | 0.427 | +0.174 |
| Ankyrin | PF12796 | II | 28.8 | 0.654 | 0.604 | +0.050 |
| Globin | PF00042 | III | 37.3 | 0.944 | 0.954 | −0.011 |
| EGF-like | PF00008 | IV | 18.2 | 0.306 | 0.316 | −0.010 |
| Immunoglobulin | PF00047 | IV | 21.9 | 0.395 | 0.366 | +0.029 |
| HMG-CoA red. | PF00505 | IV | 8.2 | 0.300 | 0.303 | −0.003 |

The cross-family results demonstrate that EFS-guided reduction is superior when structural and functional constraints diverge, equivalent when they coincide, and insufficient when no stable core exists. But a subtler observation warranted investigation. Standard RAA reduction of core positions (not the full inversion, but a within-class substitution: I → L, V → L, and so on) still preserved structure with mean TM-score of 0.850 in pilot experiments, even though it changed the identity of every core residue. If core positions are critical, why does replacing them with chemically similar amino acids leave the fold intact? This raised the hypothesis that chemical identity, not positional specificity per se, is the operative variable.

## ***Core chemistry, not positional identity, determines fold survival***

To decouple positional and chemical effects, we designed an experiment in which core positions were modified under two different chemical regimes. Thirty proteins from three families—PF00069 kinases (37% core fraction), PF01145 ferredoxin-like domains (77%), and PF07686 spectrin/VWA domains (57%)—were subjected to five conditions. In the EFS_PRESERVED condition, core positions retained native residues and non-core positions were reduced to RAA10. In ANTI_EFS_REDUCED, core positions were reduced to RAA10 (chemically conservative substitutions) while non-core residues retained native identity. In ANTI_EFS_INVERTED, core hydrophobics were replaced with charged residues (L → D, I → E, V → K, F → D), a deliberate chemical violation of the hydrophobic core; non-core residues remained native. BLIND_REDUCED applied RAA10 uniformly to all positions. All 150 variants were folded with ESMFold.

The results distinguish the two hypotheses cleanly. EFS_PRESERVED variants maintained a mean TM-score of 0.879 (90% pass rate). ANTI_EFS_REDUCED—the same core positions modified, but with conservative substitutions—achieved 0.850 (80% pass rate). Chemistry inversion of those same positions collapsed the fold to a mean TM-score of 0.380, with 60% of structures below TM = 0.50 (Table 5). Every one of the 30 proteins tested showed worse outcomes under inversion than under standard reduction. The mean gap between these two treatments of identical positions was 0.47 TM-score units (Figure 4). Controls confirming that these outcomes are insensitive to the precise core definition and to the chosen core size are reported in Supplementary Table S5.

Collapse severity tracked core size. Ferredoxin-like domains (77% core) disintegrated most completely under inversion (mean TM = 0.189), because the majority of the sequence underwent chemical violation. Kinases (37% core) retained partial structure (mean TM = 0.680), consistent with the non-core portion still providing structural information. This dose-response relationship rules out alternative explanations and establishes that the hydrophobic character of core residues—not their sequence identity or positional information content—is what holds protein folds together.

**Table 5.** Chemistry inversion experiment (n = 30 proteins, 3 families). The critical comparison is between ANTI_EFS_REDUCED and ANTI_EFS_INVERTED, which modify identical positions with different chemical outcomes.

| **Condition** | **Mean TM** | **Pass rate (≥ 0.70)** | **Collapse rate (< 0.50)** |
| --- | --- | --- | --- |
| Native | ~0.95 | 100% | 0% |
| EFS_PRESERVED | 0.879 | 90.0% | 0.0% |
| ANTI_EFS_REDUCED | 0.850 | 80.0% | 0.0% |
| BLIND_REDUCED | 0.815 | 76.7% | 0.0% |
| ANTI_EFS_INVERTED | 0.380 | 16.7% | 60.0% |

## ***Which chemical property matters: hydrophobic packing and polar geometry, not aromatic or charge identity***

The inversion experiment shows that violating core chemistry collapses the fold, but inversion changes several chemical properties at once (hydrophobicity, charge, size). To isolate the responsible property, we modified core positions along six chemical axes independently—hydrophobic size and packing, polar and hydrogen-bond character, aromatic identity, charge identity, hydrogen bonding, and backbone flexibility—each applied at four core sizes (30 proteins per condition; Supplementary Table S4). The axes separate cleanly. Violating hydrophobic size and packing, or polar hydrogen-bond character, at the smallest well-defined core (C8) degrades the fold substantially (mean TM 0.437 and 0.363, with 67% and 80% collapse respectively). Altering aromatic identity, charge, hydrogen bonding, or backbone flexibility at the same positions leaves the fold essentially intact (mean TM 0.85 to 1.00). The chemistry-inversion signal is therefore carried primarily by two properties: the packing volume of buried hydrophobic side chains and the hydrogen-bond geometry of core polar residues. Aromatic ring identity and formal charge, which the coarse alphabets also alter, contribute little to fold survival at the core. This is consistent with the ThermoMPNN finding below that within-class hydrophobic substitutions—which preserve size and packing—are near-neutral, and it explains why RAA5, which merges all hydrophobics to a single class while preserving their hydrophobic character, succeeds at non-core positions.

## ***ThermoMPNN independently validates EFS core positions as stability-critical***

The chemistry inversion experiment established that core chemical identity drives fold survival, but relied on ESMFold as the sole structural readout. To obtain an independent assessment, we performed in silico saturation mutagenesis on all 602 PF00069 kinase structures using ThermoMPNN (Dieckhaus et al., 2024), a model trained on experimental thermodynamic data rather than on structural similarity. Each non-terminal position was mutated to all 20 amino acids, yielding approximately 3.6 million ΔΔG predictions.

EFS core positions (score ≥ 20) showed a mean predicted destabilization of +1.66 kcal/mol upon mutation, compared to +0.59 kcal/mol for non-core positions—a 2.8-fold difference (Figure 5). Because ThermoMPNN's training objective is independent of ESMFold's, this confirms that EFS selects positions where amino acid identity has the greatest impact on stability.

The substitution-level predictions explained both why standard RAA reduction works and why chemistry inversion fails. Within-group substitutions characteristic of reduced alphabets—I → L (ΔΔG = −0.08 kcal/mol), V → L (−0.19), M → L (−0.12)—were predicted to be slightly stabilizing, consistent with these residues occupying similar volumes and maintaining equivalent hydrophobic contacts. The amino acids most destabilizing upon introduction at any position—proline (+1.58), aspartate (+1.25), glycine (+1.20)—are precisely those used in the chemistry inversion mapping. The convergence between ThermoMPNN predictions and folding outcomes indicates that both methods detect the same underlying constraint: disruption of hydrophobic core packing.

These data also resolved a historical anomaly. An early RAA7 alphabet contained an erroneous C → D mapping that produced results worse than RAA5—a counterintuitive outcome for a less aggressive reduction. ThermoMPNN confirmed that introducing aspartate at cysteine positions destroys the disulfide cross-link and introduces a buried negative charge, compounding two destabilizing effects. Correcting the mapping to preserve cysteine identity (as in the canonical RAA7 used throughout this paper) resolved the anomaly and motivated the design of RAA8, which separates cysteine into its own group as the first refinement above RAA7.

## ***Protein families cluster into four reducibility classes***

The 10-family dataset revealed a natural classification based on EFS-guided reduction performance (Figure 6). Three Class I families (TIM barrel, trypsin, ferredoxin) achieved TM ≥ 0.80 under EFS-guided RAA5—large, well-defined hydrophobic cores where EFS operates in its strongest regime. Three Class II families (kinase, spectrin, ankyrin) achieved TM 0.60–0.80 with EFS consistently outperforming conservation; these represent the regime where structure-aware core identification provides the greatest added value over naive conservation. Globins (Class III) represent an inherently easy case where both approaches succeed. Three Class IV families (EGF, immunoglobulin, HMG-CoA reductase) resisted reduction by any method.

Core percentage did not predict success. Ferredoxin-like domains achieve TM = 0.819 with only 14.9% of positions classified as core; trypsin achieves a comparable TM = 0.820 with 38.8% core. HMG-CoA reductase fails with 8.2% core, but it fails because the domain is intrinsically disordered, not because the core is small (Weathers et al., 2004). The critical determinant is whether the protein possesses a well-packed hydrophobic core at all—not what fraction of the sequence that core comprises.

# **Discussion**

Three decades of work on reduced amino acid alphabets have operated under a shared assumption: the same mapping is applied everywhere. Riddle et al. (1997) used five amino acids throughout the SH3 sequence, placing hydrophobics at core positions by manual design. Longo et al. (2013) built entire domains from prebiotically plausible sets. Murphy et al. (2000) optimized global alphabets against substitution matrices. The high-resolution cliff curve presented here shows that this uniform strategy hits a sharp, family-dependent limit—not at five letters, where the literature has focused, but at nine. For the five families tested, RAA10 is the minimum uniform alphabet that preserves fold integrity across all of them; the single-step reduction to RAA9 drops the cross-family mean TM-score by 0.135 and collapses trypsin-like proteases by 0.409. Position-specific reduction guided by a structure-based score breaks through this limit, restoring fold integrity at RAA5 to levels achieved only above RAA10 under uniform reduction.

The localization of the cliff to the RAA10 → RAA9 step was not anticipated. The specific chemical distinctions lost at this step—S/T versus N/Q, and A versus G—are not among the most prominent discussions in the reduced-alphabet literature, which has focused on aromatic partitioning, the hydrophobic core, and the acidic and basic groups. The S/T split separates H-bond donors from donor-acceptor amides; the A/G split separates a residue with a minimal side chain (A) from one with a unique backbone flexibility (G). Both distinctions are critical for the tight turns and hydrogen-bond networks of β-protease folds, where a single misassigned side chain can disrupt a network of specific contacts. This family-specific sensitivity suggests that future work on reduced alphabets should test on structurally diverse panels rather than optimizing against aggregate metrics, which can hide sharp family-level failures.

The chemistry inversion experiment makes the mechanism behind the EFS rescue clear. The same core positions that tolerate I → L and V → L (TM = 0.85) collapse catastrophically under L → D and V → K (TM = 0.38). Both treatments change the same number of positions. They differ only in whether the substitutions preserve hydrophobic character. The fold's survival therefore depends on the physicochemical environment of the core—not on which positions are changed or how many letters the alphabet contains. This interpretation aligns with the classical hydrophobic collapse model (Dill, 1990) and with ThermoMPNN predictions showing that cross-class substitutions are 2–3-fold more destabilizing than within-class substitutions at core positions. A per-axis dissection (Supplementary Table S4) localizes the effect further: of the chemical properties that inversion disrupts, only hydrophobic packing volume and polar hydrogen-bond geometry are load-bearing at the core, while aromatic identity and formal charge are largely dispensable—which is precisely why an alphabet that collapses all hydrophobics into one class, but keeps them hydrophobic, preserves the fold.

EFS adds modest but meaningful value over simple SASA-based burial (Jaccard = 0.865 overlap). For most families, a SASA < 25 Å² threshold would produce nearly equivalent results. The 13.5% of positions where EFS diverges from SASA are enriched in cysteines and buried aromatics—residues whose structural contributions arise from mechanisms other than volume exclusion. For disulfide-rich families like EGF domains (73% cysteine among EFS-unique positions), this distinction matters. For families with conventional hydrophobic cores, it is marginal. We do not claim EFS is always necessary; we claim it is never worse and is measurably better in chemically distinctive folds.

Two limitations bound the interpretation. First, every folding measurement here is from ESMFold and every stability estimate from ThermoMPNN; neither is an experimental observable, and reduced sequences fall outside the statistical distribution on which both models were trained. The magnitudes of the effects reported—a 0.135 ΔTM at the uniform cliff, a 0.47 ΔTM between identical position sets under conservative and inverted chemistry, a 2.8-fold ΔΔG penalty at core positions—exceed any reasonable estimate of model bias on reduced sequences, but the absolute TM-scores and the precise location of the cliff for under-sampled families warrant experimental confirmation. The natural next step is circular dichroism and thermal denaturation of EFS-guided RAA5 variants in the kinase and trypsin families, paired with crystallographic characterization of a single rescued variant, to test whether high-TM in silico predictions correspond to native-like secondary structure and cooperative unfolding in vitro. Second, the cliff scan used five proteins per family (25 total), sufficient to resolve the trypsin-specific collapse but underpowered to localize family-specific cliffs in non-protease folds; expanding this scan to 20 proteins per family is a straightforward computational extension.

The most immediate application is library design for directed evolution. Combinatorial libraries built from the natural twenty face an exponential sample-space problem that grows with each mutated position; reducing the alphabet to five letters at non-core positions while preserving native chemistry at the ~26% of positions EFS identifies compresses a million-variant library to roughly a hundred while retaining fold integrity. Two related applications follow. Ancestral and prebiotic protein reconstruction, building on Longo et al. (2013), can be guided position-by-position by EFS—keeping modern chemistry at structurally essential sites and reducing the rest to a prebiotic subset, rather than committing to a single alphabet for the entire sequence. And compressed protein representations for machine learning should treat ten letters as the floor for cross-family generalization; below RAA10, model performance will degrade preferentially on β-protease and other folds where the S/T and A/G distinctions are load-bearing. The next experimental step is to construct an EFS-guided RAA5 library for a small, fast-folding scaffold—protein L or a streptavidin variant—and measure the rescued fraction by yeast display and thermal denaturation. If the in vitro hit rate matches the predicted 85%, position-specific reduction becomes a usable tool for compressing protein design search spaces at the bench.

# **Methods**

## ***Dataset construction***

Protein structures were obtained from the AlphaFold Protein Structure Database (Varadi et al., 2022) and the Protein Data Bank (Berman et al., 2000) for 10 Pfam families: PF00008 (EGF-like), PF00042 (globin), PF00047 (immunoglobulin), PF00069 (protein kinase), PF00089 (trypsin-like serine protease), PF00121 (TIM barrel), PF00505 (HMG-CoA reductase), PF01145 (ferredoxin-like), PF07686 (spectrin/VWA), and PF12796 (ankyrin repeat). Structures with experimental resolution worse than 3.0 Å or AlphaFold pLDDT below 70 were excluded. The full dataset comprised 2,652 native structures. For the high-resolution cliff experiment, five proteins were randomly selected from each of five families (PF00069, PF00089, PF01145, PF07686, PF12796) for a total of 25 proteins tested at every alphabet size. For the 10-family EFS-versus-conservation comparison, five proteins per family were selected to span structural complexity, with two per family used in the head-to-head comparison.

## ***Reduced amino acid alphabets***

Six reduced alphabets were used. The five alphabets from RAA7 to RAA13 form a strict hierarchical refinement (each group is a subset of a group in the next-finer alphabet); RAA5 is a separate coarser baseline in which alanine is grouped with the large hydrophobics rather than with glycine, so it is not a strict parent of RAA7. RAA5 maps hydrophobics (L, I, V, M, F, Y, W, A) to L, conformationally special residues (G, P) to G, polar (S, T, N, Q, C) to S, acidic (D, E) to D, and basic (K, R, H) to K. RAA7 separates aromatics (F, Y, W → F), preserves cysteine with the polar group at this level, isolates proline, and groups small residues (A, G → A). RAA8 (introduced here) separates cysteine into its own group. RAA9 (introduced here) additionally separates methionine from the aliphatic hydrophobic group (L, I, V). RAA10 additionally splits the polar group into hydroxyl (S, T) and amide (N, Q), and separates alanine from glycine. RAA13 further distinguishes the aromatic residues (F, Y, W), the acidic pair (D, E), and several minor refinements. For every adjacent pair in the chain RAA7 ⊂ RAA8 ⊂ RAA9 ⊂ RAA10 ⊂ RAA13, every residue pair sharing a group in the finer alphabet also shares a group in the coarser one; RAA5 does not satisfy this relative to RAA7 and is treated as a distinct baseline (complete mappings in Supplementary Table S1). The chemistry inversion mapping replaces core hydrophobics with charged residues (L → D, I → E, V → K, M → R, F → D, W → E, Y → K, A → D, C → K, P → G) while leaving polar and charged residues at native identity.

## ***Essential Folding Score***

EFS was computed per residue using a discrete scoring function. Solvent-accessible surface area (SASA) was calculated with a 1.4 Å probe radius using BioPython (Cock et al., 2009). Contact density was defined as the number of Cα atoms within 8 Å. The raw score accumulated: deep burial (SASA < 20 Å², +3), moderate burial (SASA 20–40 Å², +1), high packing (contacts ≥ 8, +2), moderate packing (contacts 5–7, +1), buried hydrophobic bonus (residue in {L, I, V, F, M, W, Y, A} and SASA < 30 Å², +2), cysteine identity (+2 regardless of burial), and buried aromatic bonus (residue in {F, W, Y} and SASA < 30 Å², +1). The raw score was multiplied by 4 (scale 0–36). Positions with EFS ≥ 20 (raw ≥ 5) were classified as core.

## ***Sequence reduction and structure prediction***

For uniform (blind) reduction, every position in the input sequence was mapped to its reduced equivalent. For EFS-guided reduction, positions with EFS ≥ 20 retained native identity; remaining positions were mapped. For conservation-guided reduction, the top 25% most conserved positions (by Shannon entropy from MUSCLE alignments; Edgar, 2004) retained native identity. All variants were folded using ESMFold (Lin et al., 2023) via the TRILL pipeline (Martinez et al., 2023) with default parameters.

## ***Structure comparison***

Predicted structures were aligned to native structures using US-align (Zhang et al., 2022). TM-score ≥ 0.70 was used as the pass threshold (same fold); TM < 0.50 was classified as collapse. TM-score values reported throughout are normalized to the native structure length.

## ***ThermoMPNN saturation mutagenesis***

In silico saturation mutagenesis was performed on 602 PF00069 kinases using ThermoMPNN (Dieckhaus et al., 2024). Each non-terminal residue was mutated to all 20 amino acids and ΔΔG (kcal/mol) recorded. Approximately 3.6 million predictions were generated.

## ***Alternative core definitions***

SASA-only: positions with SASA < 25 Å². Network centrality: Cα contact graph (8 Å cutoff) via NetworkX (Hagberg et al., 2008), top 25% by betweenness. ESM attention: ESM-2 model (esm2_t33_650M_UR50D; Lin et al., 2023), top 25% by attention score. Conservation: Shannon entropy per column from MUSCLE alignment, top 25% by conservation (1 minus normalized entropy). Pairwise overlap quantified by Jaccard similarity.

## ***Data and code availability***

All reduced sequences, predicted structures, TM-scores, EFS core scores, ThermoMPNN predictions, and canonical alphabet definitions are deposited at [DATA_DOI — Zenodo deposit to be inserted here]. The high-resolution cliff experiment raw data are available as raa89_validation_results.csv (150 per-protein TM-scores) and raa89_validation_summary.csv (aggregated per family-alphabet). Analysis code is available at [GITHUB_URL — insert finalized URL here].

# **References**

1. Berman, H.M. et al. The Protein Data Bank. Nucleic Acids Res. 28, 235–242 (2000).

2. Cock, P.J.A. et al. Biopython: freely available Python tools for computational molecular biology. Bioinformatics 25, 1422–1423 (2009).

3. Dieckhaus, H. et al. Transfer learning to leverage larger datasets for improved prediction of protein stability changes. Proc. Natl. Acad. Sci. 121, e2314853121 (2024).

4. Dill, K.A. Dominant forces in protein folding. Biochemistry 29, 7133–7155 (1990).

5. Edgar, R.C. MUSCLE: multiple sequence alignment with high accuracy and high throughput. Nucleic Acids Res. 32, 1792–1797 (2004).

6. Hagberg, A.A., Schult, D.A. & Swart, P.J. Exploring network structure, dynamics, and function using NetworkX. Proc. 7th Python in Science Conf. 11–15 (2008).

7. Lin, Z. et al. Evolutionary-scale prediction of atomic-level protein structure with a language model. Science 379, 1123–1130 (2023).

8. Longo, L.M., Lee, J. & Blaber, M. Simplified protein design biased for prebiotic amino acids yields a foldable, halophilic protein. Proc. Natl. Acad. Sci. 110, 2135–2139 (2013).

9. Martinez, Z.A., Murray, R.M. & Thomson, M.W. TRILL: Orchestrating modular deep-learning workflows for democratized, scalable protein analysis and engineering. bioRxiv 2023.10.24.563881 (2023). doi:10.1101/2023.10.24.563881

10. Murphy, L.R., Wallqvist, A. & Levy, R.M. Simplified amino acid alphabets for protein fold recognition. Protein Eng. 13, 149–152 (2000).

11. Riddle, D.S. et al. Functional rapidly folding proteins from simplified amino acid sequences. Nat. Struct. Biol. 4, 805–809 (1997).

12. Solis, A.D. Reduced alphabet of prebiotic amino acids optimally encodes the conformational space of diverse extant protein folds. BMC Evol. Biol. 19, 158 (2019).

13. Varadi, M. et al. AlphaFold Protein Structure Database. Nucleic Acids Res. 50, D419–D427 (2022).

14. Weathers, E.A. et al. Reduced amino acid alphabet is sufficient to accurately recognize intrinsically disordered protein. FEBS Lett. 576, 348–352 (2004).

15. Zhang, C., Shine, M., Pyle, A.M. & Zhang, Y. US-align: universal structure alignments. Nat. Methods 19, 1109–1115 (2022).

# **Figure Legends**

**Figure 1.** High-resolution blind reduction cliff. (a) Mean TM-score across 25 proteins (5 per family, 5 families) for each of six uniform alphabets applied to every residue. The single largest drop (ΔTM = −0.135) occurs between RAA10 and RAA9, identifying the minimum viable blind alphabet. (b) Per-family decomposition of the cliff: the RAA10 → RAA9 drop is concentrated in PF00089 trypsin (−0.409), modest for ferredoxin, ankyrin, and spectrin, and nearly absent for kinases. Dashed horizontal line: TM = 0.70 pass threshold.

**Figure 2.** Pipeline schematic and EFS construction. (a) EFS scoring identifies core positions by burial, contact density, and amino acid chemistry; non-core positions are reduced to RAA5; ESMFold predicts structure; TM-score validates. (b) Overlap between EFS and SASA-only core definitions (Jaccard = 0.865); 246 EFS-unique positions are indicated. (c) Amino acid composition of EFS-unique positions, showing 23.6% cysteine enrichment across all families and 73% in EGF-like domains.

**Figure 3.** EFS-guided RAA5 reduction rescues folding. Blind RAA5 on 602 kinases (mean TM = 0.399, 69% collapse) versus EFS-guided RAA5 on the same proteins (mean TM = 0.791, 0% collapse). Dashed line at TM = 0.70 marks the pass threshold. Preserving 26% of positions at native identity doubles the mean TM-score and eliminates catastrophic failures.

**Figure 4.** Chemistry inversion collapses folds. Mean TM-scores for five conditions (n = 30 proteins, 3 families), with individual scores overlaid. The 0.47 TM-score gap between ANTI_EFS_REDUCED (0.850) and ANTI_EFS_INVERTED (0.380) demonstrates that chemical identity at core positions, not positional information, determines fold survival.

**Figure 5.** ThermoMPNN validates EFS core positions. (a) ΔΔG distributions for EFS core vs. non-core positions across 602 kinases (3.6 million predictions). Core: +1.66 kcal/mol; non-core: +0.59 kcal/mol. (b) Mean ΔΔG by substitution identity; RAA-relevant substitutions (I → L, V → L) are near zero while inversion-relevant substitutions (any → D, any → P) are highly destabilizing.

**Figure 6.** Cross-family reducibility classification. Ten families grouped by response to EFS-guided RAA5 reduction. Left: mean TM by family, colored by class (I–IV). Center: ΔTM (EFS minus conservation). Right: structural properties associated with each class.
