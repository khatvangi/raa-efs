# raa-efs — Essential Folding Score for position-specific amino acid alphabet reduction

Analysis code, figures, and result tables for:

**A five-letter amino acid alphabet suffices for protein folding when chemistry is preserved at
structurally critical positions**
Hannah Ribbeck and Boggavarapu Kiran, Department of Chemistry and Physics, McNeese State University.

## What this repository contains

| Directory | Contents |
|-----------|----------|
| `code/` | Canonical alphabet definitions and all analysis scripts (reduction, scoring, controls, visualization) |
| `figure_scripts/` | Scripts that regenerate every main and supplementary figure |
| `figures/` | Rendered figures (PNG + PDF) |
| `manuscript/` | Manuscript and Supplementary Information (Markdown source) |
| `results/` | Analysis-ready result tables (CSV) that reproduce every number in the paper |

The large primary data (reduced sequences, ~1,000s of predicted PDB structures, per-residue EFS scores,
ThermoMPNN ΔΔG tables) are archived separately on Zenodo: **[DATA_DOI — Zenodo deposit to be inserted]**.

## Core method: the Essential Folding Score (EFS)

`code/CANONICAL_ALPHABETS.py` is the single source of truth for all six reduced alphabets
(RAA5, RAA7, RAA8, RAA9, RAA10, RAA13) and the EFS definition. Run it directly to execute the
self-test (validates the strict hierarchical refinement RAA7 ⊂ RAA8 ⊂ RAA9 ⊂ RAA10 ⊂ RAA13):

```bash
python code/CANONICAL_ALPHABETS.py
```

EFS scores each residue by burial (SASA), contact density, and chemistry (bonuses for buried
hydrophobics, cysteines, and buried aromatics). Positions with EFS ≥ 20 are the structural core.

## Key result tables (`results/`)

| File | What it holds |
|------|---------------|
| `raa89_cliff_summary.csv` | High-resolution blind-reduction cliff curve, per family (Table 1, 2, Fig 1) |
| `efs_C20_summary.csv` | EFS-guided reduction at C20, per family × alphabet (Table 3, Supp S3) |
| `blind_vs_efs_RAA5.csv` | Blind vs EFS-guided RAA5, per family (Fig 3) |
| `efs_by_family_cutoff_alphabet.csv` | Full threshold × alphabet × family sweep (Supp S2) |
| `auxiliary_experiments_summary.csv` | Chemistry inversion + chemistry-axis + controls (Tables 5, Supp S4/S5) |
| `RAA_MASTER_RESULTS.csv` | Consolidated per-variant TM-scores |

## Reproducing the figures

```bash
pip install biopython pandas numpy matplotlib seaborn
python figure_scripts/fig_raa89_cliff_curve.py   # etc.
```

## Pipeline overview

1. `code/run_blind_all_alphabets.py` — blind (uniform) reduction of a FASTA to all six alphabets
2. `code/generate_efs_vs_sasa_reductions.py` — EFS-core vs SASA-core position-specific reductions
3. Fold reduced sequences with ESMFold (via the TRILL pipeline)
4. `code/score_reductions.py`, `code/score_efs_vs_sasa.py` — US-align scoring vs native
5. `figure_scripts/` — figures from the scored tables

## Citation

If you use this code or the EFS definition, please cite the paper (in preparation) and the Zenodo data
deposit.

## License

MIT (code). Data on Zenodo under CC-BY-4.0.
