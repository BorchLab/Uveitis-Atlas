# Data dictionary

## `inputs/data/metadata.csv`

One row per sequenced sample. The ingest step (`R/02_ingest_data.R`) reads it
and loads each sample from its CellRanger output directories.

| Column | Description |
|--------|-------------|
| `Sample` | Unique sample identifier (one 10x run). |
| `Subject` | Patient identifier. A subject can contribute multiple samples (eye and blood, multiple timepoints). |
| `Timepoint` | Visit / collection timepoint for the subject. |
| `Tube_ID` | Physical specimen tube identifier. |
| `Tissue_1` | Primary tissue. `Eye` (intraocular fluid) or `Blood`. Drives the eye-vs-blood contrasts. |
| `Tissue_2` | Secondary / finer tissue annotation where applicable. |
| `Cohort` | Collection cohort. |
| `scRNAseq_Submit_Date` | Submission date for the scRNA-seq run. |
| `Phenotype` | Raw clinical phenotype label. |
| `HLA_B27` | HLA-B27 genotype status. |
| `Etiology` | Specific disease etiology (e.g. VZV, HSV, sarcoid, idiopathic). |
| `Anatomy` | Anatomic site / uveitis classification. |
| `RNA_CRoutput` | Path to the CellRanger gene-expression output for this sample. |
| `TCR_CRoutput` | Path to the CellRanger TCR (VDJ) output. |
| `BCR_CRoutput` | Path to the CellRanger BCR (VDJ) output. |
| `Phenotype_2` | Two-level grouping used for the primary contrasts: `NIU` (non-infectious uveitis) vs `Viral`. |
| `Phenotype_3` | Finer phenotype grouping. |
| `Subject_Timepoint` | Convenience key combining subject and timepoint. |
| `Age` | Patient age. |
| `Sex` | Patient sex. |
| `Race` | Patient race. |
| `Ethinicity` | Patient ethnicity (column name as stored). |
| `Disease_Activity` | Clinical disease activity at collection. |
| `Symptom_Duration_Days` | Days of symptoms at collection. |
| `Disease_Duration_Years` | Years since disease onset. |

The `*_CRoutput` paths are environment specific. Remap them to your local
CellRanger outputs before running the ingest step.

## `outputs/` layout

| Path | Contents |
|------|----------|
| `outputs/objects/IntegratedSeuratObject.rds` | Full immune atlas. |
| `outputs/objects/eye/IntegratedSeuratObject.rds` | Eye sub-atlas. |
| `outputs/objects/{myeloid,bcell,tcell}/` | Per-compartment sub-atlases. |
| `outputs/tables/` | Full-atlas analysis CSVs. |
| `outputs/tables/eye/`, `outputs/tables/eye/bcell/` | Eye and per-compartment tables. |
| `outputs/tables/repertoire/` | TCR/BCR repertoire and clone-sharing tables. |
| `outputs/tables/cross_compartment/` | PC1 bridge + LIANA ligand-receptor tables. |
| `outputs/tables/stress_sensitivity/` | Sensitivity and robustness suite (see below). |
| `outputs/tables/stress_sensitivity/published_snapshot/` | Frozen copies of the nine published CSVs, written once. Every concordance metric is computed against these rather than against whatever a later pipeline run left on disk. |
| `outputs/tables/stress_sensitivity/jackknife/` | Leave-one-subject-out folds and summaries for the PC1 targets and the bridge. |
| `outputs/tables/stress_sensitivity/lodo/` | Leave-one-diagnosis-out folds for the NIU arm. |
| `outputs/tables/stress_sensitivity/floor_sensitivity/` | Bridge correlation across myeloid pseudobulk floors. |
| `outputs/tables/stress_sensitivity/pair_stability/` | Per-substate-pair significance across a permutation-count ladder. |
| `outputs/viz/` | Full-atlas figures. |
| `outputs/viz/eye/`, `outputs/viz/eye/{bcell,tcell,myeloid}/` | Eye and per-compartment figures. |
| `outputs/viz/cross_compartment/` | Myeloid-T cell coupling and LIANA figures. |
| `outputs/viz/program_contribution/` | Antiviral / autoimmune compartment-contribution panels. |

Within each target (full / eye / each compartment), figures are organized into a
single contiguous numbered subfolder scheme defined once in
`R/01_setup_utils.R::VIZ_BUCKETS` and resolved via `viz_subdir(paths, key)`:
`01_qc`, `02_integration`, `03_celltypes`, `04_markers`, `05_dge`, `06_escape`,
`07_milo`, `08_repertoire`, `09_tcr_motif`, `10_lineage_arch`, `11_pca_coupling`,
`12_composition`, `13_bcr_lineage`. A bucket folder is created only when a panel
is actually written into it (the save helpers create the parent on write), so
empty buckets never appear. `cross_compartment/` and `program_contribution/` are
cross-cutting (not per-target) and sit at the `outputs/viz/` root.
| `outputs/qc/` | QC reports. |

Per-table column schemas are documented in the header comment of the module
that writes each table (see the matching `R/NN_*.R` file).

## Sensitivity and robustness tables

Written by `R/19_stress_qc.R`, `R/46b_stress_sensitivity.R`,
`R/46c_pc1_jackknife.R` and `R/46d_diagnosis_lodo.R`. All are read-only with
respect to the published analysis.

### Per-cell scoring (`outputs/tables/eye/`)

| Table | One row per | Notes |
|-------|-------------|-------|
| `stress_ucell_per_cell.csv` | cell | UCell scores for the dissociation, inflammation-pruned, immediate-early, tissue-residency and glial panels, with `Subject`, `Phenotype_2`, `Tissue_2`, `Cohort`, `Etiology`, cluster and depth columns. |
| `stress_ucell_gene_coverage.csv` | gene set | Requested versus present gene counts, the `maxRank` used, and the missing-gene list, so any deviation from the published panel is auditable. |
| `stress_ucell_per_subject.csv` | (level, subject or sample) | Cell-count-weighted means. `level` distinguishes the pseudobulk unit from the independence unit. |
| `stress_ucell_per_pseudobulk.csv` | (target, substate, sample) | Mean score over exactly the cells forming each pseudobulk column. This is the covariate the adjusted models consume; the unit is (substate x sample), not sample. |
| `stress_ucell_group_test.csv` | score | Subject-level AUC, Hedges g, Wilcoxon p, plus the mixed-model estimate and subject ICC. `premise` records the pre-specified verdict on the reviewer's factual claim. |
| `stress_ucell_collinearity.csv` | (score, level, target, substate) | R-squared of the score on the group term, variance inflation, and a `flag` of ok / warn / abort. Above the abort threshold the adjusted arm is refused rather than reported. |
| `stress_trm_per_subject.csv`, `stress_trm_group_test.csv` | subject | Tissue-residency panel, restricted to the T cell compartment. |
| `stress_glial_per_subject.csv` | subject | Retinal-glial carryover as a fraction of cells above a within-dataset threshold, not a mean score. |

### Sensitivity outputs (`outputs/tables/stress_sensitivity/`)

| Table | One row per | Notes |
|-------|-------------|-------|
| `stress_sensitivity_verdict.csv` | (target, arm, metric) | Machine-computed MINOR / EQUIVOCAL / MATERIAL against config thresholds. `gating` marks which rows can trigger a re-run. |
| `stress_sensitivity_verdict_preamendment.csv` | as above | The first run's verdict, kept after the criteria were amended to be directional on 2026-08-01. |
| `dge_stress_gene_membership.csv` | (cluster, gene set) | Fraction of significant DEGs that are stress-panel members versus background, with a one-sided Fisher test. Makes no identifiability claim, so it is unaffected by the site/etiology collinearity. |
| `eye_dge_concordance_*.csv` | cluster | Adjusted versus published log2FC agreement, DEG Jaccard, sign-flip rate, top-50 retention. |
| `<target>_pca_concordance_*.csv` | (arm, substate) | PC1 correlation to published, Cohen's d before and after, separation retention, and the arm status and VIF actually applied. |
| `jackknife/pca_pc1_loo_*.csv` | (fold, substate) | Per-fold significance and `pc1_sign_flipped`, which detects an axis reversing on removal of one subject. |
| `jackknife/pc1_bridge_loo_*.csv` | (fold, weighting) | Per-fold partial r and its change from the full-data value, also split by whether the dropped subject was NIU or viral. |
| `lodo/lodo_*.csv` | (dropped diagnosis, substate) | Fold definitions carry a `fold_class` of informative, majority_bound or single_subject, because dropping half the arm is a power bound rather than a sensitivity test. |
| `floor_sensitivity/pc1_bridge_floor_sensitivity.csv` | (floor, weighting) | Bridge correlation with bootstrap CI and permutation p at each pseudobulk floor. |
| `pair_stability/pc1_pair_*.csv` | (pair, permutation setting) | `call_stable` is FALSE when a pair's significance call changes across permutation resolutions. |
