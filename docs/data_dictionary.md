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
