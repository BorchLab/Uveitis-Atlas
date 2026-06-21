# Data availability

This GitHub repository holds **code only**. All data are distributed through the
project **Zenodo archive** (DOI: _coming soon_). After cloning, download the
archive and unpack each part into the matching directory below. The pipeline
expects exactly this layout (paths are set in `config/config.yml`).

## Where each Zenodo part goes

The deposit is focused on data delivery (~30 GB total, one Zenodo record).

| Zenodo part | Unpack into | Approx size | Needed for |
|-------------|-------------|-------------|------------|
| Analysis objects (`.rds`) | `outputs/objects/` | ~16 GB | Analysis + figure regeneration |
| Per-sample QC reports | `outputs/qc/` | ~65 MB | QC review |
| Raw CellRanger per-sample outputs | `inputs/data/runs/` | ~13 GB | Rebuild atlases from scratch (`R/02`) |
| VDJdb reference cache | `references/vdjdb/` | ~464 MB | TCR annotation (`R/69`) |

> **Not distributed.** To keep the deposit lean:
> - **Tables (`outputs/tables/`) and figures (`outputs/viz/`)** — published in the
>   manuscript; regenerate them from `outputs/objects/` via the viz/analysis steps.
> - **Intermediate per-sample objects (`outputs/processed/`)** — regenerated from
>   `inputs/data/runs/` via `R/02` + `R/03`.
> - **Ibex/Trex encodings and escape ssGSEA chunks** (sidecar files that would
>   otherwise sit in `outputs/objects/`) — unused in the paper and bulky.

## Choosing what to download

- **Re-use / re-analyze the atlases:** download `outputs/objects/` and toggle the
  relevant `steps:` flags.
- **Regenerate figures or tables:** download `outputs/objects/` and run the
  viz/analysis steps via `config/config.run.yml`.
- **Rebuild everything from raw counts:** download `inputs/data/runs/`, recreate
  `inputs/data/metadata.csv` from the manuscript supplement, then enable ingest
  and integration. UMAP embeddings are not seeded, so from-scratch coordinates
  differ from the published ones (cluster IDs are stable). See "Reproducibility"
  in the README.

## External references not redistributed here

- **IMGT germline reference** for BCR lineage reconstruction. Obtain from IMGT
  and point `paths$imgt_dir` in the config at it. Those steps skip cleanly if it
  is absent.

## Sample metadata

The de-identified sample metadata table (`inputs/data/metadata.csv`) is published
as a **supplementary table in the manuscript**, not in the Zenodo archive. To run
the pipeline, recreate `inputs/data/metadata.csv` from that supplementary table
and remap the `*_CRoutput` path columns to your local copy of the raw CellRanger
outputs. The column schema is documented in
[`data_dictionary.md`](data_dictionary.md), and a header-only example is provided
at `inputs/data/metadata_template.csv`.
