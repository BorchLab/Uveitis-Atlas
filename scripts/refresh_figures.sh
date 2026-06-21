#!/usr/bin/env bash
# scripts/refresh_figures.sh
# ---------------------------------------------------------------------------
# Safely re-render all figures from the EXISTING objects and tables, using
# config/config.run.yml (visualization steps only; every atlas-builder and
# compute step is off). This never rebuilds or overwrites an atlas, and is
# deterministic (no UMAP recompute).
#
# What it does:
#   1. Backs up outputs/objects and outputs/viz to timestamped copies.
#   2. Clears outputs/viz so figures removed from the code do not linger.
#   3. Runs the pipeline with the viz-only run config.
#
# Run from the repository root:
#   bash scripts/refresh_figures.sh
# ---------------------------------------------------------------------------
set -euo pipefail

# Guard: must be at repo root (where R/ and config/ live).
if [[ ! -d R || ! -f config/config.run.yml ]]; then
  echo "ERROR: run from the repository root (need R/ and config/config.run.yml)." >&2
  exit 1
fi

stamp="$(date +%Y%m%d_%H%M%S)"

echo "==> Backing up objects -> outputs/objects.backup_${stamp}"
cp -a outputs/objects "outputs/objects.backup_${stamp}"

if [[ -d outputs/viz ]]; then
  echo "==> Backing up figures -> outputs/viz.backup_${stamp}"
  cp -a outputs/viz "outputs/viz.backup_${stamp}"
  echo "==> Clearing outputs/viz"
  rm -rf outputs/viz
fi
mkdir -p outputs/viz

echo "==> Running viz-only pipeline (config/config.run.yml)"
Rscript run_pipeline.R config/config.run.yml

echo "==> Done. New figures in outputs/viz/"
echo "    Object backup:  outputs/objects.backup_${stamp}"
echo "    Figure backup:  outputs/viz.backup_${stamp}"
echo "    Compare with:   diff -rq outputs/viz \"outputs/viz.backup_${stamp}\""
