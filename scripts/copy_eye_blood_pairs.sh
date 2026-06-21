#!/bin/bash

# =========================
# Configuration Parameters
# =========================

# Mounted source directory (already available locally)
SRC_BASE="/Volumes/SeqData/Eye_Blood_Pairs/AB_cellranger_v8.0.1/cellranger_v8.0.1"

# Local destination directory (change as desired)
DEST_BASE="$HOME/Eye_Blood_Pairs_local"

# =========================
# Execution
# =========================

echo "Starting local copy..."
mkdir -p "$DEST_BASE"

# Loop through each sample directory
for sample_dir in "$SRC_BASE"/*/; do
    sample_name=$(basename "$sample_dir")
    echo "Processing sample: $sample_name"

    # Define the relevant subdirectories
    COUNT_DIR="${sample_dir}outs/outs/per_sample_outs/${sample_name}/count/sample_filtered_feature_bc_matrix"
    VDJ_T_DIR="${sample_dir}outs/outs/per_sample_outs/${sample_name}/vdj_t"
    VDJ_B_DIR="${sample_dir}outs/outs/per_sample_outs/${sample_name}/vdj_b"

    # Create destination folder
    mkdir -p "${DEST_BASE}/${sample_name}"

    # Copy each if it exists
    for SRC in "$COUNT_DIR" "$VDJ_T_DIR" "$VDJ_B_DIR"; do
        if [ -d "$SRC" ]; then
            echo "  Copying $(basename "$SRC")..."
            cp -a "$SRC" "${DEST_BASE}/${sample_name}/"
        else
            echo "  ⚠️  Missing directory: $SRC (skipped)"
        fi
    done
done

echo "✅ Copy completed. Files saved under: $DEST_BASE"
