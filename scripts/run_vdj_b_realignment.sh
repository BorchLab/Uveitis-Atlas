#!/bin/bash

# =========================
# Configuration Parameters
# =========================

# Base directory (where your sample folders live)
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
BASE_DIR="${SCRIPT_DIR}/../inputs/data/runs"

# Path to IgBLAST database and germline reference
IGBLAST_DB="$HOME/share/igblast"
GERMLINE_DIR="$HOME/share/germlines/imgt/human/vdj"

# Path to 10x annotations CSV (assumed to be next to FASTA)
ANNOTATION_FILE="filtered_contig_annotations.csv"

# =========================
# Environment setup
# =========================

# Activate Change-O virtual environment
source "$HOME/BCR/bin/activate"

# =========================
# Realignment
# =========================

echo "Searching for vdj_b directories under: $BASE_DIR"
find "$BASE_DIR" -type d -name "vdj_b" | while read -r vdj_b_dir; do
    fasta_file="${vdj_b_dir}/filtered_contig.fasta"

    if [ -f "$fasta_file" ]; then
        echo "============================================"
        echo "Processing: $fasta_file"
        echo "============================================"

        # Run AssignGenes.py
        AssignGenes.py igblast \
            -s "$fasta_file" \
            -b "$IGBLAST_DB" \
            --organism human \
            --loci ig \
            --format blast

        # Construct paths for MakeDb.py
        fmt7_file="${vdj_b_dir}/filtered_contig_igblast.fmt7"
        annotation_file="${vdj_b_dir}/${ANNOTATION_FILE}"

        if [ -f "$fmt7_file" ] && [ -f "$annotation_file" ]; then
            MakeDb.py igblast \
                -i "$fmt7_file" \
                -s "$fasta_file" \
                -r "$GERMLINE_DIR" \
                --10x "$annotation_file" \
                --extended
        else
            echo "⚠️ Missing igblast output or annotation file in: $vdj_b_dir"
        fi
    else
        echo "⚠️ No FASTA file found in: $vdj_b_dir"
    fi
done

echo "✅ Realignment completed."
