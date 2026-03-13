#!/bin/bash
# batch_submit_sbcounts.sh
# Submits SBcounts pipeline for all samples in decidualization_all.csv.
# Each sample is submitted as an independent pipeline invocation.
#
# Usage:
#   bash batch_submit_sbcounts.sh             # submit all samples
#   bash batch_submit_sbcounts.sh mck_8a      # submit one sample by name
#   bash batch_submit_sbcounts.sh --dry-run   # print commands without submitting
#
# Run from: /n/eddy_lab/users/nhilgert/Slide-tags

set -euo pipefail

WORKFLOW_PATH="/n/eddy_lab/users/nhilgert/Slide-tags"
SAMPLES_CSV="$WORKFLOW_PATH/samples/decidualization_all.csv"
DRY_RUN=false
FILTER_NAME=""

# Parse args
for arg in "$@"; do
    case $arg in
        --dry-run) DRY_RUN=true ;;
        *) FILTER_NAME="$arg" ;;
    esac
done

cd "$WORKFLOW_PATH"

echo "========================================"
echo "  Slide-tags Batch SBcounts Submission  "
echo "========================================"
echo ""
[ "$DRY_RUN" = true ] && echo "  *** DRY RUN — no jobs will be submitted ***"
echo ""

submitted=0
skipped=0

tail -n +2 "$SAMPLES_CSV" | while IFS=',' read -r name rna_index sb_index puck_id species rest; do

    # Filter to single sample if specified
    if [ -n "$FILTER_NAME" ] && [ "$name" != "$FILTER_NAME" ]; then
        continue
    fi

    # Check fastq dirs exist
    gex_dir="$WORKFLOW_PATH/data/$name/fastq/$rna_index"
    sb_dir="$WORKFLOW_PATH/data/$name/fastq/$sb_index"

    if [ ! -d "$gex_dir" ] || [ ! -d "$sb_dir" ]; then
        echo "  SKIP $name — fastq dirs not set up (run setup_fastq_dirs.sh first)"
        ((skipped++)) || true
        continue
    fi

    # Check if SBcounts.h5 already exists (already done)
    h5_path="$WORKFLOW_PATH/data/$name/spatial/$rna_index/SBcounts/SBcounts.h5"
    if [ -f "$h5_path" ]; then
        echo "  SKIP $name — SBcounts.h5 already exists"
        ((skipped++)) || true
        continue
    fi

    # Write per-sample CSV
    sample_csv="$WORKFLOW_PATH/samples/${name}.csv"
    echo "name,rna_index,sb_index,puck_id,species,chemistry,vdj,cellranger,store_path" > "$sample_csv"
    echo "${name},${rna_index},${sb_index},${puck_id},${species},,,," >> "$sample_csv"

    cmd="bash workflow/slidetag_pipe.sh -bcl $name -sb -csv $sample_csv"

    if [ "$DRY_RUN" = true ]; then
        echo "  WOULD RUN: $cmd"
    else
        echo "  Submitting $name ($species, puck=$puck_id)..."
        bash workflow/slidetag_pipe.sh -bcl "$name" -sb -csv "$sample_csv" &
        sleep 2  # avoid hammering snakemake/SLURM simultaneously
    fi
    ((submitted++)) || true
done

echo ""
echo "Submitted: $submitted  Skipped: $skipped"
echo ""
echo "Monitor with: squeue -u nhilgert"
echo "Check logs:   tail data/<sample>/log/<run>/spatial_logs/SBcount_*.log"
