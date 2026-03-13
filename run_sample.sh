#!/bin/bash
# run_sample.sh — Run the Slide-tags pipeline for a single sample
#
# Usage:
#   bash run_sample.sh <SAMPLE_NAME> [OPTIONS]
#
# Examples:
#   bash run_sample.sh mck_1              # SBcounts + Spatial
#   bash run_sample.sh mck_1 -sp          # Spatial only
#   bash run_sample.sh mck_1 -sb          # SBcounts only
#   bash run_sample.sh mck_1 -f           # Force re-run
#
# Prerequisites:
#   - conda activate slidetags
#   - samples/<SAMPLE_NAME>.csv exists
#   - FASTQs symlinked under data/<SAMPLE_NAME>/fastq/
#   - Cellranger output symlinked under data/<SAMPLE_NAME>/RNAcounts/ (for -sp)

set -euo pipefail
cd "$(dirname "$0")"

# --- Parse arguments ---
if [[ $# -lt 1 ]] || [[ "$1" == "-h" ]] || [[ "$1" == "--help" ]]; then
    echo "Usage: bash run_sample.sh <SAMPLE_NAME> [OPTIONS]"
    echo ""
    echo "Runs SBcounts + Spatial Positioning for a single sample."
    echo ""
    echo "Arguments:"
    echo "  SAMPLE_NAME    Name of the sample (must match samples/<name>.csv)"
    echo ""
    echo "Options:"
    echo "  -sb            Run SBcounts only (skip spatial)"
    echo "  -sp            Run Spatial positioning only (skip SBcounts)"
    echo "  -f             Force re-run (delete existing outputs first)"
    echo "  -h, --help     Show this help message"
    echo ""
    echo "If no stage flag is given, runs both SBcounts and Spatial."
    echo ""
    echo "Examples:"
    echo "  bash run_sample.sh mck_1          # Full pipeline"
    echo "  bash run_sample.sh mck_1 -sp      # Spatial only"
    echo "  bash run_sample.sh mck_1 -sb -f   # Force rerun SBcounts"
    exit 0
fi

SAMPLE="$1"
shift

# Default: run both stages
RUN_SB=false
RUN_SP=false
FORCE=""
EXPLICIT_STAGE=false

while [[ $# -gt 0 ]]; do
    case "$1" in
        -sb) RUN_SB=true; EXPLICIT_STAGE=true; shift ;;
        -sp) RUN_SP=true; EXPLICIT_STAGE=true; shift ;;
        -f)  FORCE="-f"; shift ;;
        *)   echo "Unknown option: $1"; exit 1 ;;
    esac
done

# If no explicit stage, run both
if [[ "$EXPLICIT_STAGE" == false ]]; then
    RUN_SB=true
    RUN_SP=true
fi

# --- Validate inputs ---
CSV="samples/${SAMPLE}.csv"
if [[ ! -f "$CSV" ]]; then
    echo "ERROR: Sample CSV not found: $CSV"
    echo "Create it with: name,rna_index,sb_index,puck_id,species,chemistry,vdj,cellranger,store_path"
    exit 1
fi

DATA_DIR="data/${SAMPLE}"
if [[ ! -d "$DATA_DIR" ]]; then
    echo "ERROR: Data directory not found: $DATA_DIR"
    echo "Create it and symlink FASTQs under $DATA_DIR/fastq/<INDEX>/"
    exit 1
fi

# Build stage flags
STAGES=""
if [[ "$RUN_SB" == true ]]; then STAGES="$STAGES -sb"; fi
if [[ "$RUN_SP" == true ]]; then STAGES="$STAGES -sp"; fi

echo "=== Slide-tags: $SAMPLE ==="
echo "Stages:$STAGES ${FORCE:+(force)}"
echo "CSV:    $CSV"
echo "Data:   $DATA_DIR"
echo ""

# --- Run ---
bash workflow/slidetag_pipe.sh -bcl "$SAMPLE" $STAGES -csv "$CSV" $FORCE

# --- Check result ---
echo ""
echo "=== Results ==="
SEURAT=$(find "$DATA_DIR/spatial/" -name "seurat.qs" 2>/dev/null | head -1)
SBCOUNTS=$(find "$DATA_DIR/spatial/" -name "SBcounts.h5" 2>/dev/null | head -1)

if [[ "$RUN_SB" == true ]]; then
    if [[ -n "$SBCOUNTS" ]]; then
        echo "SBcounts.h5: OK ($(du -h "$SBCOUNTS" | cut -f1))"
    else
        echo "SBcounts.h5: MISSING"
    fi
fi

if [[ "$RUN_SP" == true ]]; then
    if [[ -n "$SEURAT" ]]; then
        echo "seurat.qs:   OK ($(du -h "$SEURAT" | cut -f1))"
    else
        echo "seurat.qs:   MISSING"
        echo ""
        echo "Check logs:"
        echo "  tail data/$SAMPLE/log/*/spatial_logs/Spatial_*.log"
        echo "  cat  data/$SAMPLE/log/*/main/spatial.err"
        exit 1
    fi
fi

echo ""
echo "Done."
