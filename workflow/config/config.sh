# config.sh — User-configurable paths for the Slide-tags pipeline
#
# This file is sourced by slidetag_pipe.sh. Edit the paths below to match
# your environment. See FASRC_SETUP.md for detailed setup instructions.

# ---------------------------------------------------------------------------
# Cluster and conda
# ---------------------------------------------------------------------------

# Path to cluster scheduler bin (leave blank for SLURM — auto-detected)
export CLUSTER_PATH=""

# Path to conda bin directory (system conda on FASRC)
export CONDA_PATH="/n/sw/Miniforge3-25.3.1-0/bin"

# Path to the slidetags conda environment
# Find yours with: conda env list | grep slidetags
export ENV_PATH="/n/home12/nhilgert/.conda/envs/slidetags"

# ---------------------------------------------------------------------------
# Pipeline paths
# ---------------------------------------------------------------------------

# Root of the Slide-tags repo (no trailing slash, no "workflow" suffix)
export WORKFLOW_PATH="/n/eddy_lab/users/nhilgert/Slide-tags"

# Where pipeline data lives. Each sample creates a subdirectory here:
#   <BASE_DATA_PATH>/<SAMPLE>/fastq/, RNAcounts/, spatial/, log/
export BASE_DATA_PATH="/n/eddy_lab/users/nhilgert/Slide-tags/data"

# Directory containing puck coordinate CSVs (filename = puck_id + .csv)
# Format: 3 columns, no header — barcode, x_um, y_um
export PUCK_PATH="/n/eddy_lab/users/nhilgert/Slide-tags/pucks"

# ---------------------------------------------------------------------------
# External software and references
# ---------------------------------------------------------------------------

# Parent directory containing cellranger-7.2.0/, cellranger-8.0.1/, etc.
# The pipeline appends the cellranger version directory automatically.
export PKG_PATH="/n/eddy_lab/users/nhilgert/software"

# Directory containing cellranger reference genomes:
#   refdata-gex-mm10-2020-A, refdata-gex-GRCh38-2024-A, etc.
# Species-to-reference mapping is in scripts/SubmitFuncs/make_input_csv.py
export REF_PATH="/n/eddy_lab/users/nhilgert/reference"

# ---------------------------------------------------------------------------
# Google Cloud / Sheets (optional — leave blank if using -csv mode)
# ---------------------------------------------------------------------------

# Path to BCL files (used as default input path for Google Sheets mode)
export BCL_MAIN_PATH=""

# Google Sheet ID for sample metadata (leave blank for -csv mode)
export GOOGLE_SHEET_ID=""

# Google Cloud Storage bucket for FASTQ/BAM upload (leave blank to skip)
export GOOGLE_CLOUD_BUCKET=""

# Path to slide-seq puck barcode file (leave blank if not used)
export PUCK_IN=""

# ---------------------------------------------------------------------------
# PATH setup (do not edit unless you know what you're doing)
# ---------------------------------------------------------------------------
export PATH="$CLUSTER_PATH:$PATH"
export PATH="$CONDA_PATH:$PATH"

# ---------------------------------------------------------------------------
# Parallelism — number of samples to process per SLURM job
# ---------------------------------------------------------------------------
# Increase these to submit fewer, larger SLURM jobs.
# Decrease to 1 for debugging (one job per sample).
export CHUNK_SIZE_MKFASTQ=1
export CHUNK_SIZE_RNACOUNTS=1
export CHUNK_SIZE_CELLBENDER=1
export CHUNK_SIZE_SBCOUNT=10
export CHUNK_SIZE_POSITION=10
