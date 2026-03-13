# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Slide-tags is a **Snakemake-based spatial transcriptomics pipeline** that processes raw BCL/FASTQ sequencing data into cell-resolved spatial coordinates. It was originally built for SGE/UGER and has been ported to **SLURM** for Harvard's FASRC Eddy partition.

Pipeline stages: `mkfastq → RNAcounts (cellranger) → [CellBender] → SBcounts (Julia) → Spatial positioning (R)`

## Running the Pipeline

```bash
# Activate the conda environment
conda activate slidetags

# Full pipeline run (cellranger + spatial)
bash workflow/slidetag_pipe.sh -bcl <BCL_NAME> -cr -sp -csv samples/<input>.csv

# SBcounts only (skip cellranger)
bash workflow/slidetag_pipe.sh -bcl <BCL_NAME> -sb -csv samples/<input>.csv

# Spatial positioning only (reuse existing input sheets)
bash workflow/slidetag_pipe.sh -bcl <BCL_NAME> -sp -us

# Batch convenience script for SBcounts across all samples
bash batch_submit_sbcounts.sh
```

Key flags: `-bcl` (BCL name, required), `-cr`/`-sb`/`-sp` (stages), `-csv` (local CSV input), `-us` (reuse input sheets), `-f` (force re-run), `-ec`/`-td` (CellBender params).

## Architecture

### Entry Point & Job Orchestration
- **`workflow/slidetag_pipe.sh`** — Main entry point. Parses flags, sources config, calls submit functions for each stage.
- **`workflow/scripts/SubmitFuncs/submit_functions.sh`** — SLURM job submission logic. Chunks samples into batches, submits via `sbatch`, monitors with `squeue`. Key functions: `submit_batch_jobs()`, `make_*_input_list()`, `submit_*_job()`.
- **`workflow/scripts/SubmitFuncs/make_input_csv.py`** — Converts simple 9-column sample CSVs into per-sample input sheets the pipeline consumes (replaces Google Sheets dependency).

### Snakemake Rules (`workflow/rules/`)
Each `.smk` file is a self-contained stage invoked independently by `submit_functions.sh` (not a single DAG):
- `read_check.smk` — Input validation
- `mkfastq_split.smk` — BCL demultiplexing via cellranger mkfastq
- `RNAcounts_split.smk` — Gene expression quantification via cellranger count/arc
- `Cellbender_split.smk` — Ambient RNA removal (optional)
- `SBcounts_split.smk` — Spatial barcode counting (Julia)
- `Spatial_split.smk` — Cell positioning via DBSCAN + KDE (R)
- `move_data.smk` — Archive results to final storage

Resource allocations per rule are in **`workflow/config/smk_config.yaml`**.

### Core Processing Scripts
- **Julia** (`workflow/scripts/SpatialCount/spatial_count.jl`) — Reads spatial FASTQ, extracts barcodes via linker sequence (`TCTTCAGCGTTCCCGAGA`), matches to puck coordinates, outputs HDF5. Uses Julia multithreading.
- **R** (`workflow/scripts/Positioning/positioning.R`) — Loads SBcounts.h5, runs DBSCAN hyperparameter search (eps, minPts) + KDE, outputs cell coordinates as `seurat.qs`.
- **Shell** (`workflow/scripts/Cellranger/`) — Wrappers for cellranger count/arc, FASTQ merging, reference genome building.

## Configuration

- **`workflow/config/config.sh`** — All user-configurable paths: conda, cellranger, data, pucks, references, chunk sizes for parallelization.
- **`workflow/config/main_env.yml`** — Conda env spec (Python 3.9, R 4.3.3, Julia 1.10.3, Snakemake 7.32.4).
- **`workflow/config/cellbender_env.yml`** — Separate conda env for CellBender.
- **`workflow/config/smk_config.yaml`** — Per-rule SLURM resources (mem, CPUs, runtime).

## Data Layout

```
data/<BCL_NAME>/
├── fastq/<INDEX>/          # Demultiplexed FASTQs
├── RNAcounts/<SAMPLE>/outs/  # Cellranger output (filtered_feature_bc_matrix.h5)
├── spatial/<SAMPLE>/
│   ├── SBcounts/SBcounts.h5  # Spatial barcode counts
│   ├── Positions/seurat.qs   # Cell coordinates
│   └── puck/<PUCK>.csv       # Puck bead positions
└── log/<RUN>/              # Logs and generated input sheets
```

Input sample CSVs go in `samples/`. Puck coordinate CSVs go in `pucks/`.

## Environment Setup

See `FASRC_SETUP.md` for detailed FASRC setup instructions, including:
- Conda environment creation
- Julia HDF5/OpenSSL artifact symlink fix
- Reference genome paths

Julia packages use a separate project at `$JULIA_PROJECT_PATH` (`~/.julia/environments/slidetags`), not the conda environment's Julia package directory.

## Key Conventions

- Each pipeline stage runs as independent Snakemake invocations (not one unified DAG). Stages are orchestrated by shell functions in `submit_functions.sh`.
- Sample parallelization is controlled by `CHUNK_SIZE_*` variables in `config.sh` — samples are batched into chunks and submitted as separate SLURM jobs.
- The `-csv` flag enables offline mode (local CSV input); without it, the pipeline reads from Google Sheets via `read_info.py`.
- Species-to-reference mapping is handled in `make_input_csv.py` (e.g., "Mouse" → `refdata-gex-mm10-2020-A`).

## Known Issues and Hard-Won Fixes

**DO NOT revert any of these — they were discovered through extensive debugging.**

### 1. `Pkg.precompile()` removed from `spatial_count.jl`
`workflow/scripts/SpatialCount/spatial_count.jl` line 4 previously had `Pkg.precompile()`.
It was removed. DO NOT add it back. It triggers fresh package downloads and OpenSSL conflicts
on every SLURM job, causing the Julia process to hang or crash silently.
The script should only call `Pkg.activate(project_path)`.

### 2. HDF5/OpenSSL symlink fix
Julia's HDF5_jll bundles OpenSSL 3.4.0 but the conda env has 3.3.x. This causes a library
conflict that kills Julia silently. Fix: symlink the HDF5_jll artifact's libssl/libcrypto to
the conda env's versions. This must be applied to TWO locations:
- `~/.conda/envs/slidetags/share/julia/artifacts/1c2814ca8577e0356a6ecf8447364c19e97b86a8/lib/`
- `~/.julia/artifacts/1c2814ca8577e0356a6ecf8447364c19e97b86a8/lib/`

If Julia jobs start failing with OpenSSL errors, reapply with:
```bash
ARTIFACT_DIR="<path_above>/lib"
CONDA_LIB="/n/home12/nhilgert/.conda/envs/slidetags/lib"
ln -sf $CONDA_LIB/libssl.so.3 $ARTIFACT_DIR/libssl.so
ln -sf $CONDA_LIB/libssl.so.3 $ARTIFACT_DIR/libssl.so.3
ln -sf $CONDA_LIB/libcrypto.so.3 $ARTIFACT_DIR/libcrypto.so
ln -sf $CONDA_LIB/libcrypto.so.3 $ARTIFACT_DIR/libcrypto.so.3
```

### 3. `tee` replaced with `&>>` in `SBcounts_split.smk`
All `| tee -a {log.SBcounts_log}` were replaced with `&>> {log.SBcounts_log}`.
`tee` causes silent failures under `set -euo pipefail` (pipefail treats a broken pipe as
an error). The symptom is a log file containing only "FAILURE" with no preceding output,
and the job completing in under 3 minutes. DO NOT reintroduce `tee` here.

### 4. `puck_in` in `params:` not `input:` in `SBcounts_split.smk`
`puck_in` is an optional empty string when no puck override is provided. Snakemake does not
allow empty strings as `input:` entries — it will fail at DAG construction. It must stay
in `params:` and referenced as `{params.puck_in}` in the shell block.

### 5. JULIA env vars hardcoded in `SBcounts_split.smk`
These must be hardcoded in the shell block — do not rely on them being inherited from the
environment, as SLURM compute nodes do not inherit the submitting shell's exports:
```bash
export JULIA_PROJECT_PATH="/n/home12/nhilgert/.julia/environments/slidetags"
export JULIA_DEPOT_PATH="/n/home12/nhilgert/.julia"
export JULIA_PACKAGES_PATH="/n/home12/nhilgert/.julia/packages"
```

### 6. `lineterminator="\n"` in `make_input_csv.py`
The CSV writer uses `lineterminator="\n"` explicitly. DO NOT change this to the default.
The default produces Windows-style `\r\n` line endings. When bash reads the CSV with
`while IFS=',' read`, the `\r` becomes part of the last field, causing `$idx` to be
e.g. `SI-NT-F11\r` which never matches the wildcard and causes silent FAILURE.

### 7. Diagnosing silent FAILURE
If a sample log contains only "FAILURE" with no preceding output:
1. Check `data/<sample>/log/<Nrun>/main/SBcounts_all.err` for snakemake stderr
2. If the job completed in <3 minutes, Julia never ran — the shell died before reaching it
3. Common cause: `set -euo pipefail` + a failing command before the Julia invocation
4. The puck copy step (`cp "$puck_csv" "$working_dir"`) runs before Julia — if the puck
   file doesn't exist, `exit 0` is called (check pucks/ directory)
5. Check if `data/<sample>/spatial/<gex_index>/SBcounts/` directory was created —
   if it exists but no .h5, Julia ran and failed; if it doesn't exist, died before Julia

### 8. CHOLMOD warning is harmless
Julia logs will always contain a CHOLMOD version incompatibility warning. This is expected
and harmless. A job is successful if and only if the log ends with "SUCCESS".

### 9. Pipeline uses `-bcl` as folder name only
We are NOT running BCL demultiplexing. The `-bcl` flag is repurposed as a sample folder
name under `data/`. FASTQs are already symlinked under `data/<sample>/fastq/<index>/`.
Always use `-csv samples/<sample>.csv` to provide input — do not rely on Google Sheets.

### 10. Cellranger outputs must be symlinked for `-sp`
Spatial positioning (`-sp`) requires cellranger output at `data/<sample>/RNAcounts/<sample>/`.
Existing cellranger outputs live at:
`/n/eddy_lab/users/nhilgert/decidualization/revisions/velocity/<sample>/cellranger/`
Symlink before running `-sp`:
```bash
mkdir -p data/<sample>/RNAcounts
ln -s /n/eddy_lab/users/nhilgert/decidualization/revisions/velocity/<sample>/cellranger \
    data/<sample>/RNAcounts/<sample>
```
Not all samples have cellranger output at this path — verify before symlinking.
## FASTQ Source Layout
FASTQs are symlinked from lanes into the pipeline's data directory via setup_fastq_dirs.sh.
Source structure:
- `/n/eddy_lab/users/nhilgert/decidualization/lanes/<sample>/gex/` → GEX FASTQs
- `/n/eddy_lab/users/nhilgert/decidualization/lanes/<sample>/spatial/` → Spatial FASTQs

These in turn point to:
- `/n/eddy_lab/users/nhilgert/decidualization/fastqs/mck_1/` (mck_1 only)
- `/n/eddy_lab/users/nhilgert/decidualization/fastqs/AE/` (mck_3-7b)
- `/n/eddy_lab/users/nhilgert/decidualization/fastqs/FG/` (mck_8a-15)
- `/n/eddy_lab/users/nhilgert/decidualization/fastqs/E/` (mck_16-20, no pucks yet)

Do NOT modify anything under `/n/eddy_lab/users/nhilgert/decidualization/` — source data only.
The pipeline reads FASTQs from `data/<sample>/fastq/<index>/` which are symlinks set up by setup_fastq_dirs.sh.

## Storage Layout
- **Home directory** (`/n/home12/nhilgert/`) — config files ONLY (.bashrc, .claude, .julia, .conda).
  Do NOT write data, logs, or pipeline outputs here.
- **Working directory** (`/n/eddy_lab/users/nhilgert/`) — all pipeline data, logs, and outputs go here.
- **Slide-tags repo** (`/n/eddy_lab/users/nhilgert/Slide-tags/`) — pipeline code and data.
- **FASTQs** (`/n/eddy_lab/users/nhilgert/decidualization/`) — source data, do not modify.
- Never write pipeline outputs, logs, or large files to `/n/home12/nhilgert/`.

## Shell Environment
User's shell config is at `/n/home12/nhilgert/.bashrc`. Read this file if you need to
understand PATH, conda setup, or environment variables. Key exports already set there:
- JULIA_PROJECT_PATH=/n/home12/nhilgert/.julia/environments/slidetags
- conda is at /n/sw/Miniforge3-25.3.1-0/bin

## Guardrails
- Never call sbatch directly — always use batch_submit_sbcounts.sh or slidetag_pipe.sh
- Never have more than 50 jobs in the queue at once — check squeue before submitting
- Always check squeue before submitting to avoid duplicate jobs for the same sample
- Do not modify workflow/config/smk_config.yaml resource limits
- Do not run cellranger (-cr flag) — outputs already exist and are symlinked
- Do not run Cellbender (-cb flag) — not configured
- Do not poll squeue more than once per minute — FASRC policy
- Do not run heavy computation directly in the terminal — submit via SLURM
- Do not write large files to home directory (/n/home12/nhilgert/)
- If a fix attempt results in the same error recurring, stop and document in claude_status.log rather than retrying indefinitely
- SLURM username is `nhilgert` — always use `squeue -u nhilgert`

## FASRC Cluster Etiquette
- Login nodes are shared — no heavy CPU/memory usage directly on them
- Always submit compute work via SLURM (sbatch/salloc)
- Keep files under 1000 per directory to avoid filesystem latency issues
- Do not use home directory for large data — use /n/eddy_lab/
- Scratch (/n/netscratch) is temporary — files deleted after 90 days, not backed up
- Job time limits: always request realistic walltime, don't over-request memory