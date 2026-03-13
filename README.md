
# Slide-tags Pipeline

A Snakemake-based pipeline for [Slide-tags](https://www.nature.com/articles/s41586-023-06837-4)
spatial transcriptomics. Processes raw FASTQs into cell-resolved spatial coordinates.
Runs on FASRC using the Eddy SLURM partition.

Ported from [thechenlab/Slide-tags](https://github.com/thechenlab/Slide-tags) (SGE/UGER).

<img width="1594" height="1328" alt="image" src="https://github.com/user-attachments/assets/2c842564-99b7-466d-9cf1-4ed5f5a92f4f" />

## What it does

```
FASTQs ──► SBcounts (Julia) ──► Spatial Positioning (R) ──► seurat.qs
                                     │
            Cellranger output ───────┘
```

| Stage | Tool | Input | Output |
|---|---|---|---|
| **SBcounts** | Julia | Spatial barcode FASTQs + puck coordinates | `SBcounts.h5` |
| **Spatial Positioning** | R (DBSCAN + KDE) | `SBcounts.h5` + cellranger `filtered_feature_bc_matrix.h5` | `seurat.qs` |

The pipeline also supports upstream stages (BCL demux via `mkfastq`, gene expression
quantification via `cellranger count`), but most users will start from existing
cellranger output.

## Prerequisites

1. **FASRC account** with access to the Eddy partition
2. **Conda environment** (`slidetags`) — see [FASRC_SETUP.md](FASRC_SETUP.md) for one-time setup
3. **Cellranger output** for each sample (existing runs can be symlinked)
4. **Puck coordinate CSVs** in `pucks/` (3 columns, no header: barcode, x, y)
5. **FASTQs** organized under `data/<SAMPLE>/fastq/<INDEX>/`

## Quick start

### 1. One-time setup

Follow [FASRC_SETUP.md](FASRC_SETUP.md) to:
- Create the `slidetags` conda environment
- Install Julia packages and apply the HDF5/OpenSSL fix
- Fill in `workflow/config/config.sh` with your paths

### 2. Prepare a sample CSV

Create a file in `samples/` with one row per sample:

```csv
name,rna_index,sb_index,puck_id,species,chemistry,vdj,cellranger,store_path
my_sample,SI-TT-A5,SI-NT-A10,Puck_230525_28,Mouse,,,,
```

| Field | Description |
|---|---|
| `name` | Sample name (used as directory name under `data/`) |
| `rna_index` | GEX FASTQ folder name (10x index, e.g. `SI-TT-A5`) |
| `sb_index` | Spatial barcode FASTQ folder name |
| `puck_id` | Puck filename without `.csv` (must exist in your `PUCK_PATH`) |
| `species` | `Mouse` or `Human` |
| Other fields | Leave blank unless needed |

### 3. Set up FASTQ and cellranger symlinks

FASTQs must be at `data/<SAMPLE>/fastq/<INDEX>/`:

```bash
# Create FASTQ directories
mkdir -p data/my_sample/fastq/SI-TT-A5 data/my_sample/fastq/SI-NT-A10

# Symlink your FASTQs
ln -s /path/to/gex/fastqs/*     data/my_sample/fastq/SI-TT-A5/
ln -s /path/to/spatial/fastqs/* data/my_sample/fastq/SI-NT-A10/
```

If cellranger has already been run, symlink the output:

```bash
mkdir -p data/my_sample/RNAcounts
ln -s /path/to/cellranger/output data/my_sample/RNAcounts/my_sample
```

The cellranger directory must contain `outs/filtered_feature_bc_matrix.h5` and
`outs/molecule_info.h5`.

### 4. Run the pipeline

**Single sample (recommended for first run):**

```bash
conda activate slidetags
bash run_sample.sh my_sample
```

This runs SBcounts + Spatial Positioning, submits SLURM jobs, and monitors them
to completion. See `run_sample.sh -h` for options.

**Using the pipeline directly:**

```bash
# SBcounts + Spatial positioning
bash workflow/slidetag_pipe.sh -bcl my_sample -sb -sp -csv samples/my_sample.csv

# Spatial positioning only (if SBcounts already done)
bash workflow/slidetag_pipe.sh -bcl my_sample -sp -csv samples/my_sample.csv

# Rerun using existing input sheets (no CSV needed)
bash workflow/slidetag_pipe.sh -bcl my_sample -sp -us
```

**Batch run:**

```bash
for s in mck_1 mck_3 mck_4; do
    bash run_sample.sh "$s" &
    sleep 10  # stagger submissions to avoid overwhelming SLURM
done
wait
```

### 5. Check outputs

```bash
# Quick check: does seurat.qs exist?
ls data/my_sample/spatial/*/Positions/seurat.qs

# Check the log (should end with "SUCCESS")
tail -5 data/my_sample/log/*/spatial_logs/Spatial_*.log
```

## Output structure

```
data/<SAMPLE>/
├── spatial/<RNA_INDEX>/
│   ├── SBcounts/
│   │   └── SBcounts.h5              # Spatial barcode count matrix
│   └── Positions/
│       ├── seurat.qs                 # Seurat object with spatial coordinates
│       ├── coords.csv                # Cell coordinate table
│       ├── summary.pdf               # QC plots (may be partial, see note below)
│       └── cb_whitelist.txt          # Cell barcode whitelist
└── log/<Nrun>/
    ├── main/spatial.{log,err}        # Snakemake stdout/stderr
    ├── spatial_logs/Spatial_*.log    # Per-sample log (ends SUCCESS or FAILURE)
    └── input/                        # Generated input sheets
```

**Loading the output in R:**

```r
library(qs)
obj <- qread("data/my_sample/spatial/SI-TT-A5/Positions/seurat.qs")

# Spatial coordinates are in reductions:
#   obj[["spatial"]]  — KDE-filtered DBSCAN positions (recommended)
#   obj[["dbscan"]]   — Raw DBSCAN positions
#   obj[["kde"]]      — KDE positions
# Also in metadata: obj$x_um, obj$y_um
```

## Configuration

All user-configurable paths are in **`workflow/config/config.sh`**:

| Variable | Description |
|---|---|
| `CONDA_PATH` | Path to conda `bin/` (system conda on FASRC: `/n/sw/Miniforge3-25.3.1-0/bin`) |
| `ENV_PATH` | Path to your `slidetags` conda environment |
| `PKG_PATH` | Parent dir containing `cellranger-7.2.0/` |
| `BASE_DATA_PATH` | Where pipeline `data/` lives |
| `WORKFLOW_PATH` | Path to this repo root (no trailing slash) |
| `PUCK_PATH` | Directory containing puck coordinate CSVs |
| `REF_PATH` | Directory containing cellranger reference genomes |

SLURM resource allocations are in **`workflow/config/smk_config.yaml`**.
Do not modify these unless you understand the implications.

Parallelism is controlled by `CHUNK_SIZE_*` variables in `config.sh`:
```bash
CHUNK_SIZE_SBCOUNT=10    # How many SBcounts samples per SLURM job
CHUNK_SIZE_POSITION=10   # How many Spatial samples per SLURM job
```

## Pipeline flags reference

| Flag | Description |
|---|---|
| `-bcl <NAME>` | **Required.** Sample folder name under `data/` |
| `-sb` | Run SBcounts (Julia spatial barcode counting) |
| `-sp` | Run Spatial Positioning (R DBSCAN + KDE) |
| `-cr` | Run cellranger count |
| `-csv <FILE>` | Use a local CSV for sample metadata |
| `-us` | Reuse input sheets from the last run |
| `-f` | Force re-run (deletes existing outputs first) |

Flags can be combined: `-sb -sp` runs SBcounts then Spatial.

## Troubleshooting

| Symptom | Likely cause | Fix |
|---|---|---|
| Log shows only "FAILURE" in < 3 min | Shell died before R/Julia | Check `data/<s>/log/<N>/main/*.err` for the real error |
| `oom_kill` in `.err` file | Sample too large for memory | Increase `--mem` in `submit_functions.sh` (currently 64G) |
| "Required files are missing: SBcounts.h5" but file exists | Corrupt HDF5 | Run `h5ls` on the file; if empty, delete and rerun `-sb` |
| UMAP/plot warnings in log | rlang 1.1.3 S4SXP bug | Harmless — plots skipped but `seurat.qs` is fine. Update rlang to >= 1.1.4 to fix |
| CHOLMOD version warning | Julia library mismatch | Harmless, ignore. Job is successful iff log ends with "SUCCESS" |

## Project structure

```
workflow/
├── slidetag_pipe.sh                  # Main entry point
├── config/
│   ├── config.sh                     # User paths — edit this
│   ├── main_env.yml                  # Conda env spec
│   └── smk_config.yaml              # SLURM resources per rule
├── rules/                            # Snakemake rules (one file per stage)
│   ├── SBcounts_split.smk
│   ├── Spatial_split.smk
│   ├── RNAcounts_split.smk
│   └── ...
└── scripts/
    ├── SpatialCount/spatial_count.jl # Julia barcode extraction + counting
    ├── Positioning/                   # R spatial positioning
    │   ├── run_spatial.R             # Entry point (called by Snakemake)
    │   ├── run_positioning.R         # Orchestrator (RNA load → positioning → save)
    │   ├── load_matrix.R            # Spatial matrix loading + QC
    │   ├── positioning.R            # DBSCAN + KDE cell positioning
    │   └── Functions/               # Shared R helper functions
    └── SubmitFuncs/
        ├── submit_functions.sh       # SLURM job submission logic
        └── make_input_csv.py         # CSV → per-sample input sheets
```
