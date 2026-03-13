# Slide-tags FASRC Setup Guide
## SLURM port for the Eddy lab cluster

This documents how to set up and run the Slide-tags pipeline on FASRC (Harvard RC)
using the Eddy partition. This is a port of the original SGE/UGER pipeline from the
[thechenlab/Slide-tags](https://github.com/thechenlab/Slide-tags) repository.

---

## Prerequisites

- Access to the Eddy partition on FASRC
- FASTQs already demultiplexed and on cluster storage
- Puck coordinate CSVs (3 columns, no header: barcode, x, y)
- Cellranger reference genomes on cluster storage

---

## Directory structure

```
<WORKFLOW_PATH>/
├── workflow/
│   ├── config/
│   │   ├── config.sh          # ← fill in your paths
│   │   ├── main_env.yml       # conda env spec
│   │   └── smk_config.yaml    # per-rule SLURM resources
│   ├── rules/                 # Snakemake rules
│   ├── scripts/
│   └── slidetag_pipe.sh       # pipeline entry point
├── data/
│   └── <BCL_NAME>/
│       ├── fastq/
│       │   ├── <RNA_INDEX>/   # GEX FASTQs
│       │   └── <SB_INDEX>/    # Spatial barcode FASTQs
│       └── RNAcounts/
│           └── <SAMPLE>/      # cellranger output (can symlink existing)
├── pucks/                     # puck coordinate CSVs
└── samples/                   # input CSVs for pipeline runs
```

---

## One-time setup

### 1. Clone the repo

```bash
git clone https://github.com/nhilgert/Slide-tags.git
cd Slide-tags
```

### 2. Build the conda environment

Do this on an **interactive compute node**, not the login node:

```bash
srun -p eddy --mem=16G --cpus-per-task=4 --time=02:00:00 --pty bash

conda env create -f workflow/config/main_env.yml --name slidetags
```

Find where it landed:

```bash
conda env list | grep slidetags
# e.g. /n/home12/<user>/.conda/envs/slidetags
```

Set a variable for the rest of setup:
```bash
export SLIDETAGS_ENV=$(conda env list | grep "^slidetags" | awk '{print $2}')
echo $SLIDETAGS_ENV  # verify
```

### 3. Install Julia packages

With the env activated on the interactive node:

```bash
conda activate slidetags

export JULIA_PROJECT_PATH="$HOME/.julia/environments/slidetags"
mkdir -p $JULIA_PROJECT_PATH

julia --project=$JULIA_PROJECT_PATH -e '
using Pkg
Pkg.add(["CSV", "HDF5", "FASTX", "CodecZlib", "IterTools",
         "StatsBase", "DataFrames", "StringViews", "Combinatorics"])
'
```

### 4. Fix HDF5/OpenSSL conflict (required on FASRC)

Julia's HDF5_jll bundles OpenSSL 3.4.0 but the conda env ships OpenSSL 3.3.x.
This causes precompilation to fail. The fix is to replace the artifact's bundled
libssl/libcrypto with symlinks to the conda versions.

First find the artifact directory:

```bash
ARTIFACT_DIR=$(find $SLIDETAGS_ENV/share/julia/artifacts -name "libssl.so" | head -1 | xargs dirname)
echo $ARTIFACT_DIR  # verify it found something
```

Then replace the bundled OpenSSL:

```bash
CONDA_LIB="$SLIDETAGS_ENV/lib"

mv $ARTIFACT_DIR/libssl.so      $ARTIFACT_DIR/libssl.so.bak
mv $ARTIFACT_DIR/libssl.so.3    $ARTIFACT_DIR/libssl.so.3.bak
mv $ARTIFACT_DIR/libcrypto.so   $ARTIFACT_DIR/libcrypto.so.bak
mv $ARTIFACT_DIR/libcrypto.so.3 $ARTIFACT_DIR/libcrypto.so.3.bak

ln -s $CONDA_LIB/libssl.so.3    $ARTIFACT_DIR/libssl.so
ln -s $CONDA_LIB/libssl.so.3    $ARTIFACT_DIR/libssl.so.3
ln -s $CONDA_LIB/libcrypto.so.3 $ARTIFACT_DIR/libcrypto.so
ln -s $CONDA_LIB/libcrypto.so.3 $ARTIFACT_DIR/libcrypto.so.3
```

Then precompile:

```bash
julia --project=$JULIA_PROJECT_PATH -e 'using Pkg; Pkg.precompile()'
```

Should complete with `2 dependencies successfully precompiled`. This is a one-time fix.

> **Note:** A CHOLMOD version mismatch warning may appear during precompile — this is
> harmless and can be ignored.

### 5. Fill in config.sh

Edit `workflow/config/config.sh` with your paths:

```bash
export CLUSTER_PATH=""
export CONDA_PATH="/n/sw/Miniforge3-25.3.1-0/bin"   # system conda on FASRC
export ENV_PATH="$SLIDETAGS_ENV"                      # from step 2
export PKG_PATH="/path/to/software"                   # parent dir containing cellranger-7.x.x/
export BASE_DATA_PATH="/path/to/data"                 # where pipeline outputs will go
export BCL_MAIN_PATH=""
export WORKFLOW_PATH="/path/to/Slide-tags"            # repo root, no trailing slash
export GOOGLE_SHEET_ID=""                             # leave blank if using -csv mode
export GOOGLE_CLOUD_BUCKET=""
export PUCK_PATH="/path/to/pucks"                     # dir containing Puck_*.csv files
export PUCK_IN=""
export REF_PATH="/path/to/reference"                  # dir containing refdata-gex-* folders
```

> **Tip:** `PKG_PATH` should be the directory *containing* `cellranger-7.2.0/`, not the
> cellranger directory itself. The pipeline appends `cellranger-7.2.0/cellranger` to this path.

### 6. Add JULIA_PROJECT_PATH to your .bashrc

The pipeline references this env variable when spawning SLURM jobs:

```bash
echo 'export JULIA_PROJECT_PATH="$HOME/.julia/environments/slidetags"' >> ~/.bashrc
source ~/.bashrc
```

---

## Running the pipeline

### Input CSV format

Create a CSV in `samples/` with one row per sample:

```csv
name,rna_index,sb_index,puck_id,species,chemistry,vdj,cellranger,store_path
mck_9b,SI-TT-F10,SI-NT-F11,Puck_240422_06,Mouse,,,,
```

| Field | Description |
|---|---|
| `name` | Sample name — used for output directories |
| `rna_index` | GEX FASTQ folder name under `data/<BCL>/fastq/` |
| `sb_index` | Spatial barcode FASTQ folder name |
| `puck_id` | Puck CSV filename without `.csv` (must exist in `PUCK_PATH`) |
| `species` | `Mouse` or `Human` (capital first letter) |
| All others | Optional, leave blank |

### FASTQ directory structure

FASTQs must be organized as:
```
$BASE_DATA_PATH/<BCL_NAME>/fastq/<INDEX>/<INDEX>_S*_L*_R[12]_001.fastq.gz
```

Symlinks work fine:
```bash
mkdir -p $BASE_DATA_PATH/mck_9b/fastq/SI-NT-F11
ln -s /path/to/original/spatial/fastqs/* \
    $BASE_DATA_PATH/mck_9b/fastq/SI-NT-F11/
```

### Skipping cellranger (using existing output)

If cellranger has already been run, symlink the output into the expected location:

```bash
mkdir -p $BASE_DATA_PATH/<BCL>/RNAcounts
ln -s /path/to/existing/cellranger/output \
    $BASE_DATA_PATH/<BCL>/RNAcounts/<SAMPLE_NAME>
```

Then start from `-sb` instead of `-cr`.

### Pipeline invocation

```bash
cd /path/to/Slide-tags

# Full run: cellranger → SBcounts → Spatial
bash workflow/slidetag_pipe.sh -bcl <BCL_NAME> -cr -csv samples/my_samples.csv

# Skip cellranger, start from SBcounts
bash workflow/slidetag_pipe.sh -bcl <BCL_NAME> -sb -csv samples/my_samples.csv

# Reuse existing pipeline input files (no CSV needed)
bash workflow/slidetag_pipe.sh -bcl <BCL_NAME> -sb -us

# Spatial positioning only
bash workflow/slidetag_pipe.sh -bcl <BCL_NAME> -sp -us
```

Monitor jobs:
```bash
squeue -u $USER
```

---

## SGE → SLURM changes summary

| File | Key changes |
|---|---|
| `submit_functions.sh` | qsub→sbatch, qstat→squeue, -l h_vmem→--mem, -pe smp→--cpus-per-task, added -p eddy |
| `slidetag_pipe.sh` | conda activate fix, SLURM detection, -csv flag for manual input mode |
| `smk_config.yaml` | Fixed YAML arithmetic strings, fixed runtime units (seconds→minutes) |
| `SBcounts_split.smk` | Removed ulimit, fixed conda PATH bug, fixed range logic |
| `Spatial_split.smk` | Removed ulimit, fixed conda PATH bug |
| `RNAcounts_split.smk` | Removed ulimit, fixed range logic |
| `mkfastq_split.smk` | Removed ulimit, fixed lanlane_nume typo |

---

## Known issues

| Issue | Fix |
|---|---|
| `defaults` channel blocked on FASRC compute nodes | Removed from `main_env.yml` |
| Build string `PackagesNotFoundError` | All build strings stripped from `main_env.yml` |
| Julia HDF5/OpenSSL 3.4 conflict | Artifact symlink fix in step 4 above |
| CHOLMOD version mismatch warning | Harmless, ignore |
| Cellbender requires separate env | Use `cellbender_env.yml`, not covered here |
