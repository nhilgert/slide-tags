#!/bin/bash
# setup_fastq_dirs.sh
# Creates data/<sample>/fastq/<gex_index>/ and data/<sample>/fastq/<sb_index>/
# symlink directories for each sample, pointing to the lanes/ FASTQ symlinks.
#
# Usage: bash setup_fastq_dirs.sh
# Run from: /n/eddy_lab/users/nhilgert/Slide-tags

set -euo pipefail

LANES_PATH="/n/eddy_lab/users/nhilgert/decidualization/lanes"
DATA_PATH="/n/eddy_lab/users/nhilgert/Slide-tags/data"
SAMPLES_CSV="/n/eddy_lab/users/nhilgert/Slide-tags/samples/decidualization_all.csv"

echo "Setting up FASTQ directories for all samples..."
echo ""

tail -n +2 "$SAMPLES_CSV" | while IFS=',' read -r name rna_index sb_index puck_id species rest; do
    lane_dir="$LANES_PATH/$name"

    if [ ! -d "$lane_dir" ]; then
        echo "  SKIP $name — no lane directory at $lane_dir"
        continue
    fi

    gex_src="$lane_dir/gex"
    sb_src="$lane_dir/spatial"

    if [ ! -d "$gex_src" ]; then
        echo "  SKIP $name — no gex dir at $gex_src"
        continue
    fi

    if [ ! -d "$sb_src" ]; then
        echo "  SKIP $name — no spatial dir at $sb_src"
        continue
    fi

    # Create fastq subdirs
    gex_dest="$DATA_PATH/$name/fastq/$rna_index"
    sb_dest="$DATA_PATH/$name/fastq/$sb_index"
    mkdir -p "$gex_dest" "$sb_dest"

    # Symlink all fastq files from lanes/gex -> data/<name>/fastq/<rna_index>/
    for f in "$gex_src"/*.fastq.gz; do
        [ -e "$f" ] || continue
        ln -sf "$f" "$gex_dest/$(basename $f)"
    done

    # Symlink all fastq files from lanes/spatial -> data/<name>/fastq/<sb_index>/
    for f in "$sb_src"/*.fastq.gz; do
        [ -e "$f" ] || continue
        ln -sf "$f" "$sb_dest/$(basename $f)"
    done

    echo "  OK $name  gex=$rna_index ($(ls $gex_dest | wc -l) files)  sb=$sb_index ($(ls $sb_dest | wc -l) files)"
done

echo ""
echo "Done. Verify with: ls $DATA_PATH"
