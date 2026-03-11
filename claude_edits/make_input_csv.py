#!/usr/bin/env python3
"""
make_input_csv.py — Manual CSV input generator for Slide-tags pipeline (FASRC/SLURM fork)

Replaces the Google Sheets dependency in read_info.py.
Takes a simple user-provided CSV and generates all input files the pipeline expects.

Usage:
    python make_input_csv.py <input_csv> <bcl> <data_path> <log_folder>

Input CSV format (with header):
    name,rna_index,sb_index,puck_id,species,chemistry,vdj,cellranger,store_path

    name        : Sample name (must be unique)
    rna_index   : RNA/GEX FASTQ index folder name (e.g. SI-TT-A1)
    sb_index    : Spatial barcode FASTQ index folder name; leave blank or X if no spatial
    puck_id     : Puck CSV filename without extension (must match file in PUCK_PATH); X if no spatial
    species     : Mouse or Human
    chemistry   : Cellranger chemistry string (default: auto)
    vdj         : VDJ chain (e.g. TR, IG); leave blank if none
    cellranger  : V7 or V8 (default: V7)
    store_path  : Path to move final results; leave blank to skip move_data step

Outputs (written to <log_folder>/input/):
    name_to_index.csv
    RNAcounts.tsv
    SBcounts.tsv
    Spatial.tsv
    split_counts/count_<rna_index>.csv    (one per RNA sample)
    split_spatial/SBcount_<sb_index>.csv  (one per SB sample)
"""

import os
import sys
import csv

# ---------------------------------------------------------------------------
# Reference genome name mapping
# ---------------------------------------------------------------------------
TRANSCRIPTOME = {
    ("Mouse", False, False): "refdata-gex-mm10-2020-A",
    ("Mouse", True,  False): "refdata-cellranger-vdj-GRCm38",
    ("Mouse", False, True):  "refdata-arc-mm10-2020-A",
    ("Human", False, False): "refdata-gex-GRCh38-2024-A",
    ("Human", True,  False): "refdata-cellranger-vdj-GRCh38",
    ("Human", False, True):  "refdata-arc-GRCh38-2020-A",
}

EMPTY = {"X", "x", "", None}


def resolve_transcriptome(species, vdj, atac=False):
    key = (species, bool(vdj and vdj not in EMPTY), atac)
    ref = TRANSCRIPTOME.get(key)
    if ref is None:
        raise ValueError(f"Unknown species '{species}' — must be Mouse or Human")
    return ref


def ensure_dir(path):
    os.makedirs(path, exist_ok=True)


def read_input_csv(input_csv):
    required = {"name", "rna_index", "sb_index", "puck_id", "species"}
    rows = []
    with open(input_csv, newline="") as f:
        reader = csv.DictReader(f)
        headers = set(reader.fieldnames or [])
        missing = required - headers
        if missing:
            raise ValueError(f"Input CSV missing required columns: {missing}\n"
                             f"  Got: {sorted(headers)}\n"
                             f"  Required: {sorted(required)}")
        for i, row in enumerate(reader, start=2):
            # strip whitespace from all values
            row = {k: (v.strip() if v else "") for k, v in row.items()}
            if not row["name"]:
                raise ValueError(f"Row {i}: 'name' column is empty")
            if not row["rna_index"] or row["rna_index"] in EMPTY:
                raise ValueError(f"Row {i} ('{row['name']}'): 'rna_index' is required")
            if not row["species"] or row["species"] not in ("Mouse", "Human"):
                raise ValueError(f"Row {i} ('{row['name']}'): 'species' must be Mouse or Human, got '{row['species']}'")
            rows.append(row)

    # Check for duplicate names
    names = [r["name"] for r in rows]
    seen = set()
    for n in names:
        if n in seen:
            raise ValueError(f"Duplicate sample name: '{n}'")
        seen.add(n)

    return rows


def generate_files(rows, log_folder, data_path, bcl):
    input_dir = os.path.join(log_folder, "input")
    split_counts_dir = os.path.join(input_dir, "split_counts")
    split_spatial_dir = os.path.join(input_dir, "split_spatial")

    for d in [input_dir, split_counts_dir, split_spatial_dir]:
        ensure_dir(d)

    # Paths for master output files
    name_to_index_path = os.path.join(input_dir, "name_to_index.csv")
    rna_counts_path    = os.path.join(input_dir, "RNAcounts.tsv")
    sb_counts_path     = os.path.join(input_dir, "SBcounts.tsv")
    spatial_path       = os.path.join(input_dir, "Spatial.tsv")

    name_to_index_rows = []
    rna_counts_rows    = []
    sb_counts_rows     = []
    spatial_rows       = []

    for row in rows:
        name        = row["name"]
        rna_index   = row["rna_index"]
        sb_index    = row.get("sb_index", "X") or "X"
        puck_id     = row.get("puck_id", "X") or "X"
        species     = row["species"]
        chemistry   = row.get("chemistry", "auto") or "auto"
        vdj         = row.get("vdj", "") or ""
        cellranger  = row.get("cellranger", "V7") or "V7"
        store_path  = row.get("store_path", "") or ""

        has_spatial = sb_index not in EMPTY and puck_id not in EMPTY

        # Validate spatial consistency
        sb_empty   = sb_index in EMPTY
        puck_empty = puck_id in EMPTY
        if sb_empty != puck_empty:
            raise ValueError(
                f"Sample '{name}': sb_index and puck_id must both be filled or both be empty/X. "
                f"Got sb_index='{sb_index}', puck_id='{puck_id}'"
            )

        # Resolve transcriptome
        transcriptome = resolve_transcriptome(species, vdj)

        # Adjust chemistry string (cellranger expects SC prefix for non-auto/ARC)
        if chemistry and chemistry not in ("auto", "ARC-v1"):
            chem_str = f"SC{chemistry}"
        else:
            chem_str = chemistry

        # Resolve VDJ
        vdj_str = vdj if vdj and vdj not in EMPTY else ""

        # name_to_index.csv row
        name_to_index_rows.append({
            "Name": name,
            "RNA Index": rna_index,
            "VDJ Chain": vdj_str or "X",
            "SB Index": sb_index if has_spatial else "X",
            "Puck ID": puck_id if has_spatial else "X",
            "Store Path": store_path or "X",
        })

        # RNAcounts.tsv (tab-separated, no header)
        rna_counts_rows.append(f"{rna_index}\t{transcriptome}")

        # split_counts/count_<rna_index>.csv
        split_count_path = os.path.join(split_counts_dir, f"count_{rna_index}.csv")
        with open(split_count_path, "w") as f:
            f.write(f"{rna_index},{transcriptome},{chem_str},{vdj_str},{cellranger}\n")

        # Spatial files
        if has_spatial:
            sb_counts_rows.append(sb_index)
            spatial_rows.append(f"{sb_index}\t{puck_id}")
            # split_spatial/SBcount_<sb_index>.csv
            split_sb_path = os.path.join(split_spatial_dir, f"SBcount_{sb_index}.csv")
            with open(split_sb_path, "w") as f:
                f.write(f"{sb_index}\n")

    # Write name_to_index.csv
    with open(name_to_index_path, "w", newline="") as f:
        fieldnames = ["Name", "RNA Index", "VDJ Chain", "SB Index", "Puck ID", "Store Path"]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(name_to_index_rows)

    # Write RNAcounts.tsv
    with open(rna_counts_path, "w") as f:
        f.write("\n".join(rna_counts_rows) + "\n")

    # Write SBcounts.tsv (blank file if no spatial samples)
    with open(sb_counts_path, "w") as f:
        if sb_counts_rows:
            f.write("\n".join(sb_counts_rows) + "\n")

    # Write Spatial.tsv (blank file if no spatial samples)
    with open(spatial_path, "w") as f:
        if spatial_rows:
            f.write("\n".join(spatial_rows) + "\n")

    # Summary
    print(f"Input sheets saved in: {input_dir}")
    print(f"RNA samples   : {len(rna_counts_rows)}")
    print(f"Spatial samples: {len(sb_counts_rows)}")
    print(f"Files written :")
    print(f"  {name_to_index_path}")
    print(f"  {rna_counts_path}")
    print(f"  {sb_counts_path}")
    print(f"  {spatial_path}")
    print(f"  {split_counts_dir}/count_*.csv  ({len(rna_counts_rows)} files)")
    if sb_counts_rows:
        print(f"  {split_spatial_dir}/SBcount_*.csv  ({len(sb_counts_rows)} files)")


def main():
    if len(sys.argv) != 5:
        print(__doc__)
        sys.exit(1)

    input_csv, bcl, data_path, log_folder = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]

    if not os.path.exists(input_csv):
        print(f"ERROR: Input CSV not found: {input_csv}")
        sys.exit(1)

    try:
        rows = read_input_csv(input_csv)
        generate_files(rows, log_folder, data_path, bcl)
    except ValueError as e:
        print(f"ERROR: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
