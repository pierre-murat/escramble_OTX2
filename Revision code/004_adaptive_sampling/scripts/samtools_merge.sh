#!/bin/bash
#BSUB -J OTX2_merge_bam_job
#BSUB -q normal
#BSUB -n 12
#BSUB -R 'select[mem>32000] rusage[mem=32000] span[hosts=1]'
#BSUB -o ./logs/OTX2_merge_bam.%J.out
#BSUB -e ./logs/OTX2_merge_bam.%J.err

# === Define working directory ===
cd /lustre/scratch126/gengen/projects_v2/escramble/001_007_OTX2/Manuscripts/OTX2/Revision_R_codes/004_adaptive_sampling

# === Configuration ===
BAM_DIR="./bam_aligned"
FINAL_DIR="./bam_merged"
COVERAGE_DIR="./bw"
FINAL_BAM="$FINAL_DIR/OTX2_as.sorted.bam"
COVERAGE_OUT="$COVERAGE_DIR/OTX2_as.1kb.bw"
THREADS=12

# === Setup ===
mkdir -p "$FINAL_DIR"
mkdir -p "$FINAL_DIR/logs"
mkdir -p "$COVERAGE_DIR"

# === Load samtools and deeptools ===
module load HGI/softpack/groups/escramble/eSCRAMBLE/9

# === Generate BAM file list, excluding any pre-existing merged output ===
find "$BAM_DIR" -maxdepth 1 -name "*.bam" ! -name "*.sorted.bam" | sort > "$FINAL_DIR/logs/bam_list.txt"

total=$(wc -l < "$FINAL_DIR/logs/bam_list.txt")
echo ">> Found ${total} BAM file(s) to merge"

if [[ $total -eq 0 ]]; then
    echo "ERROR: No BAM files found. Exiting."
    exit 1
fi

# === Merge ===
echo ">> Merging all BAM files..."
samtools merge -@ "$THREADS" -f -b "$FINAL_DIR/logs/bam_list.txt" "$FINAL_DIR/merged.bam"

# === Sort ===
echo ">> Sorting merged BAM..."
samtools sort -@ "$THREADS" -o "$FINAL_BAM" "$FINAL_DIR/merged.bam"

# === Index ===
echo ">> Indexing final BAM..."
samtools index "$FINAL_BAM"

# === Remove intermediate merged BAM ===
echo ">> Removing intermediate merged BAM..."
rm -f "$FINAL_DIR/merged.bam"

# === Compute coverage track ===
echo ">> Computing coverage track at 1000 bp resolution..."
bamCoverage \
    --bam "$FINAL_BAM" \
    --outFileName "$COVERAGE_OUT" \
    --outFileFormat bigwig \
    --binSize 1000 \
    --numberOfProcessors "$THREADS" \
    --normalizeUsing None \
    --skipNonCoveredRegions

echo ">> Finished. Output BAM: $FINAL_BAM"
echo ">> Finished. Output coverage: $COVERAGE_OUT"