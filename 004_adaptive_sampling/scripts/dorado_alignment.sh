#!/bin/bash

# === Configuration ===
DORADO_BIN="/software/hgi/envs/conda/team227/pm23/dorado/bin/dorado"
REFERENCE="/lustre/scratch125/gengen/teams_v2/parts/pm23/genomes/hg38/hg38.fa.gz"
INPUT_DIR="/lustre/scratch126/gengen/projects_v2/escramble/001_007_OTX2/Manuscripts/OTX2/Revision_R_codes/004_adaptive_sampling/bam"
OUTPUT_DIR="/lustre/scratch126/gengen/projects_v2/escramble/001_007_OTX2/Manuscripts/OTX2/Revision_R_codes/004_adaptive_sampling/bam_aligned"
BATCH_SIZE=50

# === Create output directories ===
mkdir -p "$OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR/logs"

# === Load samtools ===
module load HGI/softpack/groups/escramble/eSCRAMBLE/9

# === Find all unaligned BAM files ===
mapfile -t BAM_FILES < <(find "$INPUT_DIR" -maxdepth 1 -name "*.bam" | sort)

total=${#BAM_FILES[@]}
echo "Found ${total} BAM file(s) to align"

if [[ $total -eq 0 ]]; then
    echo "ERROR: No BAM files found. Exiting."
    exit 1
fi

# === Function to submit a batch ===
submit_batch() {
    local -n _batch=$1

    for bam in "${_batch[@]}"; do
        sample=$(basename "$bam" .bam)
        bam_out="$OUTPUT_DIR/${sample}.bam"

        bsub -q normal -n 4 \
            -R'select[mem>16000] rusage[mem=16000] span[hosts=1]' \
            -M 16000 \
            -J "dorado_align_${sample}" \
            -o "$OUTPUT_DIR/logs/${sample}.log" \
            -e "$OUTPUT_DIR/logs/${sample}.err" \
            "$DORADO_BIN aligner --threads 4 \"$REFERENCE\" \"$bam\" \
            | samtools sort -@ 2 -o \"$bam_out\" && \
            samtools index \"$bam_out\""
    done
}

# === Loop through BAM files in batches ===
for (( i=0; i<total; i+=BATCH_SIZE )); do
    batch=("${BAM_FILES[@]:i:BATCH_SIZE}")
    submit_batch batch
    echo "Submitted batch starting at file $((i+1)) of ${total}"
done

echo "All alignment jobs submitted. Monitor with: bjobs -J 'dorado_align_*'"