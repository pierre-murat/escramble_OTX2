#!/bin/bash

# === Configuration ===
POD5_DIR="/lustre/scratch125/gengen/teams_v2/parts/sequencing/260427_pm23_OTX2_adaptive_sampling"
OUTPUT_DIR="/lustre/scratch126/gengen/projects_v2/escramble/001_007_OTX2/Manuscripts/OTX2/Revision_R_codes/004_adaptive_sampling/bam"
DORADO_BIN="/software/hgi/envs/conda/team227/pm23/dorado/bin/dorado"
BATCH_SIZE=50

# === Create output directories ===
mkdir -p "$OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR/logs"

# === Find all pod5 files recursively ===
mapfile -t POD5_FILES < <(find "$POD5_DIR" -name "*.pod5" | sort)

total=${#POD5_FILES[@]}
echo "Found ${total} pod5 file(s) under ${POD5_DIR}"

if [[ $total -eq 0 ]]; then
    echo "ERROR: No pod5 files found. Exiting."
    exit 1
fi

# === Function to submit a batch ===
submit_batch() {
    local -n _batch=$1

    for pod5 in "${_batch[@]}"; do
        fname=$(basename "$pod5")
        sample="${fname%.pod5}"
        bam_out="$OUTPUT_DIR/${sample}.bam"

        bsub -q normal -n 4 -R'select[mem>128000] rusage[mem=128000] span[hosts=1]' -M128000 \
            -J "dorado_${sample}" \
            -o "$OUTPUT_DIR/logs/${sample}.log" \
            -e "$OUTPUT_DIR/logs/${sample}.err" \
            "$DORADO_BIN basecaller fast --device cpu \"$pod5\" > \"$bam_out\""
    done
}

# === Loop through files in batches ===
for (( i=0; i<total; i+=BATCH_SIZE )); do
    batch=("${POD5_FILES[@]:i:BATCH_SIZE}")
    submit_batch batch
    echo "Submitted batch starting at file $((i+1)) of ${total}"
done

echo "All jobs submitted. Monitor with: bjobs -J 'dorado_*'"