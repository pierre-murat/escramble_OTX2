#!/bin/bash

# ----------------------------------------------------------------------------------------------
# Takes in a fastq file (e.g. short-reads sampled from long-reads via engine/read_processing.py)
# maps to reference, sorts, indexes, and generates a bed file for downstream analysis
# ----------------------------------------------------------------------------------------------

# Check for correct number of arguments
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <reference_path> <input_fastq_path> <threads>"
    exit 1
fi

# Assign input arguments to variables
REFERENCE_PATH=$1
INPUT_FASTQ=$2
THREADS=$3

# Derive the output SAM and BAM paths based on the input FASTQ path
OUTPUT_SAM="${INPUT_FASTQ%.fastq}.sam"
OUTPUT_BAM="${OUTPUT_SAM%.sam}.bam"
OUTPUT_BED="${OUTPUT_SAM%.sam}.bed"
echo ${OUTPUT_SAM}
echo ${OUTPUT_BAM}


# Run BWA to generate SAM file
bwa mem -x ont2d -t $THREADS "$REFERENCE_PATH" "$INPUT_FASTQ" > "$OUTPUT_SAM"

# Sort the SAM file and convert to BAM
samtools sort -@ $THREADS -m 2G "$OUTPUT_SAM" -o "$OUTPUT_BAM"

rm $OUTPUT_SAM

# Index the sorted BAM file
samtools index "$OUTPUT_BAM"
