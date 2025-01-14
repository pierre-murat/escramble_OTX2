#!/bin/bash

input_list=$1
mode=$2

while IFS= read -r file; do
    base_name=$(basename "$file" .txt)
    dir_name=$(dirname "$file")
    echo "$dir_name/$base_name.bam"
    echo "$dir_name/$base_name.bw"
    samtools view -@ 12 -m 2G -bS -N "$file" -o "$dir_name/$base_name.bam" /lustre/scratch126/gengen/projects/escramble/Data/ONT/240731/BAM/OTX2_14_1_scramble_all_reads.bam
    samtools index -@ 12 "$dir_name/$base_name.bam"
    
    if [ "$mode" == "INV" ]; then
        bamCoverage -p 12 -bs 100 --normalizeUsing None --filterRNAstrand reverse --bam "$dir_name/$base_name.bam" -o  "$dir_name/$base_name.fwd.bw"
        bamCoverage -p 12 -bs 100 --normalizeUsing None --filterRNAstrand forward --bam "$dir_name/$base_name.bam" -o "$dir_name/$base_name.rev.bw"
    else
        bamCoverage -p 12 -bs 100 --normalizeUsing None --bam "$dir_name/$base_name.bam" -o "$dir_name/$base_name.bw"
    fi

done < "$input_list"
