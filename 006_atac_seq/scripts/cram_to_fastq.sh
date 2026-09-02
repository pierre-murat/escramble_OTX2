#!/bin/bash

# === Configuration ===
SAMPLESHEET="/lustre/scratch126/gengen/projects_v2/escramble/001_007_OTX2/Manuscripts/OTX2/Revision_R_codes/006_atac_seq/scripts/cram_to_fastq_table.csv"
OUTDIR="/lustre/scratch126/gengen/projects_v2/escramble/001_007_OTX2/Manuscripts/OTX2/Revision_R_codes/006_atac_seq/fastq"
n=12
M=16000
q="normal"

# === Create output directories ===
mkdir -p "$OUTDIR"
mkdir -p "$OUTDIR/logs"

while IFS=',' read -r SAMPLE REP PLEX PATH_CRAM; do
    J="cram_${SAMPLE}_${REP}"
    cmd="
    module load HGI/softpack/groups/escramble/eSCRAMBLE/9

    samtools sort -n -@ ${n} ${PATH_CRAM} | \
    samtools fastq - \
        -@ ${n} \
        -1 ${OUTDIR}/${SAMPLE}_rep_${REP}_R1.fastq.gz \
        -2 ${OUTDIR}/${SAMPLE}_rep_${REP}_R2.fastq.gz \
        -0 ${OUTDIR}/${SAMPLE}_rep_${REP}_unpaired.fastq.gz \
        -s ${OUTDIR}/${SAMPLE}_rep_${REP}_singleton.fastq.gz \
        -N -F 0x900
    "

    bsub -q "$q" \
         -n "$n" \
         -M "$M" \
         -R "select[mem>$M] rusage[mem=$M] span[hosts=1]" \
         -J "$J" \
         -o "$OUTDIR/logs/$J.log.out" \
         -e "$OUTDIR/logs/$J.log.err" \
         bash -c "${cmd}"

done < "$SAMPLESHEET"

echo "All CRAM conversion jobs submitted. Monitor with: bjobs -J 'cram_*'"