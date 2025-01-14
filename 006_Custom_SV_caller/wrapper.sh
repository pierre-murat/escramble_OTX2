#!/bin/bash

n=4
M=70000
J="scrambleSV"
q="long"

logdir="/lustre/scratch126/gengen/teams/parts/am86/data/7breakpoints/experiments/gates_loxpsym_inclusion/logs/"
mkdir -p $logdir

ref="/lustre/scratch126/gengen/teams/parts/pm23/Genomes/hg38/hg38.fa.gz"
breakpoints="/lustre/scratch126/gengen/teams/parts/am86/data/7breakpoints/breakpoints/7breakpoints.txt"

module load HGI/softpack/groups/escramble/eSCRAMBLE/7
source env/bin/activate

while IFS=',' read -r lr_fastq; do
        #echo "submitting job for $lr_fastq"
        lr_fastq_base=$(basename "$lr_fastq" .fastq)
        lr_fastq_end=$(basename "$(dirname "$lr_fastq")")
        output_dir=$(dirname $lr_fastq)
        cmd="python engine/read_processing.py $lr_fastq $ref $breakpoints && python utilities/read_analysis.py ${output_dir}/${lr_fastq_base}.removedups.fastq.short.grammar.vcf.clean.vcf"
        echo $cmd
        bsub -q $q -n$n -M$M -R"select[mem>$M] rusage[mem=$M]" -J $lr_fastq_end \
                -o $logdir/$lr_fastq_end.log.out -e $logdir/$lr_fastq_end.log.err \
                eval $cmd
done < "$1"
