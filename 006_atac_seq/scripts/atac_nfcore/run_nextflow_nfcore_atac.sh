#!/bin/bash
M=2000
n=1
q="long"
J="nextflow_nfcore_atac"
mkdir -p logs
bsub -q $q \
        -n$n \
        -M$M -R"select[mem>$M] rusage[mem=$M]" \
        -J $J \
        -o logs/$J.log.out -e logs/$J.log.err \
        bash nextflow_nfcore_atac.sh
echo "Submitted job $J"

