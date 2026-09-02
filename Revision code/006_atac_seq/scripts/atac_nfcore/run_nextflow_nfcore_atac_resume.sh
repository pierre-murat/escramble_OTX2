#!/bin/bash
M=2000
n=1
q="normal"
J="nextflow_nfcore_atac_resume"
mkdir -p logs
bsub -q $q \
        -n$n \
        -M$M -R"select[mem>$M] rusage[mem=$M]" \
        -J $J \
        -o logs/$J.log.out -e logs/$J.log.err \
        bash nextflow_nfcore_atac_resume.sh
echo "Submitted job $J"


