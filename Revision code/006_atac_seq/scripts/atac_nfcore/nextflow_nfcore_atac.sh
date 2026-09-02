#!/bin/bash

# Load nextflow and singularity modules
module load HGI/common/nextflow/23.10.0
module load singularityce-3.11.5/python-3.11.6

# Run nf-core atacseq pipeline
nextflow run nf-core/atacseq \
    --input /lustre/scratch126/gengen/projects_v2/escramble/001_007_OTX2/Manuscripts/OTX2/Revision_R_codes/006_atac_seq/scripts/atac_nfcore/nfcore.lookup.csv \
    --outdir /lustre/scratch126/gengen/projects_v2/escramble/001_007_OTX2/Manuscripts/OTX2/Revision_R_codes/006_atac_seq/atac_nf \
    --genome GRCh38 \
    --read_length 150 \
    -profile sanger \
    -c nextflow.config