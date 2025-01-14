# Rotation Enhancer Scramble Experiments

This repository contains scripts to process long-read nanopore sequencing data from scrambling experiments; and output a VCF file with the resulting structural variants (SVs) from the short-reads sampled from long-reads.



Workflow:

Currently structured based on combining reads across replicates for a given gate into a single FASTQ file; then running the pipeline on each gate, then jointly analyzing all gates.

path_to_long_reads.fastq is a concatenation of all replicates (e.g. for gate very_dim) together.

path/to/reference.fa.gz is the path to the reference genome (GRCh38)

breakpoints.txt is a line seperated file containing the breakpoints of interest)

    python engine/read_processing.py  [path/to/long_reads.fastq] [path/to/reference.fa.gz] [path/to/breakpoints.txt]
    python utilities/read_analysis.py [/path/to/.clean.vcf]

[wrapper.sh was used to run this on the Sanger cluster]


To generate combined .csv files that report all structural variants, run below command. replicates_lookup provides the mapping of read names back to their gates.

    python utilities/gate_sorting.py [/path/to/very_dim.csv] [/path/to/dim.csv]  [/path/to/bright.csv]  [/path/to/very_bright.csv] [/path/to/output_dir/] [/replicates_lookup_table.csv]

Example execution:

    python utilities/gate_sorting.py --very_dim /lustre/scratch126/gengen/teams/parts/am86/data/7breakpoints/experiments/gates_loxpsym_inclusion/very_dim/plots/filtered_data.csv --dim /lustre/scratch126/gengen/teams/parts/am86/data/7breakpoints/experiments/gates_loxpsym_inclusion/dim/plots/filtered_data.csv --bright /lustre/scratch126/gengen/teams/parts/am86/data/7breakpoints/experiments/gates_loxpsym_inclusion/bright/plots/filtered_data.csv --very_bright /lustre/scratch126/gengen/teams/parts/am86/data/7breakpoints/experiments/gates_loxpsym_inclusion/very_bright/plots/filtered_data.csv --output_dir /lustre/scratch126/gengen/teams/parts/am86/data/7breakpoints/experiments/gates_loxpsym_inclusion/plots/  --replicates_lookup /lustre/scratch126/gengen/projects/escramble/Notebook/002_OTX2_superenhancer_scramble/rds/OTX2_loxP7_read_summary.tsv

Final outputs are stored in --output_dir (various .csv files applying progressive filters that have reads grouped into SV architectures)

