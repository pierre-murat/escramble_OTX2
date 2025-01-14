import argparse
"""
This script processes long-read FASTQ files by chopping them into shorter reads, aligning them to a reference, and detecting structural variants (SVs).

The script performs the following steps:
1. Parses command-line arguments to get the input FASTQ file and reference file paths.
2. Chops the long reads in the input FASTQ file into shorter reads.
3. Aligns the shorter reads to the reference genome using an external alignment script.
4. Parses the aligned reads from a BED file.
5. Detects structural variants using a list of predefined breakpoints.
6. Generates a grammar for the detected structural variants.
7. Writes the detected structural variants and their grammar to a VCF file.

Arguments:
    input_fastq (str): Path to the long-read FASTQ file.
    reference (str): Path to the reference genome file.

Output:
    A VCF file containing structural variant calls for each read.
"""
import read_manipulation as read_manipulation
import variant_calling as variant_calling
import subprocess

parser = argparse.ArgumentParser(description="Process long-read FASTQ files by chopping them into shorter reads, aligning them to a reference, and detecting structural variants.")
parser.add_argument("input_fastq", type=str, help="Path to long-read fastq file")
parser.add_argument("reference", type=str, help="Path to the reference file")
parser.add_argument("breakpoints", type=str, help="Path to the breakpoints file")

# Parse the arguments
args = parser.parse_args()

input_fastq_base = args.input_fastq.rsplit(".fastq", 1)[0]

cmd = [
    "seqkit", "rmdup", "-n",
    "-o", input_fastq_base + ".removedups.fastq",
    args.input_fastq,
    "-d", input_fastq_base + ".duplicated_list.fastq"
]

# Run the command
subprocess.run(cmd, check=True)


deduped_input_fastq = input_fastq_base + ".removedups.fastq"


read_manipulation.chop_long_reads_to_short_reads(deduped_input_fastq, deduped_input_fastq + ".short.fastq", fragment_length=300, overlap_window=200)

# Call the alignment script

subprocess.run(["./scripts/alignment.sh", args.reference, deduped_input_fastq + ".short.fastq", "4"])


# TODO: make this work for multiple reads at once after ironing out bugs in test cases
with open(args.breakpoints, 'r') as f:
    breakpoints_list = [int(line.strip()) for line in f if line.strip().isdigit()]

print("breakpoints: ", breakpoints_list)
grouped_reads = read_manipulation.parse_bam_file(deduped_input_fastq + ".short.bam")

sv = variant_calling.detect_structural_variants(grouped_reads, breakpoints_list)
grammar_clean, grammar_cpx, loxpsym_dict_clean, loxpsym_dict_complex  = variant_calling.generate_sv_grammar(breakpoints_list, sv)
# Write the output to VCF
coverage = variant_calling.calculate_coverage(grouped_reads, breakpoints_list)
variant_calling.write_vcf(deduped_input_fastq + ".short.grammar.vcf", grammar_clean, grammar_cpx, loxpsym_dict_clean, loxpsym_dict_complex, breakpoints_list, coverage)
