import pysam
import sys
import os

def filter_reads_by_name(input_bam, read_names_list, output_dir):
    # Open input BAM file with an index for faster random access
    with pysam.AlignmentFile(input_bam, "rb") as in_bam:
        # Create a dictionary to store file handles for each label
        print(f"Header for output BAM: {in_bam.header}")

        output_files = {label: pysam.AlignmentFile(os.path.join(output_dir, f"{label}_{read_name}.bam"), "wb", header=in_bam.header) for read_name, label in read_names_list}
        
        
        # Iterate over reads only if the BAM file is indexed
        for read in in_bam.fetch(until_eof=True):  # until_eof improves speed for unsorted, unindexed BAMs
            for architecture, label in read_names_list:
                #print("architecture, label :", architecture, label)
                if label in read.query_name:
                    print("found occurent of label in read.query_name", label, read.query_name)
                    print(f"Read info: {read.query_name}, Flag: {read.flag}, Aligned Pairs: {len(read.get_aligned_pairs())}")  # Debug read details
                    print(f"Writing to {output_files[label]}")  # Debug print for write action
                    output_files[label].write(read)
    
    # Close all output BAM files
    for out_bam in output_files.values():
        out_bam.close()
    print("exiting")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python visualize.py <input_bam> <read_names_file> <output_dir>")
        sys.exit(1)

    input_bam = sys.argv[1]
    read_names_file = sys.argv[2]
    output_dir = sys.argv[3]
    
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    with open(read_names_file, 'r') as f:
        read_names_list = [line.strip().split('\t') for line in f]
    print("read names list: ", read_names_list)
    
    # Run filtering function
    filter_reads_by_name(input_bam, read_names_list, output_dir)

