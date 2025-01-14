
from collections import defaultdict
import pysam


def chop_long_reads_to_short_reads(input_fastq, output_fastq, fragment_length=300, overlap_window=200):
    with open(input_fastq, 'r') as infile, open(output_fastq, 'w') as outfile:
        while True:
            # Read one FASTQ entry (four lines)
            read_id_line = infile.readline().strip()
            if not read_id_line:
                break  # EOF reached
            sequence = infile.readline().strip()
            plus_line = infile.readline().strip()  # separator
            qualities = infile.readline().strip()

            # Split the read ID and take only the main identifier part
            primary_read_id = read_id_line.split()[0]

            # Get the length of the read
            read_length = len(sequence)

            # Chop the read into overlapping fragments
            i = 0
            fragment_index = 0
            while i < read_length:
                # Define the fragment start and end
                chopped_sequence = sequence[i:i + fragment_length]
                chopped_qualities = qualities[i:i + fragment_length]

                # Check if the last fragment is smaller than 1/3 of fragment_length
                if len(chopped_sequence) < fragment_length / 3:
                    break

                new_read_id = f"{primary_read_id}_frag{fragment_index}"

                # Write the short-read fragment to the output FASTQ file
                outfile.write(f"{new_read_id}\n")
                outfile.write(f"{chopped_sequence}\n")
                outfile.write(f"{plus_line}\n")
                outfile.write(f"{chopped_qualities}\n")

                # Move to the next fragment start position with the specified overlap
                i += fragment_length - overlap_window
                fragment_index += 1


# Function to parse a BAM file and group reads by their original long read name
def parse_bam_file(bam_file):
    grouped_reads = defaultdict(list)
    with pysam.AlignmentFile(bam_file, 'rb') as bam:
        for read in bam:
            chrom = bam.get_reference_name(read.reference_id)
            start = read.reference_start
            end = read.reference_end
            read_name = read.query_name
            query_start = read.query_alignment_start
            query_end = read.query_alignment_end
            strand = '-' if read.is_reverse else '+'
            sequence = read.query_sequence
            supplementary = read.is_supplementary
            secondary = read.is_secondary
            read_length = read.query_length
            unmapped = read.is_unmapped
            
            # Extract the original read name by stripping the "_fragX" suffix
            original_read_name = read_name.rsplit('_', 1)[0]
            
            # Store the read information
            grouped_reads[original_read_name].append({
                'chrom': chrom,
                'start': start,
                'end': end,
                'query_start': query_start,
                'query_end': query_end,
                'supplementary': supplementary,
                'secondary': secondary,
                'unmapped': unmapped,
                'read_length': read_length,
                'read_name': read_name,
                'strand': strand,
                'sequence': sequence
            })
            
    
    return grouped_reads
