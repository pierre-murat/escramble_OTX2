def sort_reads(reads):
    # Define the sorting key
    def sorting_key(read):
        # Extract frag number from the read name (assuming the format '..._fragX')
        frag_number = int(read['read_name'].rsplit('_frag', 1)[1])
        start_position = read['start']
        # Return tuple for sorting
        return (frag_number, start_position)

    # Sort using the sorting key
    return sorted(reads, key=sorting_key)


# Function to find the closest breakpoint to a given position
# Modified function to find the closest single breakpoint to a given position
def find_nearest_breakpoint(read, breakpoints, event_distance_threshold=1000):
    """Return the closest breakpoint to a given read's start and end position.
    If the closest breakpoint is more than 1000 distance from both positions, return None.
    Also return whether most of the read is to the left or right of the found breakpoint."""
    start = read['start']
    end = read['end']
    
    closest_bp_start = min(breakpoints, key=lambda bp: abs(bp - start))
    closest_bp_end = min(breakpoints, key=lambda bp: abs(bp - end))
    
    # Choose the closest breakpoint out of the two positions
    if abs(closest_bp_start - start) <= abs(closest_bp_end - end):
        closest_bp = closest_bp_start
    else:
        closest_bp = closest_bp_end
    
    if abs(closest_bp - start) > event_distance_threshold and abs(closest_bp - end) > event_distance_threshold:
        return None, None
    
    
    # Determine whether most of the read is on the left or right side of the found breakpoint
    read_midpoint = (start + end) / 2
    direction = 'right' if read_midpoint >= closest_bp else 'left'
    return closest_bp, direction



# Function to detect structural variants between known breakpoints
from collections import defaultdict

def calculate_coverage(grouped_reads, breakpoints):
    breakpoints.insert(0, breakpoints[0] - 40)  # Add item to the start
    breakpoints.append(breakpoints[-1] + 40)  # Add item to the end
    coverage_results = defaultdict(lambda: {'forward': {}, 'reverse': {}})


    for original_read_name, reads in grouped_reads.items():
        # Sort the reads by their start position
        reads = sort_reads(reads)

        # Initialize coverage calculation
        for i in range(len(breakpoints) - 1):
            start_bp = breakpoints[i]
            end_bp = breakpoints[i + 1]
            span = end_bp - start_bp
            total_coverage_forward = 0
            total_coverage_reverse = 0

            for read in reads:
                
                if read["unmapped"] or read["secondary"]:
                    continue
         
                if read['start'] <= end_bp and read['end'] >= start_bp:
                    overlap_start = max(read['start'], start_bp)
                    overlap_end = min(read['end'], end_bp)
                    coverage = overlap_end - overlap_start + 1
                    if read['strand'] == '+':
                        total_coverage_forward += coverage
                    else:
                        total_coverage_reverse += coverage

            average_coverage_forward = total_coverage_forward / span if span > 0 else 0
            average_coverage_reverse = total_coverage_reverse / span if span > 0 else 0
            coverage_results[original_read_name]['forward'][(start_bp, end_bp)] = average_coverage_forward
            coverage_results[original_read_name]['reverse'][(start_bp, end_bp)] = average_coverage_reverse

    return coverage_results


from Levenshtein import distance as levenshtein_distance
import sys

def detect_loxpysm_sequence(read_a, read_b):

    # Define the loxpsym components
    loxpsym = "ATAACTTCGTATAATGTACATTATACGAAGTTAT"
    left_flank = "ATAACTTCGTATA"
    core = "ATGTACAT"
    right_flank = "TATACGAAGTTAT"

    # Function to check approximate match for any part of loxpsym site
    def find_approximate_match(target, query, max_mismatches=2):
        target_length = len(target)
        for i in range(len(query) - target_length + 1):
            substring = query[i:i + target_length]
            mismatches = levenshtein_distance(target, substring)
            if mismatches <= max_mismatches:
                return True, substring, mismatches
        return False, None, None

    # Function to check if each loxpsym component is contained within the input sequence
    def check_loxpsym_sites(input_sequence):
        # Check for each part of the loxpsym sequence (left_flank, core, right_flank)
        results = {}
        for site_name, site_sequence in [("loxpsym", loxpsym), ("left_flank", left_flank), ("core", core), ("right_flank", right_flank)]:
            max_mismatches = 2 if site_name == "loxpsym" else 1
            found, _, _ = find_approximate_match(site_sequence, input_sequence, max_mismatches)
            results[site_name] = found
        return results
    

    # Run the function to check for loxpsym sites in the combined sequence
    results_a = check_loxpsym_sites(read_a["sequence"])
    results_b = check_loxpsym_sites(read_b["sequence"])
    combined_results = {}
    for site in results_a.keys():
        combined_results[site] = results_a[site] or results_b[site]
    return combined_results


def detect_structural_variants(grouped_reads, breakpoints):
    """
    # Deletion signature; no overlap, same strand, different direction (earlier bp left, later bp right), copy number change (decrease)
    # Inv signature; no overlap, different strand, same direction, no copy number change
    # Simple duplication; no overlap, same strand, different direction (earlier bp right, later bp left), copy number change (increase)
    # Inverted duplication; overlap, different strand, same direction, copy number change (increase)
    """ 

    sv_results = defaultdict(list)

    total_reads = len(grouped_reads)
    processed_reads = 0

    for original_read_name, reads in grouped_reads.items():
        processed_reads += 1
        if processed_reads % 10000 == 0 or processed_reads == total_reads:
            print(f"Processed {processed_reads}/{total_reads} reads", file=sys.stderr)
        # Sort the reads by their start position
        reads = sort_reads(reads)

        sv_results[original_read_name] = []
        
        # Initialize structural variant detection
        for i in range(len(reads) - 1):
            current_read = reads[i]
            next_read = reads[i + 1]

            # if any short-read doesn't align, less confident and should skip such cases
            if current_read['secondary'] or next_read['secondary']:

                print("ERROR: secondary alignment detected")
                sv_results[original_read_name].append({
                    'type': 'secondary_alignment',
                    'chrom': current_read['chrom'],
                    'start': current_read['start'],
                    'end': next_read['start'],
                    'read_name': original_read_name,
                    'loxpsym': None
                })
                continue   

            # if any short-read doesn't align, less confident and should skip such cases
            if current_read['unmapped'] or next_read['unmapped']:

                print("ERROR: unmapped read detected")
                sv_results[original_read_name].append({
                    'type': 'unmapped_read',
                    'chrom': current_read['chrom'],
                    'start': current_read['start'],
                    'end': next_read['start'],
                    'read_name': original_read_name,
                    'loxpsym': None
                })
                continue  


            if current_read['chrom'] != next_read['chrom']:

                print("ERROR: cross chromosomal event detected", current_read['chrom'], next_read['chrom'])
                sv_results[original_read_name].append({
                    'type': 'cross_chrom',
                    'chrom': current_read['chrom'],
                    'start': current_read['start'],
                    'end': next_read['start'],
                    'read_name': original_read_name,
                    'loxpsym': None
                })
                continue
        
            if current_read['chrom'] != 'chr14' or next_read['chrom'] != 'chr14':
                print("ERROR: outside of chr14", current_read['chrom'], next_read['chrom'])
                sv_results[original_read_name].append({
                    'type': 'event_outside_chr14',
                    'chrom': current_read['chrom'],
                    'start': current_read['start'],
                    'end': next_read['start'],
                    'read_name': original_read_name,
                    'loxpsym': None
                })
                continue    
            
            
            # check for overlap
            overlap = False
            strand_flip = current_read['strand'] != next_read['strand']
             
            overlap = (current_read['start'] <= next_read['end']) and (next_read['start'] <= current_read['end'])

            # TODO: iron this out
            if current_read['read_name'] != next_read['read_name']:
                continue

            # If reads overlap and there's no strand flip, continue to next read
            if overlap and not strand_flip:
                continue


            loxpsym_results = detect_loxpysm_sequence(current_read, next_read)

            
            """
                # Deletion signature; no overlap, same strand, different direction (earlier bp left, later bp right), copy number change (decrease)
                # Simple duplication; no overlap, same strand, different direction (earlier bp right, later bp left), copy number change (increase)
            """ 
            # Case 1: Deletion + Duplication (no overlap, same strand)
            if not overlap and not strand_flip:
                print("DETECTED DELETION")
                # Find breakpoints near current_read end and next_read start
                del_start_bp, start_direction = find_nearest_breakpoint(current_read, breakpoints)
                del_end_bp, end_direction = find_nearest_breakpoint(next_read, breakpoints)
                
                if not (del_start_bp and del_end_bp):
                    print("ERROR: DEL/DUP far from breakpoint")
                    sv_results[original_read_name].append({
                    'type': 'DEL_DUP_far_from_breakpoint',
                    'chrom': current_read['chrom'],
                    'start': current_read['start'],
                    'end': next_read['start'],
                    'read_name': original_read_name,
                    'loxpsym': loxpsym_results
                    })
                    continue

                earlier = min(del_start_bp, del_end_bp)

                # Deletion
                if del_start_bp == earlier:
                    if start_direction == 'left' and end_direction == 'right': 
                        sv_results[original_read_name].append({
                            'type': 'DEL',
                            'chrom': current_read['chrom'],
                            'start': del_start_bp,
                            'end': del_end_bp,
                            'read_name': original_read_name,
                            'loxpsym': loxpsym_results
                        })
                    elif start_direction == 'right' and end_direction == 'left':
                        sv_results[original_read_name].append({
                            'type': 'DUP',
                            'chrom': current_read['chrom'],
                            'start': del_start_bp,
                            'end': del_end_bp,
                            'read_name': original_read_name,
                            'loxpsym': loxpsym_results
                        })
                    else:
                        sv_results[original_read_name].append({
                            'type': 'DEL_DUP_breakpoints_direction_mismatch',
                            'chrom': current_read['chrom'],
                            'start': del_start_bp,
                            'end': del_end_bp,
                            'read_name': original_read_name,
                            'loxpsym': loxpsym_results
                        })
                else:
                    if start_direction == 'left' and end_direction == 'right': 
                        sv_results[original_read_name].append({
                            'type': 'DUP',
                            'chrom': current_read['chrom'],
                            'start': del_start_bp,
                            'end': del_end_bp,
                            'read_name': original_read_name,
                            'loxpsym': loxpsym_results
                        })
                    elif start_direction == 'right' and end_direction == 'left':
                        sv_results[original_read_name].append({
                            'type': 'DEL',
                            'chrom': current_read['chrom'],
                            'start': del_start_bp,
                            'end': del_end_bp,
                            'read_name': original_read_name,
                            'loxpsym': loxpsym_results
                        })
                    else:
                        sv_results[original_read_name].append({
                            'type': 'DEL_DUP_breakpoints_direction_mismatch',
                            'chrom': current_read['chrom'],
                            'start': del_start_bp,
                            'end': del_end_bp,
                            'read_name': original_read_name,
                            'loxpsym': loxpsym_results
                        })

                """
                # Inv signature; no overlap, different strand, same direction, no copy number change
                # Inverted duplication; overlap, different strand, same direction, copy number change (increase)
                """ 
            # Case: Inversion (no overlap, with strand flip)
            elif not overlap and strand_flip:
                print("DETECTED INVERSION")
                
                # Find breakpoints near current_read end and next_read start

                inv_start_bp, start_direction = find_nearest_breakpoint(current_read, breakpoints)
                inv_end_bp, end_direction = find_nearest_breakpoint(next_read, breakpoints)
                
                if not (inv_start_bp and inv_end_bp):
                    print("ERROR: INV far from breakpoint")
                    sv_results[original_read_name].append({
                    'type': 'INV_far_from_breakpoint',
                    'chrom': current_read['chrom'],
                    'start': current_read['start'],
                    'end': next_read['start'],
                    'read_name': original_read_name,
                    'loxpsym': loxpsym_results
                    })
                    continue

                if start_direction == end_direction:
                    sv_results[original_read_name].append({
                        'type': 'INV',
                        'chrom': current_read['chrom'],
                        'start': inv_start_bp,
                        'end': inv_end_bp,
                        'read_name': original_read_name,
                        'loxpsym': loxpsym_results
                    })
                else:
                    sv_results[original_read_name].append({
                        'type': 'INV_breakpoints_direction_mismatch',
                        'chrom': current_read['chrom'],
                        'start': inv_start_bp,
                        'end': inv_end_bp,
                        'read_name': original_read_name,
                        'loxpsym': loxpsym_results
                    })

            # Case: Overlap and strand flip (potential complex event)
            elif overlap and strand_flip:
                # Inverted duplication; overlap, different strand, same direction, copy number change (increase)
                # Find breakpoints near current_read end and next_read start
                inv_start_bp, start_direction = find_nearest_breakpoint(current_read, breakpoints)
                inv_end_bp, end_direction = find_nearest_breakpoint(next_read, breakpoints)
                
                if not (inv_start_bp and inv_end_bp):
                    print("ERROR: INVDUP far from breakpoint")
                    sv_results[original_read_name].append({
                    'type': 'INVDUP_far_from_breakpoint',
                    'chrom': current_read['chrom'],
                    'start': current_read['start'],
                    'end': next_read['start'],
                    'read_name': original_read_name,
                    'loxpsym': loxpsym_results
                    })
                    continue
                if start_direction == end_direction:
                    sv_results[original_read_name].append({
                        'type': 'INVDUP',
                        'chrom': current_read['chrom'],
                        'start': inv_start_bp,
                        'end': inv_end_bp,
                        'read_name': original_read_name,
                        'loxpsym': loxpsym_results
                    })
                else:
                    sv_results[original_read_name].append({
                        'type': 'INVDUP_breakpoints_direction_mismatch',
                        'chrom': current_read['chrom'],
                        'start': inv_start_bp,
                        'end': inv_end_bp,
                        'read_name': original_read_name,
                        'loxpsym': loxpsym_results
                    })
            else:
                print("COMPLEX SIGNAL")
                # Handle or log unhandled cases as necessary 
                sv_results[original_read_name].append({
                    'type': 'COMPLEX_SV',
                    'chrom': current_read['chrom'],
                    'start': current_read['start'],
                    'end': next_read['start'],
                    'read_name': original_read_name,
                    'loxpsym': loxpsym_results
                })     
    return sv_results

    
def generate_sv_grammar(breakpoints, sv_results):
    # Start with sequential labels A, B, C, etc., for breakpoints
    grammar_labels = ['A', 'B', 'C', 'D', 'E', 'F', 'G']
    
    # Map breakpoints to grammar labels
    breakpoint_labels = {bp: grammar_labels[i] for i, bp in enumerate(breakpoints)}

    # Dictionary to store grammar for each original read name
    grammar_dict_clean = {}
    grammar_dict_cpx = {}
    loxpsym_dict_clean = {}
    loxpsym_dict_complex = {}


    for original_read_name, sv_calls in sv_results.items():
        cpx = False
        grammar = ""
        loxpsym = ""
        # going through all SV calls for a given read; and stitching together a grammar
        for sv in sv_calls:
            # Add grammar for inversion
            if sv['type'] == 'INV':
                grammar += f"{breakpoint_labels[sv['start']]}^{breakpoint_labels[sv['end']]},"
            # Add grammar for deletion
            elif sv['type'] == 'DEL':
                grammar += f"{breakpoint_labels[sv['start']]}_{breakpoint_labels[sv['end']]},"
            # Add grammar for duplication
            elif sv['type'] == 'DUP':
                grammar += f"{breakpoint_labels[sv['start']]}+{breakpoint_labels[sv['end']]},"
            # Add grammar for inverted duplication
            elif sv['type'] == 'INVDUP':
                grammar += f"{breakpoint_labels[sv['start']]}^+{breakpoint_labels[sv['end']]},"

            elif sv['type'] == 'DEL_DUP_breakpoints_direction_mismatch':
                cpx = True
                grammar += f"DEL_DUP_bp_direction_mismatch_start_{sv['start']}_end_{sv['end']},"
            elif sv['type'] == 'INV_breakpoints_direction_mismatch':
                cpx = True
                grammar += f"INV_bp_direction_mismatch_start_{sv['start']}_end_{sv['end']},"

            # Mark these breakpoints as used
            elif sv['type'] == 'COMPLEX_SV':
                cpx = True
                grammar += f"CPX_chrom{sv['chrom']}_start_{sv['start']}_end_{sv['end']},"
            elif sv['type'] == 'cross_chrom':
                cpx = True
                grammar += f"cross_chrom{sv['chrom']}_start_{sv['start']}_end_{sv['end']},"
            elif sv['type'] == 'event_outside_chr14':
                cpx = True
                grammar += f"outside_chrom14_chrom{sv['chrom']}_start_{sv['start']}_end_{sv['end']},"
            elif sv['type'] == 'DEL_DUP_far_from_breakpoint':
                cpx = True
                grammar += f"DEL_DUP_far_bp_start_{sv['start']}_end_{sv['end']},"
            elif sv['type'] == 'INV_far_from_breakpoint':
                cpx = True
                grammar += f"INV_far_bp_start_{sv['start']}_end_{sv['end']},"
            elif sv['type'] == 'INVDUP_far_from_breakpoint':
                cpx = True
                grammar += f"INVDUP_far_bp_start_{sv['start']}_end_{sv['end']},"
            elif sv['type'] == 'INVDUP_breakpoints_direction_mismatch':
                cpx = True
                grammar += f"INVDUP_bp_direction_mismatch_start_{sv['start']}_end_{sv['end']},"
            elif sv['type'] == 'secondary_alignment':
                cpx = True
                grammar += f"secondary_alignment_chrom{sv['chrom']}_start_{sv['start']}_end_{sv['end']},"
            elif sv['type'] == 'unmapped_read':
                cpx = True
                grammar += f"unmapped_read_chrom{sv['chrom']}_start_{sv['start']}_end_{sv['end']},"

            if cpx:
                loxpsym += ","
            else:
                loxpsym += f"{int(sv['loxpsym']['loxpsym'])}{int(sv['loxpsym']['left_flank'])}{int(sv['loxpsym']['core'])}{int(sv['loxpsym']['right_flank'])},"
         

        if cpx is True:
            grammar_dict_cpx[original_read_name] = grammar
            loxpsym_dict_complex[original_read_name] = loxpsym
        else:
            grammar_dict_clean[original_read_name] = grammar
            loxpsym_dict_clean[original_read_name] = loxpsym

    return grammar_dict_clean, grammar_dict_cpx, loxpsym_dict_clean, loxpsym_dict_complex

def write_vcf(output_file, grammar_dict_clean, grammar_dict_cpx, loxpsym_dict_clean, loxpsym_dict_complex, breakpoints_list, coverage_results):
    def write_vcf_file(output_file, grammar_dict, loxpsym_dict, coverage_results):
        with open(output_file, 'w') as vcf_file:
            # Write VCF headers
            vcf_file.write("##fileformat=VCFv4.2\n")
            vcf_file.write("##INFO=<ID=GRAMMAR,Number=1,Type=String,Description=\"Structural variant grammar notation\">\n")
            vcf_file.write("##INFO=<ID=LOXPSYM=,Number=1,Type=String,Description=\"Loxpsym site location\">\n")
            vcf_file.write("##INFO=<ID=SUPPORT_READS,Number=.,Type=String,Description=\"Supporting read names for SV\">\n")
            vcf_file.write("##INFO=<ID=BREAKPOINTS,Number=.,Type=String,Description=\"List of breakpoints used\">\n")
            vcf_file.write("##INFO=<ID=COVERAGE_DETAILS_FORWARD,Number=.,Type=String,Description=\"Coverage details between breakpoints on forward strand\">\n")
            vcf_file.write("##INFO=<ID=COVERAGE_DETAILS_REVERSE,Number=.,Type=String,Description=\"Coverage details between breakpoints on reverse strand\">\n")
            vcf_file.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            
            breakpoints_str = ",".join(map(str, breakpoints_list))
            
            for original_read_name, grammar in grammar_dict.items():
                
                chrom = "chr14"
                pos = 0  # Placeholder position, adjust as needed
                support_reads = original_read_name
                
                # Get coverage details
                coverage_details_forward = []
                coverage_details_reverse = []
                if original_read_name in coverage_results:
                    for (start_bp, end_bp), coverage in coverage_results[original_read_name]['forward'].items():
                        if pos == 0:
                            pos = start_bp
                        coverage_details_forward.append(f"{start_bp}-{end_bp}:{coverage}")
                    for (start_bp, end_bp), coverage in coverage_results[original_read_name]['reverse'].items():
                        coverage_details_reverse.append(f"{start_bp}-{end_bp}:{coverage}")
                
                coverage_info_forward = "+".join(coverage_details_forward)
                coverage_info_reverse = "+".join(coverage_details_reverse)

                loxpsym = loxpsym_dict[original_read_name]
                
                # Build the INFO field
                info = f"GRAMMAR={grammar};LOXPSYM={loxpsym};SUPPORT_READS={support_reads};BREAKPOINTS={breakpoints_str};COVERAGE_DETAILS_FORWARD={coverage_info_forward};COVERAGE_DETAILS_REVERSE={coverage_info_reverse}"
                
                # Write to VCF
                vcf_file.write(f"{chrom}\t{pos}\t.\tN\t<SV>\t.\tPASS\t{info}\n")

    # Write clean grammar VCF
    write_vcf_file(output_file + ".clean.vcf", grammar_dict_clean, loxpsym_dict_clean, coverage_results)
    # Write complex grammar VCF
    write_vcf_file(output_file + ".cpx.vcf", grammar_dict_cpx, loxpsym_dict_complex, coverage_results)




# Function to load breakpoints from a list or a file
def load_breakpoints(breakpoints_list):
    return sorted(breakpoints_list)

