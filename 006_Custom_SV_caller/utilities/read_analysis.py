import subprocess
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import re
import argparse
import os
import logging



parser = argparse.ArgumentParser(description="generate various plots related to a computed VCF file")
parser.add_argument("clean_vcf_file", type=str, help="Path to clean events VCF file")
parser.add_argument("plot_directory", type=str, nargs='?', help="Path to the directory where plots will be saved", default=None)
args = parser.parse_args()

# If plot_directory is not supplied, create a "plots" folder in the same directory as clean_vcf_file
if not args.plot_directory:
    plot_directory = os.path.join(os.path.dirname(args.clean_vcf_file), "plots")
    os.makedirs(plot_directory, exist_ok=True)
else:
    plot_directory = args.plot_directory
    
# Parse the arguments
args = parser.parse_args()

# Set up logging
logging.basicConfig(filename= plot_directory + '/read_analysis_stats.log', level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def parse_vcf(vcf_file):
    data = []

    with open(vcf_file, 'r') as file:
        for line in file:
            if line.startswith("#"):
                continue

            columns = line.strip().split('\t')
            info = columns[7]
            # Extract GRAMMAR, BREAKPOINTS, COVERAGE_DETAILS_FORWARD, and COVERAGE_DETAILS_REVERSE
            grammar_match = re.search(r'GRAMMAR=([^;]*)', info)
            grammar = grammar_match.group(1) if grammar_match else ""
            loxpsym_match = re.search(r'LOXPSYM=([^;]*)', info)
            loxpsym = loxpsym_match.group(1) if loxpsym_match else ""

            breakpoints = re.search(r'BREAKPOINTS=([^;]+)', info).group(1)
            reads = re.search(r'SUPPORT_READS=([^;]+)', info).group(1)


            coverage_details_forward = re.search(r'COVERAGE_DETAILS_FORWARD=([^;]+)', info).group(1)
            coverage_details_reverse = re.search(r'COVERAGE_DETAILS_REVERSE=([^;]+)', info).group(1)

            # Parse COVERAGE_DETAILS_FORWARD and COVERAGE_DETAILS_REVERSE into dictionaries
            coverage_dict_forward = {}
            coverage_dict_reverse = {}

            for pair in coverage_details_forward.split('+'):
                range_, coverage = pair.split(':')
                coverage_dict_forward[range_] = float(coverage)

            for pair in coverage_details_reverse.split('+'):
                range_, coverage = pair.split(':')
                coverage_dict_reverse[range_] = float(coverage)

            data.append({
                "chromosome": columns[0],
                "position": columns[1],
                "grammar": grammar,
                "loxpsym": loxpsym,
                "breakpoints": breakpoints,
                "support_reads": reads,
                "coverage_details_forward": coverage_dict_forward,
                "coverage_details_reverse": coverage_dict_reverse
            })

    df = pd.DataFrame(data)
    return df


def assign_coverage_vectors(df):
    def assign_vector(coverage_dict, breakpoints):
        vector = []
        for i in range(len(breakpoints) - 1):
            range_key = f"{breakpoints[i]}-{breakpoints[i + 1]}"
            coverage = coverage_dict.get(range_key, 0.0)
            assigned_value = round(coverage / 3.0)
            vector.append(assigned_value)
        return vector
    def assign_vector_raw(coverage_dict, breakpoints):
        vector = []
        for i in range(len(breakpoints) - 1):
            range_key = f"{breakpoints[i]}-{breakpoints[i + 1]}"
            coverage = coverage_dict.get(range_key, 0.0)
            assigned_value = coverage
            vector.append(assigned_value)
        return vector
    
    df['forward_vector'] = df.apply(lambda row: assign_vector(row['coverage_details_forward'], row['breakpoints'].split(',')), axis=1)
    df['reverse_vector'] = df.apply(lambda row: assign_vector(row['coverage_details_reverse'], row['breakpoints'].split(',')), axis=1)
    df['forward_vector_raw'] = df.apply(lambda row: assign_vector_raw(row['coverage_details_forward'], row['breakpoints'].split(',')), axis=1)
    df['reverse_vector_raw'] = df.apply(lambda row: assign_vector_raw(row['coverage_details_reverse'], row['breakpoints'].split(',')), axis=1)

def coverage_analysis(df):
    logging.info(df[['forward_vector', 'reverse_vector']])

    # Filter rows where forward strand coverage is greater than 1 between the first 2 breakpoints and the last 2 breakpoints
    logging.info("filtered rows based on coverage anchoring before: %s", df.shape)
    rows_before = df.shape[0]
    df_filtered = df[df.apply(lambda row: row['forward_vector_raw'][0] > 0 and row['forward_vector_raw'][-1] > 0, axis=1)]
    rows_after = df_filtered.shape[0]
    logging.info("filtered rows based on coverage anchoring after: %s", df_filtered.shape)
    logging.info("removed %s reads due to coverage anchoring, reads leftover %s", rows_before - rows_after, rows_after)

    # Remove the first 2 breakpoints and last 2 breakpoints from the vectors for further analysis
    df_filtered.loc[:, 'forward_vector'] = df_filtered['forward_vector'].apply(lambda x: x[1:-1])
    df_filtered.loc[:, 'reverse_vector'] = df_filtered['reverse_vector'].apply(lambda x: x[1:-1])

    # Combine forward and reverse vectors into a single tuple and convert to string for easier handling
    df_filtered['vector_assignment'] = df_filtered.apply(lambda row: (tuple(row['forward_vector']), tuple(row['reverse_vector'])), axis=1)
    df_filtered['vector_assignment'] = df_filtered['vector_assignment'].astype(str)  # Convert tuples to strings
    logging.info(df_filtered)

    # Count the occurrences of each unique vector assignment
    vector_counts = df_filtered['vector_assignment'].value_counts().reset_index()
    vector_counts.columns = ['vector_assignment', 'count']
    logging.info(vector_counts)

    # Plot the distribution of vector assignments on a logarithmic scale
    plt.figure(figsize=(12, 6))
    sns.barplot(x='vector_assignment', y='count', data=vector_counts)
    plt.yscale('log')  # Set the y-axis to log scale
    plt.xlabel('Vector Assignment (Forward, Reverse)')
    plt.ylabel('Number of Reads (log scale)')
    plt.title('Distribution of Reads with Each Type of Vector Assignment')
    # Save the vector counts to a CSV file
    vector_counts.to_csv(plot_directory + "/vector_assignment_counts.csv", index=False)
    plt.xticks(rotation=90)
    plt.subplots_adjust(bottom=0.3)  # Adjust bottom margin
    plt.savefig(plot_directory + "/vector_assignment_distribution.png")
    coverage_boxplots(df_filtered, strand_key='coverage_details_forward')
    coverage_boxplots(df_filtered, strand_key='coverage_details_reverse')
    return df_filtered

def coverage_boxplots(df, strand_key='coverage_details_forward'):
    coverage_data = []

    for index, row in df.iterrows():
        for range_, coverage in row[strand_key].items():
            coverage_data.append({'range': range_, 'coverage': coverage})

    # Create a DataFrame from the coverage data
    coverage_df = pd.DataFrame(coverage_data)

    # Plot the distribution of coverage values
    plt.figure(figsize=(12, 6))
    sns.boxplot(x='range', y='coverage', data=coverage_df)
    plt.xlabel('Breakpoint Range')
    plt.ylabel('Coverage')
    plt.title('Distribution of Coverage Values Between Breakpoints: ' + f' ({strand_key})')
    plt.xticks(rotation=90)
    plt.subplots_adjust(bottom=0.3)  # Adjust bottom margin
    plt.savefig(plot_directory + "/" + strand_key + ".png")

    
# Example usage
vcf_file = args.clean_vcf_file
df = parse_vcf(vcf_file)
logging.info("initial number of reads: %s", df.shape[0])


assign_coverage_vectors(df)
df_filtered = coverage_analysis(df)

# Function to normalize grammar patterns by sorting elements and treating equivalent patterns as the same
def normalize_grammar(grammar, loxpsym):
    parts = grammar.split(",")  # Split by comma
    loxpsym_parts = loxpsym.split(",")
    normalized_parts = []
    for i, part in enumerate(parts):
        if not part:
            continue
        presence_loxpsym_site = loxpsym_parts[i]
        if not (presence_loxpsym_site[0] == "1" or presence_loxpsym_site[1] == "1" or presence_loxpsym_site[3] == "1"): 
            return "loxpsym_unsupported"
            
        if "^+" in part:
            elements = part.split("^+")
            normalized_parts.append("^+".join(sorted(elements)))
        elif "+" in part:
            elements = part.split("+")
            normalized_parts.append("+".join(sorted(elements)))
        elif "^" in part:
            elements = part.split("^")
            normalized_parts.append("^".join(sorted(elements)))
        elif "_" in part:
            elements = part.split("_")
            normalized_parts.append("_".join(sorted(elements)))
        else:
            normalized_parts.append(part.strip())

     # Collapse consecutive identical elements
    collapsed_parts = []
    for part in normalized_parts:
        if not collapsed_parts or collapsed_parts[-1] != part:
            collapsed_parts.append(part)
    normalized_parts = collapsed_parts

    return ",".join((normalized_parts))


# Apply normalization to grammar column
df_filtered["normalized_grammar"] = df_filtered.apply(lambda row: normalize_grammar(row["grammar"], row["loxpsym"]), axis=1)
df_filtered['normalized_grammar'] = df_filtered['normalized_grammar'].fillna('WT') # WT entries have NaN values as grammar


# Drop entries with "normalized_grammar" field == "loxpsym_unsupported"
logging.info("before dropping unsupported: %s", df_filtered.shape)
rows_before = df_filtered.shape[0]
logging.info("reads being dropped due to loxpsym support: %s", df_filtered[df_filtered["normalized_grammar"] == "loxpsym_unsupported"])
df_filtered = df_filtered[df_filtered["normalized_grammar"] != "loxpsym_unsupported"]
rows_after = df_filtered.shape[0]
logging.info("after dropping loxpsym unsupported: %s", df_filtered.shape)
logging.info("removed %s reads due to loxpsym support", rows_before - rows_after)

logging.info("final number of reads: %s", rows_after)

# Save the filtered DataFrame to a CSV file for further analysis in a Jupyter notebook
df_filtered.to_csv(plot_directory + "/filtered_data.csv", index=False)