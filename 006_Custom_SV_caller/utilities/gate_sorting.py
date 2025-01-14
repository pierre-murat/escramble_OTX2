import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import argparse
import os
import logging



def load_and_label_data(file_path, label):
    data = pd.read_csv(file_path)
    data['source'] = label
    return data

def augment_with_replicates(data, replicates_lookup):
    replicates_df = pd.read_csv(replicates_lookup, sep='\t')
    
    logging.info("original data shape: %s", data.shape)
    # Select only the columns needed
    logging.info("original replicates size: %s", replicates_df.shape)
    replicates_df = replicates_df[['read.id', 'replicate']]
    logging.info("number of duplicated rows: %s", replicates_df['read.id'].duplicated().sum())  # Counts duplicated `read.id` entries

    # Identify duplicated rows in replicates_df
    duplicated_replicates = replicates_df[replicates_df['read.id'].duplicated(keep=False)]
    logging.info("duplicated_replicates: %s", duplicated_replicates.shape)

    # Drop duplicated rows from replicates_df
    replicates_df = replicates_df.drop_duplicates(subset='read.id', keep=False)
    logging.info("replicates_df after dropping duplicates: %s", replicates_df.shape)

    logging.info("data shape before dropping duplicates: %s", data.shape)

    # Drop corresponding rows in data

    rows_before = data.shape[0]
    data = data[~data['support_reads'].isin(duplicated_replicates['read.id'])]
    rows_after = data.shape[0]
    logging.info("data after dropping duplicates: %s", data.shape)
    logging.info("removed %s reads due to duplicate reads", rows_before - rows_after)

    # Perform the merge and ensure it adds only the 'replicate' column where matches exist
    data = data.merge(replicates_df, left_on='support_reads', right_on='read.id', how='left')
    logging.info("data shape after merge: %s", data.shape)
    
    # Drop the redundant 'read.id' column from the merge
    #data = data.drop(columns=['read.id'])
    
    return data


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot percentage of read counts for each vector_assignment from CSV files.")
    parser.add_argument("--very_dim", required=True, help="Path to very_dim CSV file")
    parser.add_argument("--dim", required=True, help="Path to dim CSV file")
    parser.add_argument("--bright", required=True, help="Path to bright CSV file")
    parser.add_argument("--very_bright", required=True, help="Path to very_bright CSV file")
    parser.add_argument("--output_dir", required=True, help="Directory to save the output plots")
    parser.add_argument("--replicates_lookup", required=True, help="Path to replicates lookup TSV file")

    args = parser.parse_args()
    # Set up logging
    logging.basicConfig(filename=args.output_dir + '/gate_analysis_stats.log', level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

    # Load data with labels
    very_dim = load_and_label_data(args.very_dim, 'very_dim')
    logging.info("very_dim: %s", very_dim)
    dim = load_and_label_data(args.dim, 'dim')
    logging.info("dim: %s", dim)
    bright = load_and_label_data(args.bright, 'bright')
    logging.info("bright: %s", bright)
    very_bright = load_and_label_data(args.very_bright, 'very_bright')
    logging.info("very_bright: %s", very_bright)

    # Concatenate all data into a single DataFrame
    data = pd.concat([very_dim, dim, bright, very_bright])
    logging.info("data: %s, %s", data, data.shape)
    logging.info("initiial # of reads: %s", data.shape[0])


    # Augment data with replicates
    data = augment_with_replicates(data, args.replicates_lookup)
    logging.info("augmented data: %s, %s", data, data.shape)

    data['normalized_grammar'] = data['normalized_grammar'].fillna('WT') # WT entries have NaN values as grammar

    data.to_csv(os.path.join(args.output_dir, 'augmented_data.csv'), index=False)
    grouped_df = data.groupby(['vector_assignment', 'normalized_grammar', 'source', 'replicate']).agg(
    read_count=('read.id', 'count'),
    read_support=('read.id', lambda x: list(x))).reset_index()

    grouped_df['percentage'] = grouped_df.groupby(['source', 'replicate'])['read_count'].transform(lambda x: (x / x.sum()))

    pivot_df = grouped_df.pivot(
    index=['vector_assignment', 'normalized_grammar'],
    columns=['source', 'replicate'],
    values=['percentage', 'read_count', 'read_support'])

    # Flatten the MultiIndex columns
    pivot_df.columns = ['_'.join(map(str, col)).strip() for col in pivot_df.columns.values]
    pivot_df = pivot_df.reset_index()
    pivot_df = pivot_df.fillna(0)
    read_count_columns = [col for col in pivot_df.columns if 'read_count' in col]

    pivot_df['total_read_count'] = pivot_df[read_count_columns].sum(axis=1)
    pivot_df

    pivot_df['mean_expression_score_1'] = (
        -1 * pivot_df['percentage_very_dim_1'] +
        -0.5 * pivot_df['percentage_dim_1'] +
        0.5 * pivot_df['percentage_bright_1'] +
        1 * pivot_df['percentage_very_bright_1']
    ) / (
        pivot_df['percentage_very_dim_1'] +
        pivot_df['percentage_dim_1'] +
        pivot_df['percentage_bright_1'] +
        pivot_df['percentage_very_bright_1']
    )

    pivot_df['mean_expression_score_2'] = (
        -1 * pivot_df['percentage_very_dim_2'] +
        -0.5 * pivot_df['percentage_dim_2'] +
        0.5 * pivot_df['percentage_bright_2'] +
        1 * pivot_df['percentage_very_bright_2']
    ) / (
        pivot_df['percentage_very_dim_2'] +
        pivot_df['percentage_dim_2'] +
        pivot_df['percentage_bright_2'] +
        pivot_df['percentage_very_bright_2']
    )

    logging.info("WT: %s", pivot_df[(pivot_df["normalized_grammar"] == "WT") & (pivot_df["vector_assignment"] == "((1, 1, 1, 1, 1, 1), (0, 0, 0, 0, 0, 0))")])
    WT_mean_expression_score_1 = pivot_df.loc[(pivot_df["normalized_grammar"] == "WT") & (pivot_df["vector_assignment"] == "((1, 1, 1, 1, 1, 1), (0, 0, 0, 0, 0, 0))"), "mean_expression_score_1"].values[0]
    WT_mean_expression_score_2 = pivot_df.loc[(pivot_df["normalized_grammar"] == "WT") & (pivot_df["vector_assignment"] == "((1, 1, 1, 1, 1, 1), (0, 0, 0, 0, 0, 0))"), "mean_expression_score_2"].values[0]

    logging.info("full del: %s", pivot_df[(pivot_df["normalized_grammar"] == "A_G") & (pivot_df["vector_assignment"] == "((0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0))")])
    full_del_mean_expression_score_1 = pivot_df.loc[(pivot_df["normalized_grammar"] == "A_G") & (pivot_df["vector_assignment"] == "((0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0))"), "mean_expression_score_1"].values[0]
    full_del_mean_expression_score_2 = pivot_df.loc[(pivot_df["normalized_grammar"] == "A_G") & (pivot_df["vector_assignment"] == "((0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0))"), "mean_expression_score_2"].values[0]
    
    pivot_df['mean_expression_score_1_scaled'] = ((pivot_df['mean_expression_score_1'] - full_del_mean_expression_score_1) /  (WT_mean_expression_score_1 - full_del_mean_expression_score_1)) 
    pivot_df['mean_expression_score_2_scaled'] = ((pivot_df['mean_expression_score_2'] -  full_del_mean_expression_score_2) /  (WT_mean_expression_score_2 - full_del_mean_expression_score_2))
    pivot_df['mean_expression_score_avg'] = pivot_df[['mean_expression_score_1_scaled', 'mean_expression_score_2_scaled']].mean(axis=1)

    pivot_df["event_name"] = pivot_df["vector_assignment"] + " " + pivot_df["normalized_grammar"]
    pivot_df["forward_coverage"] = pivot_df["vector_assignment"].apply(lambda x: eval(x)[0])
    pivot_df["reverse_coverage"] = pivot_df["vector_assignment"].apply(lambda x: eval(x)[1])
    pivot_df["sum_coverage"] = pivot_df.apply(lambda row: tuple(map(sum, zip(row["forward_coverage"], row["reverse_coverage"]))), axis=1)

    pivot_df = pivot_df.sort_values(by='total_read_count', ascending=False)
    
    pivot_df.to_csv(os.path.join(args.output_dir, 'grouped_data_by_architectures.csv'), index=False)
    pivot_df.drop(columns=[col for col in pivot_df.columns if 'read_support' in col]).to_csv(os.path.join(args.output_dir, 'grouped_data_by_architectures_remove_read_ids.csv'), index=False)

        # Create a mask where read_count is greater than 0
    read_count_rep_1_columns = [col for col in pivot_df.columns if col.startswith('read_count') and col.endswith('_1')]
    read_count_rep_2_columns = [col for col in pivot_df.columns if col.startswith('read_count') and col.endswith('_2')]

    rows_before = pivot_df.shape[0]
    non_zero_read_count_mask = (((pivot_df[read_count_rep_1_columns] > 0).sum(axis=1) >= 2) | ((pivot_df[read_count_rep_2_columns] > 0).sum(axis=1) >= 2))

    # Filter out rows where the condition is not met
    logging.info("number of cases with < 2 columns with non-zero coverage: %s, %s", pivot_df[~non_zero_read_count_mask], pivot_df[~non_zero_read_count_mask].shape)
    pivot_df = pivot_df[non_zero_read_count_mask]
    rows_after = pivot_df.shape[0]
    logging.info("filtered based on < 2 gates (across replicates) with non-zero read counts: %s", rows_before - rows_after)
    pivot_df.to_csv(os.path.join(args.output_dir, 'grouped_data_by_architectures_read_across_gates.csv'), index=False)
    pivot_df.drop(columns=[col for col in pivot_df.columns if 'read_support' in col]).to_csv(os.path.join(args.output_dir, 'grouped_data_by_architectures_read_across_gates_remove_read_ids.csv'), index=False)

    rows_before = pivot_df.shape[0]
    pivot_df = pivot_df[(abs(pivot_df["mean_expression_score_1_scaled"] - pivot_df["mean_expression_score_2_scaled"]) < 0.25)]
    rows_after = pivot_df.shape[0]
    logging.info("filtered based on replicates concordance 0.25: %s", rows_before - rows_after)

    pivot_df.to_csv(os.path.join(args.output_dir, 'grouped_data_by_architectures_replicates_0.25_read_across_gates.csv'), index=False)
    pivot_df.drop(columns=[col for col in pivot_df.columns if 'read_support' in col]).to_csv(os.path.join(args.output_dir, 'grouped_data_by_architectures_replicates_0.25_read_across_gates_remove_read_ids.csv'), index=False)

