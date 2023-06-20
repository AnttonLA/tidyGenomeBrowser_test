import argparse
import os
import pandas as pd

"""
T&his script takes a .gwas file and creates a new txt file that contains only its chromosome and minumum and maximum
positions.
"""


def extract_window(raw_args=None):
    parser = argparse.ArgumentParser(description='Extract the genomic window from a gwas file.')
    parser.add_argument('gwas_file', type=str, help='The gwas file to be curated.')
    parser.add_argument('-c', '--chromosome', type=str, default='chromosome',
                        help='Name of column containing the chromosome. Default: chromosome')
    parser.add_argument('-p', '--position', type=str, default='position',
                        help='Name of column containing the position. Default: position')
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Name of the output file.')

    args = parser.parse_args(raw_args)

    # Check that we can find the gwas file
    if not os.path.exists(args.gwas_file):
        raise ValueError(f"GWAS file {args.gwas_file} not found.")

    # Read the gwas file
    df = pd.read_csv(args.gwas_file, sep='\t', header=0)

    # Make sure the chromosome and position columns exist
    user_columns = [args.chromosome, args.position]
    for column in user_columns:
        if column not in df.columns:
            raise ValueError(f"Column {column} not found in the gwas file.")

    # Get the values
    chromosome = df[args.chromosome].iloc[0]
    start = df[args.position].min()
    end = df[args.position].max()

    # Write the output file
    with open(args.output_file, 'w') as f:
        f.write(f"{chromosome}\t{start}\t{end}")


if __name__ == '__main__':
    extract_window()
