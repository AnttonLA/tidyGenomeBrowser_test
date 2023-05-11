import argparse
import os
import sys
import pandas as pd

"""
This script is used to curate the gwas file before the plotting.
See format: https://software.broadinstitute.org/software/igv/GWAS

The script can apply a phenotype filter, a p-value filter and a region filter to the given file.
It can also simplify the names of the phenotypes to be more human-readable.

The phenotype filter and the mapping make use of external text files. If non are provided but a 'resources' directory is
given, the script will create them. The user can then edit those files and use them in the next run of the script.
"""

parser = argparse.ArgumentParser(description='Curate the gwas file before plotting.')
parser.add_argument('gwas_file', type=str, help='The gwas file to be curated.')
parser.add_argument('-f', '--phenotype_file', type=str,
                    help='File listing the phenotypes to be included in the output. If not provided, all phenotypes '
                         'will be included. Note that these should be the original names, not the simplified ones.')
parser.add_argument('-fm', '--phenotype_map', type=str,
                    help='File containing the mapping between the original phenotype names and the simplified ones. '
                         'If non is given but a resources directory is specified , a new file with the original names '
                         'repeated in two columns will be created.')
parser.add_argument('-fc', '--phenotype_column', type=str, default="phenotype",
                    help='Name of column containing the phenotype. Default: phenotype')
parser.add_argument('-r', '--resources_directory', type=str,
                    help='Resources directory. Additional files created by the script are stored here.')
parser.add_argument('-p', '--p_value', type=float, default=None,
                    help='P-value threshold. Entries with p-values larger than the threshold will be dropped. '
                         'Default: None')
parser.add_argument('-pc', '--p_value_column', type=str, default="pval",
                    help='Name of column containing the p-value. Default: pval')
parser.add_argument('-s', '--new_start', type=int, default=None,
                    help='New leftmost limit for the genomic window. Positions smaller than the limit will be '
                         'dropped. Default: None')
parser.add_argument('-e', '--new_end', type=int, default=None,
                    help='New rightmost limit for the genomic window. Positions larger than the limit will be '
                         'dropped. Default: None')
parser.add_argument('-bpc', '--position_column', type=str, default="position",
                    help='Name of column containing the end position. Default: end')
parser.add_argument('-o', '--output_file', type=str, required=True,
                    help='Name of the curated gwas file. It is recommended to use the same name as the input file, '
                         'with the suffix "_curated".')

args = parser.parse_args()

# Create output directory if it does not exist
if not os.path.exists(args.output_directory):
    os.makedirs(args.output_directory)

# If a resources directory was specified but does not exist, create resources directory
if args.resources_directory is not None:
    if not os.path.exists(args.resources_directory):
        os.makedirs(args.resources_directory)

# Read the gwas file
df = pd.read_csv(args.gwas_file, sep='\t', header=0)

# Check that we can find all the provided columns
user_columns = [args.phenotype_column, args.p_value_column, args.position_column]
for column in user_columns:
    if column not in df.columns:
        raise ValueError("Column {} not found in the gwas file.".format(column))

# Apply phenotype name map if a file was given
if args.phenotype_map is not None:
    if not os.path.exists(args.phenotype_map):  # Check that the file exists
        raise ValueError("Phenotype map file {} not found.".format(args.phenotype_map))

    # Read the phenotype map file
    phenotype_map_df = pd.read_csv(args.phenotype_map, sep='\t', header=None)
    if phenotype_map_df.empty:  # Check that it is not empty
        raise ValueError("Phenotype map file {} is empty.".format(args.phenotype_map))
    if phenotype_map_df.shape[1] != 2:  # Check that it has two columns
        raise ValueError(f"Phenotype map file {args.phenotype_map} does not have two columns. Please ensure that the "
                         f"map file is a tab-separated file with two columns, the first one containing the original "
                         f"phenotype names and the second one containing the simplified names.")

    # Check that the phenotypes are present in the gwas file
    for phenotype in phenotype_map_df[0].tolist():
        if phenotype not in df[args.phenotype_column].tolist():
            sys.stdout.write(f"WARNING Phenotype '{phenotype}' in the phenotype map file {args.phenotype_map} could "
                             f"not be found in the gwas file. It will be ignored.\n")
            phenotype_map_df = phenotype_map_df[phenotype_map_df[0] != phenotype]
    # Make a dictionary with the mapping
    phenotype_map = dict(zip(phenotype_map_df[0].tolist(), phenotype_map_df[1].tolist()))
    # Create a new column with the simplified names
    df['simplified_phenotype'] = df[args.phenotype_column].map(phenotype_map)
    # Fill the simplified names with the original names if no mapping was found
    df['simplified_phenotype'] = df['simplified_phenotype'].fillna(df[args.phenotype_column])

else:  # if no map file given but resources directory was specified, create a new map file with the original names
    if args.resources_directory is not None:
        phenotype_map_df = pd.DataFrame(df[args.phenotype_column].unique(), columns=['original'])
        phenotype_map_df['simplified'] = phenotype_map_df['original']
        phenotype_map_df.to_csv(os.path.join(args.resources_directory,
                                             f"{os.path.basename(args.gwas_file)}_phenotype_map.tsv"),
                                sep='\t', index=False, header=False)

# Apply phenotype filter if a file was given
if args.phenotype_file is not None:
    if not os.path.exists(args.phenotype_file):  # Check that the file exists
        raise ValueError("Phenotype file {} not found.".format(args.phenotype_file))

    # Read the phenotype file
    phenotype_df = pd.read_csv(args.phenotype_file, sep='\t', header=None)
    if phenotype_df.empty:  # Check that it is not empty
        raise ValueError("Phenotype file {} is empty.".format(args.phenotype_file))
    if phenotype_df.shape[1] != 1:  # Check that it has only one column
        raise ValueError(f"Phenotype file {args.phenotype_file} has more than one column. '"
                         f"Please provide a bare list without tab characters.")

    # Check if all the phenotypes in the file are present in the gwas file
    phenotype_list = phenotype_df[0].tolist()
    for phenotype in phenotype_list:
        if phenotype not in df[args.phenotype_column].tolist():
            sys.stdout.write(f"WARNING: Phenotype '{phenotype}' not found in the gwas file. It will be ignored.\n")
            phenotype_list.remove(phenotype)

    df = df[df[args.phenotype_column].isin(phenotype_list)]

else:  # No phenotype file was given. We will not filter, but we will create a new file with the phenotypes
    if args.resources_directory is not None:
        phenotype_df = df[args.phenotype_column].drop_duplicates().to_frame()
        phenotype_df.to_csv(os.path.join(args.resources_directory,
                                         f"{os.path.basename(args.gwas_file)}_phenotype_list.txt"),
                            sep='\t', index=False, header=False)

# Apply p-value filter if a value was given
if args.p_value is not None:
    df = df[df[args.p_value_column] <= args.p_value]

# Apply region filter if values were given
if args.new_start is not None:
    df = df[df[args.position_column] >= args.new_start]
if args.new_end is not None:
    df = df[df[args.position_column] <= args.new_end]

df.to_csv(args.output_file, sep='\t', index=False, header=True)
