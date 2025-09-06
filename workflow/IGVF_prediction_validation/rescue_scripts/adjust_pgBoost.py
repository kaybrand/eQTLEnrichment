"""
Add Percentile column to pgBoost predictions
(USE eQTL venv to run!)
July 15, 2025
Kayla Brand
"""

import gzip
import os
import sys
import pandas as pd
from scipy.stats import rankdata
import argparse
import datetime

def calculate_percentiles(series):
    """Calculates percentiles of a Pandas Series using rankdata."""
    ranks = rankdata(series, method='ordinal')  # Or choose a different method if needed
    percentiles = 100 * (ranks - 1) / (len(series) - 1)  # Ensure 0-100 range
    return percentiles

def add_precentile_suffix(input_file_path):
    # Get the directory and base filename
    dir_name = os.path.dirname(input_file_path)
    base_name = os.path.basename(input_file_path)
    
    # Create output filename with '_absolute' appended
    if base_name.endswith("_percentile_percentile.e2g.tsv.gz"):
        # Too long
        output_name = base_name[:-len("_percentile_percentile.e2g.tsv.gz")] + '_percentile.e2g.tsv.gz'
    elif base_name.endswith("_percentile.e2g.tsv.gz"):
        # Just right
        output_name = base_name
    elif base_name.endswith('.e2g.tsv.gz'):
        # Too short
        output_name = base_name[:-len('.e2g.tsv.gz')] + '_percentile.e2g.tsv.gz'
    else:
        output_name = base_name + '_percentile.tsv.gz'
    
    output_file_path = os.path.join(dir_name, output_name)
    return output_file_path

def add_percentile_column(input_file, output_file, sig_fig=5, source_col='Score', percentile_col='Percentile'):
    """
    Adds a 'Percentile' column to a gzipped TSV file containing the percentile
    of each value in the 'Score' column. Uses pd.read_csv for efficiency.

    Args:
        input_file (str): Path to the gzipped TSV input file.
        output_file (str): Chosen name for output file
        Score (str):  Name of column to convert to percentiles
        sig_fig (int): Number of decimal places to keep in percentile score

    Returns:
        output_file_name
    """
    # Check if input path is valid
    if not os.path.isfile(input_file):
        raise FileNotFoundError(f"No such file: {input_file}")

    # Check if this file was created already
    if os.path.exists(output_file) and (os.path.getmtime(output_file) > os.path.getmtime(input_file)):
        print(f"ALREADY EXISTS: {output_file}")
        return output_file

    # Create the file
    comments = []  # Store the comments

    try:
        # Save the comments
        with gzip.open(input_file, 'rt') as f:
            for line in f:
                line = line.strip()
                if line.startswith('#'):
                    comments.append(line + '\n') # Keep the newline
                else:
                    break  # Stop reading comments when data starts
         # Read in the data       
        df = pd.read_csv(input_file, sep='\t', comment='#', compression='gzip')
    except FileNotFoundError:
        raise FileNotFoundError(f"Input file not found: {input_file}")
    except Exception as e:
        raise ValueError(f"Error reading input file: {e}")

    if source_col not in df.columns:
        raise ValueError(f"The input file must have a {source_col} column.")
    
    # if ('Percentile' in df.columns) or ('percentile' in df.columns):
    #     print(f"Percentile column already exists in {input_file}, exiting...")
    #     return input_file

    # Convert 'Score' column to numeric (important for percentile calculation)
    df[source_col] = pd.to_numeric(df[source_col], errors='coerce')

    # Calculate percentiles
    df[percentile_col] = calculate_percentiles(df[source_col])

    # Convert Percentile to string, with the same number of decimal places as the original file
    df[percentile_col] = df[percentile_col].round(sig_fig).astype(str)

    output_file_name = output_file or input_file
    
    # Write the comments back to the beginning of the file
    with gzip.open(output_file_name, 'wt') as f:
        f.writelines(comments) # Write all comments
        df.to_csv(f, sep='\t', index=False, header=True) # Write data

    return output_file_name


def main(args):
    if (args.input_file == None):
            print (parser.usage)
            exit(0)
    else:
        if (args.output_file == None):
            output_file = add_precentile_suffix(args.input_file)
        else:
            output_file = args.output_file

        print(f"Adding Percentile column to pgBoost prediction file: {args.input_file}")
        add_percentile_column(args.input_file, output_file, source_col=args.source_col, percentile_col=args.percentile_col)
        print(f"Find updated file at: {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "Find percentile of pgBoost scores")
    parser.add_argument('-i', '--input', dest='input_file', help = 'path to pgBoost prediction file in IGVF format')
    parser.add_argument('-o', '--output', dest='output_file', help = 'name of validated pgBoost file')
    parser.add_argument('-s', '--source-score-col', dest='source_col', default='Score', help = 'column to take percentile of')
    parser.add_argument('-p', '--percentile-score-col', dest='percentile_col', default='Percentile', help = 'name of column to save percentile in')
    args = parser.parse_args()
    main(args)
    
