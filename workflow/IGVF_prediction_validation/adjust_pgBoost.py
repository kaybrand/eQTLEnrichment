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
from optparse import OptionParser

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

def add_percentile_column(input_file, output_file, sig_fig=5, score_col='Score'):
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
    if os.path.exists(output_file):
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

    if score_col not in df.columns:
        raise ValueError(f"The input file must have a {score_col} column.")
    
    # if ('Percentile' in df.columns) or ('percentile' in df.columns):
    #     print(f"Percentile column already exists in {input_file}, exiting...")
    #     return input_file

    # Convert 'Score' column to numeric (important for percentile calculation)
    df[score_col] = pd.to_numeric(df[score_col], errors='coerce')

    # Calculate percentiles
    df['Percentile'] = calculate_percentiles(df[score_col])

    # Convert Percentile to string, with the same number of decimal places as the original file
    df['Percentile'] = df['Percentile'].round(sig_fig).astype(str)

    output_file_name = output_file or input_file
    
    # Write the comments back to the beginning of the file
    with gzip.open(output_file_name, 'wt') as f:
        f.writelines(comments) # Write all comments
        df.to_csv(f, sep='\t', index=False, header=True) # Write data

    return output_file_name


def main():
    parser = OptionParser()
    parser.add_option('-i', '--input', dest='input_file', type='string', help = 'path to pgBoost prediction file in IGVF format')
    parser.add_option('-o', '--output', dest='output_file', type='string', help = 'name of validated pgBoost file')
    parser.add_option('-s', '--score-col', dest='raw_score_col', type='string', default='Score', help = 'column to take percentile of')
    # parser.add_option('-p', '--percentile-name', dest='percentile_col', type='string', default='Percentile', help = 'name of column to save percentile in')
    (options, args) = parser.parse_args()
    if (options.input_file == None):
            print (parser.usage)
            exit(0)
    else:
        if (options.output_file == None):
            output_file = add_precentile_suffix(options.input_file)
        else:
            output_file = options.output_file
        print(f"Adding Percentile column to pgBoost prediction file: {options.input_file}")
        add_percentile_column(options.input_file, output_file, score_col=options.raw_score_col)
        print(f"Find updated file at: {output_file}")

if __name__ == "__main__":
    main()
    # if len(sys.argv) != 2:
    #     print("Usage: python script.py <input_file.tsv.gz> <output_file.tsv.gz>")
    #     sys.exit(1)

    # input_file = sys.argv[1]
    # output_file = sys.argv[2]
    
