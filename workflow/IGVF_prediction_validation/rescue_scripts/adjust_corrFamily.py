"""
Take absolute value of correlation-based scores (Kendall, ArchR, Cicero, Signac)
(Use eQTL venv to run)
July 15, 2025
Kayla Brand
"""
import gzip
import os
from pathlib import Path
import sys
import pandas as pd
# import numpy as np
import re
import argparse
import datetime

def add_absolute_suffix(input_file_path):
    # Get the directory and base filename
    dir_name = os.path.dirname(input_file_path)
    base_name = os.path.basename(input_file_path)
    
    # Create output filename with '_absolute' appended
    if base_name.endswith('_absolute.e2g.tsv.gz'):
        return input_file_path
    elif base_name.endswith('.e2g.tsv.gz'):
        output_name = base_name[:-len('.e2g.tsv.gz')] + '_absolute.e2g.tsv.gz'
        output_file_path = os.path.join(dir_name, output_name)
        return output_file_path
    else:
        raise ValueError(f"Unexpected file name: {base_name}")

def take_absolute_value_of_corr(input_file_path, output_file_path, corr_score_col="Score", abs_score_col="abs_Score"):
    """
    Process a gzipped TSV file to:
    1. Extract and preserve comment lines (starting with #) as key-value pairs
    2. Verify that "# ScoreType: divergent" or "# ScoreType: Divergent" exists
    3. Update the ScoreType line to "# ScoreType: abs(correlation)"
    4. Process the data, taking the absolute value of the 'Score' column
    5. Save the result to a new file with '_absolute' appended to the name
    
    Args:
        input_file_path (str): Path to the input .tsv.gz file
        
    Returns:
        str: Path to the output file
    """
    if not os.path.isfile(input_file_path):
        raise FileNotFoundError(f"No such file {input_file_path}")

    # If this file was already produced and was modified more recently than the input, exit
    if os.path.exists(output_file_path) and (os.path.getmtime(output_file_path) > os.path.getmtime(input_file_path)):
        print(f"ALREADY EXISTS: {output_file_path}")
        return output_file_path
    
    # Storage for comment lines
    comments = []
    
    # Read the input file
    with gzip.open(input_file_path, 'rt') as f:
        # Process comments and header
        for line in f:
            if line.startswith('#'):
                comments.append(line.strip())
            else:
                break
    
    # Check for required ScoreType comment and update it
    score_type_found = False
    modified_comments = []
    
    for comment in comments:
        if re.match(r'# ScoreType:\s*[dD]ivergent', comment):
            modified_comments.append("# ScoreType: abs(correlation)")
            score_type_found = True
        elif ("ScoreType" in comment):
            score_type = comment.split(" ")[-1] 
            modified_comments.append(f"# ScoreType: abs({score_type})")
        else:
            modified_comments.append(comment)
    
    if not score_type_found:
        print("Error: Missing required ScoreType comment with 'divergent' or 'Divergent' value")
        print("Ignoring issue and forging blithely ahead...")
    
    # Read the data
    df = pd.read_csv(input_file_path, 
                     sep='\t', 
                     comment='#', 
                     header=0)
    
    # Take absolute value of Score column
    if corr_score_col in df.columns:
        df[abs_score_col] = df[corr_score_col].abs()
        print(f"Saved the absolute value of the correlation score as the {abs_score_col} column")
    else:
        raise ValueError(f"Error: {corr_score_col} column not found in the data")
    
    # Write out the processed file
    with gzip.open(output_file_path, 'wt') as f:
        # Write modified comments
        for comment in modified_comments:
            f.write(comment + '\n')
        
        # Write the data
        df.to_csv(f, sep='\t', index=False)
    
    return output_file_path

def main(input_file, output_file, corr_score_col, abs_score_col):
    print(f"Taking absolute value of correlation in {corr_score_col} column: {input_file}")
    take_absolute_value_of_corr(input_file, output_file, corr_score_col, abs_score_col)
    print(f"Find updated file at: {output_file}")

if __name__ == "__main__":  
    parser = argparse.ArgumentParser(description = "Adjust a correlation score prediction file for eQTL pipeline")
    parser.add_argument("input_file", help="path to prediction file")
    parser.add_argument("output_file", help="path to be given to output file")
    parser.add_argument("-s", "--source-score-col", action='store', dest='source_col', help="column to take absolute value of")
    parser.add_argument("-t", "--target-score-col", action='store', dest='target_col', help='name of column to be created')
    args = parser.parse_args()

    # Pull file paths from Namespace
    main(args.input_file, args.output_file, args.source_col, args.target_col)
    
