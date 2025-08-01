"""
Kayla Brand
July 15, 2025
Calculate boolean FDR > 0.001, Spearman corr > 0.1, z-score > 0.5 for SCARlink score
"""

import pandas as pd
import gzip
from pathlib import Path
import sys

def process_e2g_file(input_file_path, 
                    score_col='Score', 
                    fdr_col='Score_fdr', 
                    rho_col='Score_rho',
                    score_threshold=0.5,
                    fdr_threshold=0.05,
                    rho_threshold=0.5,
                    chunk_size=10000):
    """
    Process e2g file to add boolean score column based on thresholds.
    
    Parameters:
    -----------
    input_file_path : str
        Path to input e2g.tsv.gz file
    score_col : str
        Name of the score column (default: 'Score')
    fdr_col : str
        Name of the FDR column (default: 'Score_fdr')
    rho_col : str
        Name of the rho column (default: 'Score_rho')
    score_threshold : float
        Threshold for score column (default: 0.5)
    fdr_threshold : float
        Threshold for FDR column (default: 0.05)
    rho_threshold : float
        Threshold for rho column (default: 0.5)
    chunk_size : int
        Number of rows to process at a time (default: 10000)
    
    Returns:
    --------
    str : Path to the output file
    """
    
    # Convert to Path object for easier manipulation
    input_path = Path(input_file_path)
    
    # Check if file exists
    if not input_path.exists():
        raise FileNotFoundError(f"Input file does not exist: {input_file_path}")
    
    # Check if already processed (filename ends with _boolean_score.e2g.tsv.gz)
    if input_path.name.endswith('_boolean_score.e2g.tsv.gz'):
        print(f"File already has boolean score suffix: {input_file_path}")
        return str(input_path)
    
    # Create output file path
    if input_path.name.endswith('.e2g.tsv.gz'):
        output_name = input_path.name.replace('.e2g.tsv.gz', '_boolean_score.e2g.tsv.gz')
    else:
        raise ValueError(f"Input file must end with '.e2g.tsv.gz', got: {input_path.name}")
    
    output_path = input_path.parent / output_name
    
    # Read header comments and first few rows to check structure
    header_comments = []
    
    print(f"Reading file: {input_file_path}")
    
    # Read header comments
    with gzip.open(input_path, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                header_comments.append(line.rstrip())
            else:
                # This is the first data line, put it back by resetting file position
                break
    
    print(f"Found {len(header_comments)} header comment lines")
    
    # Read first chunk to validate columns
    try:
        first_chunk = pd.read_csv(input_path, sep='\t', comment='#', nrows=1000)
    except Exception as e:
        raise ValueError(f"Error reading file as TSV: {e}")
    
    # Check if isE2GLink column already exists
    if 'isE2GLink' in first_chunk.columns:
        print("isE2GLink column already exists")
        return str(input_path)
    
    # Validate required columns exist
    required_cols = [score_col, fdr_col, rho_col]
    missing_cols = [col for col in required_cols if col not in first_chunk.columns]
    
    if missing_cols:
        available_cols = list(first_chunk.columns)
        raise ValueError(f"Missing required columns: {missing_cols}. "
                        f"Available columns: {available_cols}")
    
    # Validate that columns contain numeric values
    for col in required_cols:
        if not pd.api.types.is_numeric_dtype(first_chunk[col]):
            # Try to convert to numeric, which will raise error if not possible
            try:
                pd.to_numeric(first_chunk[col], errors='coerce')
            except Exception as e:
                raise ValueError(f"Column '{col}' does not contain numeric values: {e}")
    
    print(f"Validated columns: {required_cols}")
    print(f"Using thresholds - Score: {score_threshold}, FDR: {fdr_threshold}, Rho: {rho_threshold}")
    
    # Process file in chunks
    print(f"Processing file in chunks of {chunk_size} rows...")
    
    # Initialize output file with header comments
    with gzip.open(output_path, 'wt') as output_file:
        # Write header comments
        for comment in header_comments:
            output_file.write(comment + '\n')
        
        # Process chunks
        first_chunk_written = False
        
        for chunk_num, chunk in enumerate(pd.read_csv(input_path, sep='\t', comment='#', 
                                                     chunksize=chunk_size)):
            
            # Convert columns to numeric, handling any conversion issues
            for col in required_cols:
                chunk[col] = pd.to_numeric(chunk[col], errors='coerce')
            
            # Create boolean column based on thresholds
            # isE2GLink is True if ALL conditions are met:
            # - Score >= score_threshold
            # - Score_fdr <= fdr_threshold  
            # - Score_rho >= rho_threshold
            # If any value is NaN, that condition is considered False
            
            chunk['isE2GLink'] = (
                (chunk[score_col] >= score_threshold) & 
                (chunk[fdr_col] <= fdr_threshold) & 
                (chunk[rho_col] >= rho_threshold)
            )
            
            # Write chunk to output file
            chunk.to_csv(output_file, sep='\t', index=False, 
                        header=(not first_chunk_written))
            
            first_chunk_written = True
            
            if (chunk_num + 1) % 10 == 0:
                print(f"Processed {(chunk_num + 1) * chunk_size} rows...")
    
    print(f"Successfully created: {output_path}")
    
    # Verify the output file
    try:
        verify_chunk = pd.read_csv(output_path, sep='\t', comment='#', nrows=100)
        if 'isE2GLink' not in verify_chunk.columns:
            raise ValueError("Output file verification failed: isE2GLink column not found")
        
        true_count = verify_chunk['isE2GLink'].sum()
        print(f"Verification: Found {true_count} True values in first 100 rows")
        
    except Exception as e:
        print(f"Warning: Could not verify output file: {e}")
    
    return str(output_path)


# Example usage function with error handling
def add_boolean_score_to_SCARlink(file_path, **kwargs):
    """
    Wrapper function with comprehensive error handling.
    
    Parameters:
    -----------
    file_path : str
        Path to the e2g file
    **kwargs : dict
        Additional arguments to pass to process_e2g_file
    
    Returns:
    --------
    str : Path to output file
    """
    try:
        return process_e2g_file(file_path, **kwargs)
    
    except FileNotFoundError as e:
        print(f"File not found error: {e}")
        raise
    
    except ValueError as e:
        print(f"Data validation error: {e}")
        raise
    
    except pd.errors.EmptyDataError:
        print("Error: The file appears to be empty or corrupted")
        raise
    
    except Exception as e:
        print(f"Unexpected error processing file: {e}")
        raise

def main(input_file):
    print(f"Adding boolean isE2GLink to SCARlink predictions: {input_file}")
    add_boolean_score_to_SCARlink(input_file, output_file)
    print(f"Find updated file at: {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage adjust_SCARlink.py <input_file_path> <output_file_path>")
        sys.exit(1)
 
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    main(input_file, output_file)