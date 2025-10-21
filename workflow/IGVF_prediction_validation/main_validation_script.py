"""
main_script.py
Kayla Brand

Checks all file paths to confirm these are valid prediction files that will be accepted by the process predictions step in the eQTL pipeline.
Can make specific adjustments for known models, for example, calculating percentile of scores for pgBoost predictions.

To add new model fixes, define a new validator inheriting from the base.py abstract class, import it to this document, and 
indicate whch models use this validator in the validator map.  Unrecognized models will be processed with the default validator, whose
only requirement for valid scores is that they be interpretable as a float.

Requires -p, --prediction-config path/to/eQTL/prediction/paths/config/table
Requires -m, --method-config path/to/
"""

import logging
import pandas as pd
import numpy as np
import os
import click
from pathlib import Path

from validation_orchestrator import process_prediction_file
# Using the __init__.py trick
from validators import (
    pgBoostValidator,
    SCARlinkValidator,
    AbsCorrelationValidator,
    DefaultValidator,
    scE2GValidator
)

VALIDATOR_MAP = {
    "pgBoost": pgBoostValidator,
    "SCARlink": SCARlinkValidator,
    # "Signac": AbsCorrelationValidator,
    # "Cicero": AbsCorrelationValidator,
    # "ArchR": AbsCorrelationValidator,
    "SCENT": DefaultValidator, # p-values
    "scE2G_multiome": scE2GValidator,
    "scE2G_ATAC": scE2GValidator,
    "scABC": scE2GValidator,
    "Kendall": scE2GValidator,
    "ABC_distanceToTSS": scE2GValidator,
} # FigR, EPCOT will also use default validators

DEFAULT_NA_VALUES = pd.io.parsers.readers.STR_NA_VALUES
# This set includes: {'', 'NA', 'NaN', 'NULL', 'n/a', 'nan', 'null', ...and many more}

def is_intentional_blank(value) -> bool:
    """
    Checks if a value is NaN, None, or a common string representation of a blank/NA.

    This is used to differentiate user-intended blanks from actual file path typos.
    - Returns True for: np.nan, None, "", " ", "NA", "nan", "NULL"
    - Returns False for: "data/my_file.txt", "some_other_string", "0"
    """
    # 1. The most robust check for any NaN-like object (np.nan, None, pd.NaT).
    # This will catch values that pd.read_csv converted automatically.
    if pd.isna(value):
        return True
    
    # 2. Check if it's a string that matches one of Pandas' default NA values.
    # We strip whitespace to handle cells containing only spaces (" ").
    if isinstance(value, str) and value.strip() in DEFAULT_NA_VALUES:
        return True
    
    return False

def process_table(df_table, config, validator_map, check_all_rows=False, verbose=False):
    """
    Iterates through an eQTL Benchmark prediction configuration table and 
    verifies each prediction file for pipeline compatibility.

    Args:
        df_table (pd.DataFrame): DataFrame of prediction file paths (rows=biosamples, cols=models).
        config (pd.DataFrame): DataFrame with configuration for each model.
        check_all_rows (bool): Flag to validate all rows or just the head.

    Returns:
        pd.DataFrame: A new DataFrame with paths to validated/rescued files or NaN.
    """
    validated_data = {}  # Store validated file paths. keys are biosamples. cols are method names

    for biosample in df_table.index: # iterate over the table's index. The rows indicate the number of samples to process in the eQTL pipeline
        validated_data[biosample] = {}
        row = df_table.loc[biosample]

        for model_name, cell_value in row.items():
            # Verify the model name and cell contents are valid

            # Don't validate GTEx tissue mappings; Handle special, non-path metadata columns explicitly
            if model_name == 'GTExTissue': 
                if is_intentional_blank(cell_value):
                    validated_data[biosample][model_name] = "" # Store as clean empty string
                else:
                    validated_data[biosample][model_name] = str(cell_value).strip() # Keep original value
                continue

            # Standardize empty/NaN cells
            if is_intentional_blank(cell_value): # includes blanks and a range of NA options
                validated_data[biosample][model_name] = np.nan
                continue

            # This is intended as a file path. Save the cell contents as a file Path object
            file_path = Path(cell_value.strip())
            
            # --- Dynamic Validator Configuration ---
            try:
                model_config = config.loc[model_name]
                score_col_name = model_config['score_col']
                original_score_col = model_config.get('original_score_col', None) 
            except KeyError:
                logging.error(f"Method '{model_name}' not found in configuration file. Skipping {file_path}")
                validated_data[biosample][model_name] = np.nan
                continue

            ValidatorClass = validator_map.get(model_name, DefaultValidator)

            try:
                if ValidatorClass in [AbsCorrelationValidator, pgBoostValidator, SCARlinkValidator]:
                    validator_to_use = ValidatorClass(
                        score_col_name=score_col_name,
                        original_score_col=original_score_col
                    )
                # Add elif for other validators with special configs
                else: # For simple validators like Default
                    validator_to_use = ValidatorClass(score_col_name=score_col_name)

            except TypeError as e:
                logging.error(f"Could not initialize validator for {model_name}. "
                            f"Check config parameters. Error: {e}")
                validated_data[biosample][model_name] = np.nan
                continue

            if validator_to_use is DefaultValidator:
                logging.info(f"No specific validator found for model '{model_name}'. "
                            f"Applying default validation rules to {file_path.name}.")

            # --- Run Validation ---
            result = process_prediction_file(file_path, validator_to_use, check_all_rows, verbose)

            if result.status.startswith('FAILED'):
                validated_data[biosample][model_name] = np.nan
            elif result.status.startswith('PASSED'):
                validated_data[biosample][model_name] = result.final_path
                if verbose and "ORIGINAL" not in result.status:
                    print(f"Replacing {os.path.basename(result.original_path)} with {os.path.basename(result.final_path)}")
    
    # Convert validated data to DataFrame
    validated_df = pd.DataFrame.from_dict(validated_data, orient='index', columns=df_table.columns)
    validated_df.index.name = "biosample"
    return validated_df


@click.command()
@click.option( # or generate this table automatically 
    '-p', '--prediction-config', 'table_file',
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    required=True,
    help='Path to the prediction config TSV containing predictions for each biosample and model.'
)
@click.option(
    '-m', '--method-config', 'config_file',
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    required=True,
    help='Path to the method config TSV containing information on each model (e.g., score column names).'
)
@click.option(
    '-c', '--check-all-rows',
    is_flag=True,  # This makes it a boolean flag, e.g., --check-all-rows
    default=False,
    help='If set, confirm all rows in prediction files are valid scores, else check the head only.'
)
@click.option(
    '-v', '--verbose',
    is_flag=True,  # This makes it a boolean flag
    default=False,
    help='If set, print row by row information for debug'
)
def main(table_file: Path, config_file: Path, check_all_rows: bool, verbose: bool):
    """
    Validates a table of prediction files from various models against core format
    and model-specific score requirements, rescuing files where possible.
    """

    # Configure logging for the entire application, ONCE.
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    )

    logging.info("Starting validation process...")
    logging.info(f"Reading prediction table: {table_file}")
    logging.info(f"Reading method configuration: {config_file}")
    if check_all_rows:
        logging.info("Full file validation enabled (checking all rows).")

    # Read in the model prediction file paths table and the configruation of score column names
    try:
        prediction_table = pd.read_csv(table_file, sep='\t', dtype=str, index_col="biosample")
        method_config = pd.read_csv(config_file, sep='\t', index_col="method", dtype=str)
    except FileNotFoundError as e:
        logging.error(f"File not found: {e}")
        exit(1)

    # 2. DELEGATE the core logic to the helper function
    validated_df = process_table(prediction_table, method_config, VALIDATOR_MAP, check_all_rows, verbose)
    
    # 3. Handle final output
    output_dir = table_file.parent
    output_filename = f"{table_file.stem}_validated{table_file.suffix}"
    validated_table_file = output_dir / output_filename
    
    validated_df.to_csv(validated_table_file, sep='\t')
    logging.info(f"Validation process finished. Validated table saved to: {validated_table_file}")

if __name__ == "__main__":
    main()
