"""
Helper functions for the validation script
Use eQTLEnv
Kayla Brand
June 24, 2025
"""
import gzip
import csv
import logging
import re
import numpy as np
import pandas as pd
import os
import sys
import collections
import subprocess

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# A simple object to return structured results
ValidationResult = collections.namedtuple('ValidationResult', ['status', 'final_path'])

def process_prediction_file(file_path, validator_strategy):
    """
    Orchestrates the validation and rescue process for a single file.

    Args:
        file_path (str): The path to the prediction file to validate.
        validator_strategy (BaseValidator): An object that knows how to validate 
                                            and rescue a file for a specific model.

    Returns:
        ValidationResult: An object with a status ('PASSED', 'RESCUED', 'FAILED')
                          and the path to the valid file (or None).
    """
    # if validator_strategy is just reformating, add a reformatting option HERE

    # 1. Fast check: Do the core columns exist and have the right format?
    if not _validate_core_format(file_path):
        logging.error(f"Failed core validation for {file_path}. Skipping score check and rescue.")
        return ValidationResult(status='FAILED_CORE', final_path=None)

    # 2. Model-specific check: Is the score column valid?
    if validator_strategy.is_score_valid(file_path):
        logging.info(f"File {file_path} passed all checks.")
        return ValidationResult(status='PASSED_ORIGINAL', final_path=file_path)

    # 3. Score is not valid. Attempt rescue.
    logging.warning(f"Score validation failed for {file_path}. Attempting rescue...")
    rescued_file_path = validator_strategy.rescue(file_path)

    if not rescued_file_path:
        logging.error(f"Rescue attempt failed for {file_path}. The file will be removed from config.")
        return ValidationResult(status='FAILED_RESCUE', final_path=None)
        
    # 4. A new file was created. We MUST re-validate its core integrity.
    #    This prevents a broken rescue script from poisoning the pipeline.
    logging.info(f"Rescue created new file: {rescued_file_path}. Verifying its core format...")
    if _validate_core_format(rescued_file_path):
        logging.info(f"Rescued file {rescued_file_path} passed core validation.")
        # Does not check score column again
        return ValidationResult(status='PASSED_RESCUED', final_path=rescued_file_path)
    else:
        logging.error(f"Rescued file {rescued_file_path} FAILED core validation. The rescue script is faulty.")
        # Optional: Clean up the bad rescued file
        # os.remove(rescued_file_path)
        return ValidationResult(status='FAILED_POST_RESCUE', final_path=None)

def _validate_core_format(file_path, check_all_rows=False, verbose=False):
    """
    Validates only the core columns (ElementChr, Start, End, GeneSymbol) 
    and their data types, as expected by bedtools. This is the reusable, generic part of the validation.
    Flexible column order.

    Args:
        file_path (str): Path to the gzipped TSV file.
        check_all_rows (bool): If True, checks all rows for validity (slower). Default: False.
        verbose (bool): If True, print notices for valid files as well as invalid files

    Returns:
        bool: True if the file is valid, False otherwise. Prints error messages if invalid.
    """

    try:
        with gzip.open(file_path, 'rt') as f:
            # Read comments and extract header.
            header = None
            for line in f:
                line = line.strip()
                if not line.startswith("#"):
                    header = line.split('\t')
                    break

            if not header:
                logging.error(f"[Core Check] No valid header in {file_path}")
                return False

            required_columns = ["ElementChr", "ElementStart", "ElementEnd", "GeneSymbol"]
            if not all(col in header for col in required_columns):
                logging.error(f"[Core Check] Missing required columns in {file_path}. ")
                logging.error(f"Required: {required_columns}, Found: {header}")
                return False
            
            # Create index lookups once
            col_indices = {col: header.index(col) for col in required_columns}

            # Data validation
            csv_reader = csv.reader(f, delimiter='\t')
            row_count = 0
            for i, row in enumerate(csv_reader):
                row_count += 1
                #Only check limited rows if requested.
                if not check_all_rows and i >= 5:
                    break

                # Basic row integrity check
                if len(row) < len(required_columns+1):  # Check for proper row width
                    logging.error(f"[Core Check] Row {i+1} in {file_path} has wrong number of columns. Found {len(row)}, required at least {len(required_columns)} columns")
                    return False

                # Access the right columns based on their positions in the header:
                element_chr = row[col_indices["ElementChr"]]
                element_start = row[col_indices["ElementStart"]]
                element_end = row[col_indices["ElementEnd"]]
                gene_symbol = row[col_indices["GeneSymbol"]]

                # ElementChr: Check for chromosome format (e.g., chr1, chrX, chrM)
                if not re.match(r"^chr([0-9]+|[XYM])$", element_chr):
                    logging.error(f"Invalid ElementChr value in row {i+1} of {file_path}: {element_chr}")
                    return False

                # ElementStart and ElementEnd: Check for positive integers, disallowing ".0"
                try:
                    if "." in element_start or "." in element_end:
                         raise ValueError("Start or End coordinates contain '.'.  Must be integers.")

                    element_start = int(element_start)
                    element_end = int(element_end)
                    if element_start < 0 or element_end < 0:
                        raise ValueError("Start or end cannot be negative") #Custom error

                except ValueError as e:
                    logging.error(f"Invalid ElementStart/End in row {i+1} of {file_path}: {element_start}, {element_end}. Error: {e}")
                    return False

                # GeneSymbol: Check for alphanumeric (this is a simple check, adjust as needed)
                if not re.match(r"^[A-Za-z0-9\\-]+$", gene_symbol):
                    logging.error(f"Invalid GeneSymbol in row {i+1} of {file_path}: {gene_symbol}")
                    return False

        if verbose:
            logging.info(f"[Core Check] File {file_path} is valid (checked {row_count} rows).")
        return True

    except FileNotFoundError:
        logging.error(f"[Core Check] File not found: {file_path}")
        return False
    except gzip.BadGzipFile:
        logging.error(f"[Core Check] File is not a valid gzip file: {file_path}")
        return False
    except Exception as e:
        logging.error(f"[Core Check] An unexpected error occurred: {e}")
        return False
    


"""
is_scE2G_original_format = len([col for col in ['chr', 'start', 'end', 'TargetGene'] if col not in header]) > 0
"""