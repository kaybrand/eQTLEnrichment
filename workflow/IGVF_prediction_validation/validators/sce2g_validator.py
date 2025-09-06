# In validators/scE2g_validator.py

"""
Note: Aug 26, 2025
If you integrate the metadata for a proper IGVF format document, you will also be able to store the reformatted file in the correct $SCRATCH folder, saving space in $OAK
"""

import logging
from pathlib import Path
import subprocess
import pandas as pd
import os
from .base import BaseValidator
from validation_orchestrator import _validate_core_format # Import the core validator helper

logger = logging.getLogger(__name__)

THIS_SCRIPT_DIR = Path(__file__).resolve().parent
RESCUE_SCRIPTS_DIR = THIS_SCRIPT_DIR.parent / 'rescue_scripts'

class scE2GValidator(BaseValidator):
    def __init__(self, score_col_name = "Score", original_score_col = "Score"):
        """
        Args:
            score_col_name (str): The name of the TARGET column
            original_score_col (str | None): The name of the SOURCE column to read from
                                             during rescue
        """
        self.score_col_name = score_col_name
        self.original_score_col = original_score_col

    def pre_validate_and_reformat(self, file_path: Path) -> tuple[Path | None, bool]:
        """
        For scE2G, we check core format first. If it fails, we assume it's in the
        original format and attempt to reformat it.

        Prootype version uses placeholders for file header metadata (IGVF format)
        Currently makes one file intended to be used for all members of scE2G family (excl. scATAC)
        THIS WILL NEED TO FIND ALL INFORMATION NEEDED IF SERIOUS REFORMAT EXPECTED
        """
        logger.info(f"Running scE2G pre-validation hook for {file_path.name}...")
        
        # 1. First, check if the file is ALREADY valid.
        if _validate_core_format(file_path):
            logger.info(f"{file_path.name} is already in standard format. No reformatting needed.")
            return (file_path, False) # Return original path, indicate no change

        # 2. If it's not valid, it likely needs reformatting.
        logger.warning(f"{file_path.name} failed initial core check. Attempting reformatting...")

        # Generate output path
        if file_path.name.endswith('.e2g.tsv.gz'):
            output_filename = file_path.name.replace('.e2g.tsv.gz', '_reformated.e2g.tsv.gz')
            output_path = file_path.with_name(output_filename)
        elif file_path.name == 'encode_e2g_predictions.tsv.gz':
            output_path = file_path.with_name('scE2G_predictions_for_eQTL_reformated.e2g.tsv.gz')
        else:
            logger.error(f"Cannot rescue: Input file must end with '.e2g.tsv.gz', got: {file_path.name}")
            return None

        # check if a reformated version aalready exists (should be triggered often because these files are reused)
        if os.path.exists(output_path) and _validate_core_format(output_path):
            logger.info(f"{output_path.name} has already been created.  Using this.")
            return (output_path, True) # Return output path
        
        # This is the script that fixes chr->ElementChr, start->ElementStart, etc.
        rescue_script_path = RESCUE_SCRIPTS_DIR / 'adjust_format_scE2GFamily.R'
        if not rescue_script_path.is_file():
            logger.error(f"Rescue script not found at the expected path: {rescue_script_path}")
            return None

        command = [
            "Rscript", str(rescue_script_path), 
            "-i", str(file_path), 
            "-o", str(output_path),
            # While part of the official IGVF format, these fields are not needed for the eQTL pipeline:
            # "-c", "placeholder", # Use cell cluster type from file path or CellType; use as key to fill in metadata
            # "-d", "placeholder",
            # "-s", "placeholder",
            # "-m", "scE2G", # take model type from file path
            # "-v", "1.0",
            "-a"
            ]
        # print(f'Will run {" ".join(command)}')

        try:
            subprocess.run(command, check=True, capture_output=True, text=True)
            logger.info(f"Successfully reformatted {file_path.name} to {output_path.name}.")
            return (output_path, True) # Return NEW path, indicate it was changed
        except Exception as e:
            logger.error(f"scE2G reformatting script failed for {file_path.name}. Stderr: {e.stderr if hasattr(e, 'stderr') else e}")
            return (None, False) # Signal catastrophic failure

    def is_score_valid(self, file_path: Path, check_all_rows: bool = False) -> bool:
        # This method is now simple. It assumes the core format is correct and just
        # checks the score column.
        try:
            # Using dtype='float' is the most efficient way to validate numeric content.
            # Pandas will raise a ValueError if the column is missing or contains non-numeric text.
            chunk_iterator = pd.read_csv(
                file_path,
                sep='\t',
                comment='#',
                compression='gzip',
                usecols=[self.score_col_name],
                dtype={self.score_col_name: 'float'},
                chunksize=10000
            )

            # We just need to consume the iterator to trigger the validation.
            # The work is done by pandas' parser.
            for _ in chunk_iterator:
                # If we're not checking all rows, the first valid chunk is enough.
                if not check_all_rows:
                    break
            
            return True

        except (ValueError, KeyError) as e:
            # This is the expected failure for a non-compliant file.
            logger.warning(f"[scE2G Family Check] Score validation failed for {file_path.name}. Reason: {e}")
            return False
        except FileNotFoundError:
            logger.error(f"[scE2G Family Check] File not found during score validation: {file_path}")
            return False
        except Exception as e:
            logger.error(f"[scE2G Family Check] An unexpected error occurred reading {file_path}: {e}", exc_info=True)
            return False

    def rescue(self, file_path: Path) -> Path | None:
        # For scE2G, there is no separate rescue for the score. The reformatting
        # is all-or-nothing. So this does nothing.
        logger.warning(f"scE2G model has no specific score rescue. If score is invalid after reformatting, the file will fail.")
        return None