# In validators/scarlink_validator.py

import logging
import subprocess
from pathlib import Path  # Use pathlib for robust path handling
import pandas as pd
from .base import BaseValidator

logger = logging.getLogger(__name__)

class SCARlinkValidator(BaseValidator):
    def __init__(self):
        self.score_col_name = "isE2GLink"

    def is_score_valid(self, file_path: Path, check_all_rows: bool = False) -> bool:
        """
        Checks if the isE2GLink score column exists and contains only boolean-like values.
        
        Args:
            file_path (Path): Path to the gzipped TSV file.
            check_all_rows (bool): If True, validates all rows; otherwise, checks the first chunk.

        Returns:
            bool: True if the score column is valid, False otherwise.
        """
        try:
            # Use the 'boolean' dtype for robust checking. Pandas will raise an error
            # if the column is missing (from usecols) or contains non-boolean values.
            chunk_iterator = pd.read_csv(
                file_path,
                sep='\t',
                comment='#',
                compression='gzip',
                usecols=[self.score_col_name],
                dtype={self.score_col_name: 'boolean'}, # This is the key improvement!
                chunksize=100000 if check_all_rows else 1000
            )

            # Iterating through the chunk(s) implicitly triggers the validation.
            # We only need to consume the iterator.
            for chunk in chunk_iterator:
                # If we got here, the chunk is valid. If we're not checking all rows, we're done.
                if not check_all_rows:
                    break
            
            # If we looped through without an error, the file is valid.
            return True

        except (ValueError, KeyError) as e:
            # This will catch:
            # 1. Missing 'isE2GLink' column (from usecols).
            # 2. Values that cannot be converted to boolean (from dtype).
            # This is the EXPECTED failure path for a file that needs rescuing.
            logger.info(f"'{self.score_col_name}' not found or invalid in {file_path.name}. "
                        f"Proceeding to rescue.")
            return False
        except FileNotFoundError:
            logger.error(f"File not found during score validation: {file_path}")
            return False
        except Exception as e:
            logger.error(f"An unexpected error occurred reading {file_path}: {e}", exc_info=True)
            return False

    def rescue(self, file_path: Path) -> Path | None:
        """Calls an external script to add the 'isE2GLink' column."""
        # Use pathlib's `with_name` for safe and correct path construction.
        # This keeps the new file in the same directory as the old one.
        if file_path.name.endswith('.e2g.tsv.gz'):
            output_filename = file_path.name.replace('.e2g.tsv.gz', '_boolean_score.e2g.tsv.gz')
            output_path = file_path.with_name(output_filename)
        else:
            # This check is good, but log it as an error for better tracking.
            logger.error(f"Cannot rescue: Input file must end with '.e2g.tsv.gz', got: {file_path.name}")
            return None

        # Best practice: ensure the path to the script is reliable.
        # This could be made configurable later if needed.
        rescue_script_path = "workflow/IGVF_prediction_validation/adjust_SCARlink.py"

        command = [
            "python",
            rescue_script_path,
            str(file_path),  # Convert Path object to string for subprocess
            str(output_path)
        ]
        
        try:
            result = subprocess.run(command, check=True, capture_output=True, text=True)
            logger.info(f"Successfully ran SCARlink rescue script for {file_path.name}.")
            logger.debug(f"Rescue script stdout: {result.stdout}")
            return output_path
        except FileNotFoundError:
            logger.error(f"Rescue script not found at '{rescue_script_path}'. Cannot rescue {file_path.name}.")
            return None
        except subprocess.CalledProcessError as e:
            logger.error(f"SCARlink rescue script failed for {file_path.name}. Stderr: {e.stderr}")
            return None