import logging
import subprocess
from pathlib import Path  # Use pathlib for robust path handling
import pandas as pd
from .base import BaseValidator

logger = logging.getLogger(__name__)

# Define paths at the module level for reliability
THIS_SCRIPT_DIR = Path(__file__).resolve().parent
RESCUE_SCRIPTS_DIR = THIS_SCRIPT_DIR.parent / 'rescue_scripts'

class AbsCorrelationValidator(BaseValidator):
    def __init__(self, score_col_name = "abs_Score", original_score_col = "Score"):
        """
        Args:
            score_col_name (str): The name of the TARGET column (e.g., 'abs_Score').
            original_score_col (str | None): The name of the SOURCE column to read from
                                             during rescue (e.g., 'Score'). If None,
                                             defaults to 'Score'.
        """
        self.score_col_name = score_col_name
        self.original_score_col = original_score_col # The column we read from to create abs_Score

    def is_score_valid(self, file_path: Path, check_all_rows: bool = False) -> bool:
        """
        Checks if the 'abs_Score' column exists and is valid (numeric, 0 <= x <= 1).
        This will typically fail on the first pass, triggering the rescue.
        """
        try:
            # Use the 'boolean' dtype for robust checking. Pandas will raise an error
            # if the column is missing (from usecols) or contains non-boolean values.
            chunk_iterator = pd.read_csv(
                file_path,
                sep='\t',
                comment='#',
                compression='gzip',
                usecols=[self.score_col_name],  # Will fail if abs_Score doesn't exist
                dtype={self.score_col_name: 'float'},
                chunksize=100000
            )

            for chunk in chunk_iterator:
                # Confirm that all score values are plausible abs(correlations)
                if not chunk[self.score_col_name].between(0, 1).all():
                    logger.warning(f"Validation failed for {file_path.name}: "
                                   f"Score contains negative values or values > 1.")
                    return False
                if not check_all_rows:
                    break
            return True

        except (ValueError, KeyError) as e:
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
        """Calls an external script to take the absolute value of the 'Score' column."""
        if file_path.name.endswith('.e2g.tsv.gz'):
            output_filename = file_path.name.replace('.e2g.tsv.gz', '_absolute.e2g.tsv.gz')
            output_path = file_path.with_name(output_filename)
        else:
            logger.error(f"Cannot rescue: Input file must end with '.e2g.tsv.gz', got: {file_path.name}")
            return None

        # Construct the full, absolute path to the rescue script.
        rescue_script_path = RESCUE_SCRIPTS_DIR / 'adjust_corrFamily.py'
        if not rescue_script_path.is_file():
            logger.error(f"Rescue script not found at the expected path: {rescue_script_path}")
            return None

        command = [
            "python", rescue_script_path,
            str(file_path),
            str(output_path),
            "--source-score-col", self.original_score_col,
            "--target-score-col", self.score_col_name
        ]
        
        try:
            result = subprocess.run(command, check=True, capture_output=True, text=True)
            logger.info(f"Successfully added '{self.score_col_name}' column for {file_path.name}.")
            logger.debug(f"Rescue script stdout: {result.stdout}")
            return output_path
        except FileNotFoundError:
            logger.error(f"Rescue script not found at '{rescue_script_path}'. Cannot rescue {file_path.name}.")
            return None
        except subprocess.CalledProcessError as e:
            logger.error(f"abs(correlation) rescue script failed for {file_path.name}. Stderr: {e.stderr}")
            return None