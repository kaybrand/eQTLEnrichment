import logging
import subprocess
import sys
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
        # Accept a bare *.tsv.gz as well as *.e2g.tsv.gz. Open4Gene ships its
        # standard-format predictions as *_Open4Gene.tsv.gz with no .e2g infix, and
        # requiring the infix made every such file unrescuable.
        if file_path.name.endswith('.e2g.tsv.gz'):
            output_filename = file_path.name.replace('.e2g.tsv.gz', '_absolute.e2g.tsv.gz')
            output_path = file_path.with_name(output_filename)
        elif file_path.name.endswith('.tsv.gz'):
            output_filename = file_path.name.replace('.tsv.gz', '_absolute.e2g.tsv.gz')
            output_path = file_path.with_name(output_filename)
        else:
            logger.error(f"Cannot rescue: Input file must end with '.tsv.gz', got: {file_path.name}")
            return None

        # Construct the full, absolute path to the rescue script.
        rescue_script_path = RESCUE_SCRIPTS_DIR / 'adjust_corrFamily.py'
        if not rescue_script_path.is_file():
            logger.error(f"Rescue script not found at the expected path: {rescue_script_path}")
            return None

        command = [
            # Use the interpreter running this validator, not a bare "python": on this
            # cluster /usr/bin/python is Python 2.7, which cannot parse the rescue
            # scripts' f-strings, so every rescue died with a SyntaxError.
            sys.executable, rescue_script_path,
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

class DivergentAbsValidator(AbsCorrelationValidator):
    """abs() of a signed score that is NOT a correlation, so is not bounded by 1.

    Same rescue path as AbsCorrelationValidator -- adjust_corrFamily.py, which
    already rewrites "# ScoreType: divergent" to "# ScoreType: abs(correlation)"
    and writes an abs column. Only the validity test differs: a correlation's
    absolute value must lie in [0, 1], but Open4Gene's Score is a hurdle-model
    z-statistic spanning roughly -5.7 .. 8.0, so |z| legitimately exceeds 1 and
    the inherited [0, 1] check would reject a perfectly good rescued file.

    Used for the divergent arm of the Open4Gene comparison (Open4Gene_abs), which
    scores |z| against the positive-only arm that thresholds z directly.
    """

    def is_score_valid(self, file_path: Path, check_all_rows: bool = False) -> bool:
        """Valid if the target column exists, is numeric, and is non-negative."""
        try:
            chunk_iterator = pd.read_csv(
                file_path,
                sep='\t',
                comment='#',
                compression='gzip',
                usecols=[self.score_col_name],  # raises if the abs column is absent
                dtype={self.score_col_name: 'float'},
                chunksize=100000
            )
            for chunk in chunk_iterator:
                # No upper bound: this is an abs(z), not an abs(correlation).
                if (chunk[self.score_col_name] < 0).any():
                    logger.warning(f"Validation failed for {file_path.name}: "
                                   f"{self.score_col_name} contains negative values, "
                                   f"which an absolute value cannot.")
                    return False
                if not check_all_rows:
                    break
            return True
        except (ValueError, KeyError):
            # EXPECTED path for a file that has not been abs-transformed yet.
            logger.info(f"'{self.score_col_name}' not found or invalid in {file_path.name}. "
                        f"Proceeding to rescue.")
            return False
        except FileNotFoundError:
            logger.error(f"File not found during score validation: {file_path}")
            return False
        except Exception as e:
            logger.error(f"An unexpected error occurred reading {file_path}: {e}", exc_info=True)
            return False
