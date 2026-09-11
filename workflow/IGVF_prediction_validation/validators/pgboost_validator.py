# In validators/pgboost_validator.py

import logging
import subprocess
import sys
import pandas as pd
from pathlib import Path
from .base import BaseValidator # This will work when you run as a module

# Get a logger for this module
logger = logging.getLogger(__name__)

# Define paths at the module level for reliability
THIS_SCRIPT_DIR = Path(__file__).resolve().parent
RESCUE_SCRIPTS_DIR = THIS_SCRIPT_DIR.parent / 'rescue_scripts'

# pgBoost saturates at Score == 1 for a large block of links (18-19% of rows in the
# 2026 submission). The uploader asked for those to be excluded from the benchmark.
EXCLUDED_PGBOOST_SCORE = 1.0

class pgBoostValidator(BaseValidator):
    """Validates pgBoost predictions: 0-1 Percentile, with Score == 1 links excluded.

    pgBoost diverges from GEEES here. Both rank-transform a score into a 0-1
    Percentile, but only pgBoost drops its saturated Score == 1 links first. See
    GEEESValidator for the un-filtered variant; the two are kept as separate classes
    so a change to pgBoost's handling cannot silently alter GEEES.
    """

    def __init__(self, score_col_name = "Percentile", original_score_col = "Score"):
        """
        Args:
            score_col_name (str): The name of the TARGET column
            original_score_col (str | None): The name of the SOURCE column to read from
                                             during rescue
        """
        self.score_col_name = score_col_name
        self.original_score_col = original_score_col

    def is_score_valid(self, file_path, check_all_rows=False):
        """
        Checks that Percentile exists, is numeric and lies in [0, 1], AND that no
        link retains the saturated Score == 1.
        """
        try:
            # Read the percentile column and the raw score together: a submitted file
            # can have a perfectly valid Percentile while still containing the
            # Score == 1 links that must be removed.
            usecols = [self.score_col_name]
            if self.original_score_col and self.original_score_col != self.score_col_name:
                usecols.append(self.original_score_col)

            chunk_iterator = pd.read_csv(
                file_path,
                sep='\t',
                comment='#',
                compression='gzip',
                usecols=usecols, # This will raise ValueError if a column is missing
                chunksize=100000,
                low_memory=False
            )

            for i, chunk in enumerate(chunk_iterator):
                # Ensure the column is numeric. `to_numeric` will raise an error on non-numeric values.
                score_col = pd.to_numeric(chunk[self.score_col_name])

                # Vectorized check is faster than `all()` on a Python loop.
                # Range is 0-1, matching the submitted pgBoost files and the 0.95
                # threshold in config_methods.tsv. Accepting 0-100 here (as this
                # previously did) would pass a file on the wrong scale, and a 0.95
                # threshold against 0-100 values admits nearly every pair.
                if not ((score_col >= 0) & (score_col <= 1)).all():
                    logger.warning(f"Validation failed for {file_path}: Percentiles found outside [0, 1] range.")
                    return False

                # Saturated links must already have been removed. Because they are all
                # tied at the ceiling, leaving them in also distorts every other
                # percentile: Percentile >= 0.95 admits ~22% of pairs rather than 5%.
                if self.original_score_col in chunk.columns:
                    raw = pd.to_numeric(chunk[self.original_score_col])
                    if (raw == EXCLUDED_PGBOOST_SCORE).any():
                        logger.info(f"{file_path.name} still contains "
                                    f"{self.original_score_col} == {EXCLUDED_PGBOOST_SCORE} links. "
                                    f"Proceeding to rescue to exclude them.")
                        return False

                # If we're only checking the first part of the file, break after one chunk.
                if not check_all_rows:
                    break

            # If all chunks passed (or the first one passed on a quick check), the score is valid.
            return True

        # This will catch both the missing column and non-numeric score values.
        except (ValueError, KeyError) as e:
            # This is the EXPECTED failure path for a file that needs rescuing.
            logger.info(f"'{self.score_col_name}' not found or invalid in {file_path.name}. "
                        f"Proceeding to rescue.")
            return False
        except FileNotFoundError:
            logger.error(f"File not found during score validation: {file_path}")
            return False
        # A broader exception for unexpected pandas errors
        except Exception as e:
            logger.error(f"An unexpected error occurred reading {file_path}: {e}", exc_info=True)
            return False

    def rescue(self, file_path):
        """Calls an external script to add the 'Percentile' column."""
        if file_path.name.endswith('.e2g.tsv.gz'):
            output_filename = file_path.name.replace('.e2g.tsv.gz', '_percentile.e2g.tsv.gz')
            output_path = file_path.with_name(output_filename)
        else:
            logger.error(f"Cannot rescue: Input file must end with '.e2g.tsv.gz', got: {file_path.name}")
            return None
        
        rescue_script_path = RESCUE_SCRIPTS_DIR / 'adjust_pgBoost.py'
        if not rescue_script_path.is_file():
            logger.error(f"Rescue script not found at the expected path: {rescue_script_path}")
            return None

        command = [
            # Use the interpreter running this validator, not a bare "python": on this
            # cluster /usr/bin/python is Python 2.7, which cannot parse the rescue
            # scripts' f-strings, so every rescue died with a SyntaxError.
            sys.executable, rescue_script_path,
            "--input", str(file_path),
            "--output", str(output_path),
            "--source-score-col", self.original_score_col,
            "--percentile-score-col", self.score_col_name,
            # Drop the saturated links, then rank what remains. Percentiles are
            # recomputed on the retained rows so the ranking describes the set
            # actually being benchmarked. GEEESValidator omits this flag.
            "--exclude-score-equal", str(EXCLUDED_PGBOOST_SCORE)
        ]

        try:
            # Use the logger instance we defined at the top of the file
            subprocess.run(command, check=True, capture_output=True, text=True)
            logger.info(f"Successfully ran pgBoost rescue script. New file at {output_path}")
            return output_path
        except FileNotFoundError:
            logger.error(f"Rescue script not found at {command[1]}. Cannot rescue {file_path}.")
            return None
        except subprocess.CalledProcessError as e:
            logger.error(f"pgBoost rescue script failed for {file_path}. Stderr: {e.stderr}")
            return None