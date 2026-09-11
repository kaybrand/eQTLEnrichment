# In validators/geees_validator.py

import logging
import subprocess
import sys
import pandas as pd
from pathlib import Path
from .base import BaseValidator

# Get a logger for this module
logger = logging.getLogger(__name__)

# Define paths at the module level for reliability
THIS_SCRIPT_DIR = Path(__file__).resolve().parent
RESCUE_SCRIPTS_DIR = THIS_SCRIPT_DIR.parent / 'rescue_scripts'


class GEEESValidator(BaseValidator):
    """Rank-transforms a raw score into a 0-1 Percentile column.

    GEEES publishes a raw positive Score with no Percentile column and no
    established score cutoff, so it is thresholded on its own rank instead
    (Percentile >= 0.95 keeps the top 5% of each biosample's predictions).

    This is the behaviour pgBoostValidator had before pgBoost diverged: pgBoost
    now additionally excludes its Score == 1 links, which GEEES must not do.
    Keep the two separate rather than parameterizing one class, so a change to
    pgBoost's handling cannot silently alter GEEES.
    """

    def __init__(self, score_col_name="Percentile", original_score_col="Score"):
        """
        Args:
            score_col_name (str): The name of the TARGET column
            original_score_col (str | None): The name of the SOURCE column to read from
                                             during rescue
        """
        self.score_col_name = score_col_name
        self.original_score_col = original_score_col

    def is_score_valid(self, file_path, check_all_rows=False):
        """Checks that the percentile column exists, is numeric, and lies in [0, 1]."""
        try:
            chunk_iterator = pd.read_csv(
                file_path,
                sep='\t',
                comment='#',
                compression='gzip',
                usecols=[self.score_col_name],  # raises ValueError if column is missing
                chunksize=100000,
                low_memory=False
            )

            for chunk in chunk_iterator:
                score_col = pd.to_numeric(chunk[self.score_col_name])

                # 0-1, matching the submitted pgBoost files and the 0.95 threshold in
                # config_methods.tsv. A 0-100 file thresholded at 0.95 would admit
                # nearly every pair, so reject it and let the rescue recompute.
                if not ((score_col >= 0) & (score_col <= 1)).all():
                    logger.warning(f"Validation failed for {file_path}: "
                                   f"Percentiles found outside [0, 1] range.")
                    return False

                if not check_all_rows:
                    break

            return True

        except (ValueError, KeyError):
            # EXPECTED failure path for a file that needs rescuing: GEEES files ship
            # with no Percentile column at all.
            logger.info(f"'{self.score_col_name}' not found or invalid in {file_path.name}. "
                        f"Proceeding to rescue.")
            return False
        except FileNotFoundError:
            logger.error(f"File not found during score validation: {file_path}")
            return False
        except Exception as e:
            logger.error(f"An unexpected error occurred reading {file_path}: {e}", exc_info=True)
            return False

    def rescue(self, file_path):
        """Calls the shared rank-percentile script to add the 'Percentile' column."""
        if file_path.name.endswith('.e2g.tsv.gz'):
            output_filename = file_path.name.replace('.e2g.tsv.gz', '_percentile.e2g.tsv.gz')
            output_path = file_path.with_name(output_filename)
        else:
            logger.error(f"Cannot rescue: Input file must end with '.e2g.tsv.gz', "
                         f"got: {file_path.name}")
            return None

        rescue_script_path = RESCUE_SCRIPTS_DIR / 'adjust_pgBoost.py'
        if not rescue_script_path.is_file():
            logger.error(f"Rescue script not found at the expected path: {rescue_script_path}")
            return None

        # No --exclude-score-equal here: every GEEES link is ranked. Bare "python"
        # would resolve to the system Python 2.7 on this cluster and fail to parse
        # the rescue script's f-strings.
        command = [
            sys.executable, str(rescue_script_path),
            "--input", str(file_path),
            "--output", str(output_path),
            "--source-score-col", self.original_score_col,
            "--percentile-score-col", self.score_col_name
        ]

        try:
            subprocess.run(command, check=True, capture_output=True, text=True)
            logger.info(f"Successfully ran GEEES percentile rescue. New file at {output_path}")
            return output_path
        except FileNotFoundError:
            logger.error(f"Rescue script not found at {rescue_script_path}. Cannot rescue {file_path}.")
            return None
        except subprocess.CalledProcessError as e:
            logger.error(f"GEEES percentile rescue script failed for {file_path}. Stderr: {e.stderr}")
            return None
