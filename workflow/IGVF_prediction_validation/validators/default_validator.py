"""
default validator confirms that the score is numeric
no rescue function
"""

# In validators/default_validator.py

import logging
from pathlib import Path
import pandas as pd
from .base import BaseValidator

logger = logging.getLogger(__name__)

class DefaultValidator(BaseValidator):
    def __init__(self, score_col_name: str = "Score"):
        """
        A default validator for models with no specific rules.
        
        Args:
            score_col_name (str): The name of the score column to validate.
        """
        self.score_col_name = score_col_name
        logger.debug(f"DefaultValidator initialized for score column '{self.score_col_name}'")

    def is_score_valid(self, file_path: Path, check_all_rows: bool = False) -> bool:
        """
        Checks that the score column exists and contains only numeric data.
        """
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
                chunksize=100000
            )

            # We just need to consume the iterator to trigger the validation.
            # The work is done by pandas' parser.
            for chunk in chunk_iterator:
                # If we're not checking all rows, the first valid chunk is enough.
                if not check_all_rows:
                    break
            
            return True

        except (ValueError, KeyError) as e:
            # This is the expected failure for a non-compliant file.
            logger.warning(f"[Default Check] Score validation failed for {file_path.name}. Reason: {e}")
            return False
        except FileNotFoundError:
            logger.error(f"[Default Check] File not found during score validation: {file_path}")
            return False
        except Exception as e:
            logger.error(f"[Default Check] An unexpected error occurred reading {file_path}: {e}", exc_info=True)
            return False

    def rescue(self, file_path: Path) -> Path | None:
        """
        The DefaultValidator does not support a rescue operation.
        This method will always fail.
        """
        logger.warning(
            f"No rescue operation available for {file_path.name} "
            f"under the DefaultValidator. File will be marked as failed."
        )
        return None