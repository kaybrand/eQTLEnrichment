# In validators/scE2g_validator.py

import logging
from pathlib import Path
import subprocess
from .base import BaseValidator
from validation_orchestrator import _validate_core_format # Import the core validator helper

logger = logging.getLogger(__name__)

class scE2GValidator(BaseValidator):
    def __init__(self, score_col_name: str = "Score"):
        self.score_col_name = score_col_name

    def pre_validate_and_reformat(self, file_path: Path) -> tuple[Path | None, bool]:
        """
        For scE2G, we check core format first. If it fails, we assume it's in the
        original format and attempt to reformat it.

        Prootype version uses placeholders for file header metadata (IGVF format)
        THIS WILL NEED TO FIND ALL INFORMATION NEEDED IF SERIOUS REFORMAT EXPECTED
        """
        logger.info(f"Running scE2G pre-validation hook for {file_path.name}...")
        
        # 1. First, check if the file is ALREADY valid.
        if _validate_core_format(file_path):
            logger.info(f"{file_path.name} is already in standard format. No reformatting needed.")
            return (file_path, False) # Return original path, indicate no change

        # 2. If it's not valid, it likely needs reformatting.
        logger.warning(f"{file_path.name} failed initial core check. Attempting reformatting...")
        output_path = file_path.with_name(file_path.name.replace('.e2g.tsv.gz', '_reformatted.e2g.tsv.gz'))
        
        # This is your script that0 fixes chr->ElementChr, start->ElementStart, etc.
        reformat_script = "workflow/IGVF_prediction_validation/adjust_format_scE2GFamily.R" 
        command = [
            "Rscript", reformat_script, 
            "-i", str(file_path), 
            "-o", str(output_path),
            "-c", "placeholder", # Use cell cluster type from file path or CellType; use as key to fill in metadata
            "-d", "placeholder",
            "-s", "placeholder",
            "-m", "scE2G", # take model type from file path
            "-v", "1.0"
            ]

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
        # This is essentially the same as your DefaultValidator's check.
        try:
            pd.read_csv(file_path, sep='\t', compression='gzip', usecols=[self.score_col_name], dtype='float', nrows=10)
            return True
        except (ValueError, KeyError):
            return False

    def rescue(self, file_path: Path) -> Path | None:
        # For scE2G, there is no separate rescue for the score. The reformatting
        # is all-or-nothing. So this does nothing.
        logger.warning(f"scE2G model has no specific score rescue. If score is invalid after reformatting, the file will fail.")
        return None