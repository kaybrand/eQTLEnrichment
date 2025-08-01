import abc
from pathlib import Path

class BaseValidator(abc.ABC):
    """Abstract Base Class for all model-specific validators."""
    
    @abc.abstractmethod
    def is_score_valid(self, file_path):
        """Check if the score column is valid for this specific model."""
        pass

    @abc.abstractmethod
    def rescue(self, file_path):
        """Attempt to fix the file and return the path to the new, fixed file."""
        pass

    def pre_validate_and_reformat(self, file_path: Path) -> tuple[Path | None, bool]:
        """
        Optional hook for models that need reformatting BEFORE standard validation.
        
        This method is called once at the very beginning of the validation process.
        
        - If a file is already in the correct format, this should do nothing and
          return the original path and a status indicating no change was made.
        - If it needs reformatting, it should run its script, create a new file,
          and return the path to that new file and a status indicating it was reformatted.
        - If reformatting fails, it should return (None, False).

        Returns:
            A tuple of (Path | None, bool):
            - The path to the file to be used for subsequent validation.
            - A boolean indicating if the file was modified (True if changed, False if not).
        """
        # Default implementation: do nothing, the file was not changed.
        return (file_path, False)

