"""
main_script.py
"""

import logging
from validation_orchestrator import process_prediction_file
# Using the __init__.py trick
from validators import (
    pgBoostValidator,
    SCARlinkValidator,
    AbsCorrelationValidator,
    DefaultValidator,
    scE2GValidator
)

# Add click options module Python functionality
def main():
    # Configure logging for the entire application, ONCE.
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        # Optional: Log to a file as well as the console
        # handlers=[
        #     logging.FileHandler("validation.log"),
        #     logging.StreamHandler()
        # ]
    )

    logging.info("Starting validation process...")
    
    VALIDATOR_MAP = {
        "pgBoost": pgBoostValidator(),
        "SCARlink": SCARlinkValidator(),
        "Signac": AbsCorrelationValidator(),
        "Cicero": AbsCorrelationValidator(),
        "ArchR": AbsCorrelationValidator(),
        "SCENT": AbsCorrelationValidator(),
        "scE2G_multiome": scE2GValidator(),
        "scE2G_ATAC": scE2GValidator(),
        "scABC": scE2GValidator(),
        "Kendall": scE2GValidator(),
        "ABC_distanceToTSS": scE2GValidator(),
    } # FigR, EPCOT will use default validators

    default_validator = DefaultValidator()

    # *** ... your main loop that reads the config and calls process_prediction_file ...
    model_name = None
    file_path = None
    # 4. Use dict.get() to retrieve the validator.
    # If model_name is in the map, it returns the specific validator.
    # If not, it returns the `default_validator` instance.
    validator_to_use = VALIDATOR_MAP.get(model_name, default_validator)

    if validator_to_use is default_validator:
        logging.info(f"No specific validator found for model '{model_name}'. "
                     f"Applying default validation rules to {file_path.name}.")

    # The rest of your code remains unchanged!
    result = process_prediction_file(file_path, validator_to_use)

    # Act on the results (drop row or update path)
    if result.status.startswith('FAILED'):
        # ... drop row ...
        pass
    elif result.status == 'PASSED_RESCUED':
        # ... update path ...
        pass


    logging.info("Validation process finished.")

if __name__ == "__main__":
    main()