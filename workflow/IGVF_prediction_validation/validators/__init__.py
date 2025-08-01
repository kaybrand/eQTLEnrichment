from .base import BaseValidator
from .pgboost_validator import pgBoostValidator
from .default_validator import DefaultValidator
from .correlation_validator import AbsCorrelationValidator
from .scarlink_validator import SCARlinkValidator
from .sce2g_validator import scE2GValidator

# Add a new line here every time you create a new validator file

# Optional: You can also define what `from validators import *` would do
__all__ = [
    'DefaultValidator',
    'pgBoostValidator',
    'AbsCorrelationValidator',
    'SCARlinkValidator',
    'scE2GValidator'
]