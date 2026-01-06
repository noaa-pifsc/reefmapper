import importlib
import logging
import glob
import os
from pathlib import Path

logger = logging.getLogger(__name__)

def get_step_function(step_number: int):
    """
    Dynamically import and return the 'run' function for the given step number.
    Assumes step files follow the pattern: step_XX_description.py
    """
    # Always resolve relative to this package directory
    steps_dir = Path(__file__).parent
    pattern = str(steps_dir / f"step_{step_number:02d}_*.py")
    matches = glob.glob(pattern)

    if not matches:
        raise ImportError(f"No module found for step {step_number}")

    module_filename = os.path.basename(matches[0])
    module_name = f"{__package__}.{module_filename[:-3]}"  # e.g. steps.step_01_initialize

    try:
        module = importlib.import_module(module_name)
        return getattr(module, "run")
    except (ModuleNotFoundError, AttributeError) as e:
        raise ImportError(f"Could not load step {step_number}: {e}")
