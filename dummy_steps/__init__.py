import importlib
import logging

logger = logging.getLogger(__name__)

def get_step_function(step_number: int):
    """
    Dynamically import and return the run function for the given dummy step.
    This is for testing without using real Metashape steps.
    """
    module_name = f"steps_dummy.step_{step_number:02d}"
    try:
        module = importlib.import_module(module_name)
        return module.run
    except ImportError:
        logger.error(f"Step {step_number} not found in steps_dummy/")
        return None
    except AttributeError:
        logger.error(f"Step {step_number} does not have a 'run' function.")
        return None
