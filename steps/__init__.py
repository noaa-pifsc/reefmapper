import importlib
import logging
import os
import glob

logger = logging.getLogger(__name__)

def get_step_function(step_number: int):
    """Dynamically import and return the run function for the given step."""
    pattern = f"steps/step_{step_number:02d}_*.py"
    matches = glob.glob(pattern)

    if not matches:
        raise ImportError(f"No module found for step {step_number}")

    module_filename = os.path.basename(matches[0])
    module_name = f"steps.{module_filename[:-3]}"  # Remove .py extension

    try:
        module = importlib.import_module(module_name)
        return getattr(module, "run")
    except (ModuleNotFoundError, AttributeError) as e:
        raise ImportError(f"Could not load step {step_number}: {e}")


def run_metashape_pipeline(job, db, Metashape, html2text, shutil, os, math):
    folder_name = job.site_id
    try:
        doc = Metashape.app.document
        chunk = None
        opf = None

        for step_number in range(job.start_step, job.end_step + 1):
            logger.info(f"Starting Step {step_number}")

            run_func = get_step_function(step_number)

            if step_number == 1:
                doc, chunk, opf = run_func(job, doc, db, logger, os, shutil, Metashape)
            elif step_number == 2:
                chunk, opf = run_func(job, doc, chunk, db, logger, opf, quality, root_path, folder_name)
            else:
                chunk, opf = run_func(job, doc, chunk, db, logger, opf)
            logger.info(f"Completed Step {step_number} successfully")

        db.update_job_status(job.id, "completed")
        logger.info(f"Job {job.id} ({folder_name}) processing completed successfully.")

    except Exception as e:
        error_msg = f"Job {job.id} ({folder_name}) failed at step {step_number}: {str(e)}"
        logger.error(error_msg)
        db.log_step(job.id, step_number, "error", error_msg)
        db.update_job_status(job.id, "failed", error_msg)
        raise
