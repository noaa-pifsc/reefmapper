import os
import sys
import logging
from datetime import datetime
from database.db_manager_mod import MyDatabaseManager, ProcessingJob
from steps import get_step_function
from steps.utils import parse_marker_pairs
import Metashape

# --- Setup working directory ---
script_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(script_dir)
sys.path.insert(0, script_dir)
print("Working directory:", os.getcwd())
print("Command-line args:", sys.argv)

# --- Initialize logger ---
logger = logging.getLogger("reefmapper")
logger.setLevel(logging.INFO)
ch = logging.StreamHandler()
ch.setFormatter(logging.Formatter("%(asctime)s - %(levelname)s - %(message)s"))
logger.addHandler(ch)

# --- Initialize database manager ---
db = MyDatabaseManager()


MAX_STEP = 7
def MetashapeProcess(job: ProcessingJob):
    valid_markers = parse_marker_pairs(job)
    logger.info(f"Processing job: {job.site_id}")

    # --- Setup paths ---
    products_dir = os.path.join(job.project_path, "Products_automation")
    os.makedirs(products_dir, exist_ok=True)

    log_path = os.path.join(products_dir, f"{job.site_id}_log.txt")
    psx_path = os.path.join(products_dir, f"{job.site_id}.psx")

    job.project_file = psx_path        # derived attribute for .psx file
    job.output_log_path = log_path

    warnings = []

    # Normalize start/end to ensure forward progress
    if job.end_step is None:
        job.end_step = 0
        db.set_end_step(job.id, job.end_step)
    if job.start_step is None or job.start_step <= job.end_step:
        job.start_step = job.end_step + 1
        db.set_start_step(job.id, job.start_step)

    db.update_job_status(job.id, "running")

    # --- Load images ---
    image_files = [
        os.path.join(job.project_path, f)
        for f in os.listdir(job.project_path)
        if f.lower().endswith(".jpg") and not f.startswith('.') and not f.startswith('._')
    ]
    if not image_files:
        raise FileNotFoundError(f"No JPG images found in {job.project_path}")

    # --- Load or create Metashape project ---
    doc = Metashape.Document()
    if os.path.exists(job.project_file):
        logger.info(f"Opening existing project: {job.project_file}")
        doc.open(job.project_file)
        if doc.read_only:
            logger.warning("Project opened in read-only mode. Check for lock files or concurrent access.")
        chunk = doc.chunk
    else:
        logger.info("Creating new project and adding photos...")
        doc.clear()
        chunk = doc.addChunk()
        chunk.addPhotos(image_files)
        doc.save(job.project_file)

    # --- Run steps ---
    try:
        with open(log_path, "a") as opf:
            # If user set start_step > end_step treat it as an explicit single-step override
            # (useful for testing: re-run a completed step once even if DB marks it completed).
            explicit_override = False
            if job.start_step is not None and job.end_step is not None and job.start_step > job.end_step:
                loop_end = job.start_step
                explicit_override = True
                logger.info(f"Explicit override detected: will run step {job.start_step} once (start>{job.end_step})")
            else:
                loop_end = min(MAX_STEP, job.end_step if job.end_step is not None else MAX_STEP)

            while job.start_step <= loop_end:
                step_number = job.start_step
                step_status = db.get_step_status(job.id, step_number)

                # Skip terminal states unless we're explicitly overriding to re-run this step
                if not explicit_override and step_status in ("completed", "skipped"):
                    logger.info(f"Step {step_number} already {step_status}. Skipping.")
                    job.end_step = max(job.end_step, step_number)
                    db.set_end_step(job.id, job.end_step)
                    job.start_step = step_number + 1
                    db.set_start_step(job.id, job.start_step)
                    continue

                # Resolve the step function
                try:
                    step_func = get_step_function(step_number)
                except Exception as e:
                    error_msg = f"Step {step_number} resolution/import failed: {e}"
                    logger.error(error_msg)
                    db.log_step(job.id, step_number, "error", error_msg)
                    db.update_job_status(job.id, "failed", error_msg)
                    raise

                # Mark running
                db.log_step(job.id, step_number, "running", f"Step {step_number} in progress")
                db.update_job_status(job.id, "running")

                step_result = None
                try:
                    # Execute the step
                    step_func(
                        job, doc, chunk, db, logger, opf, log_path, products_dir, job.project_file
                    )

                    # Success
                    db.log_step(job.id, step_number, "completed", f"Step {step_number} completed successfully")
                    doc.save(job.project_file)
                    step_result = "completed"

                except Exception as e:
                    error_msg = f"Step {step_number} failed: {e}"
                    logger.error(error_msg)

                    if step_number == 4:
                        # Explicitly mark skipped with a reason and continue
                        skip_msg = (
                            "Step 4 failed but will be skipped to continue pipeline."
                            f" Reason: {error_msg}"
                        )
                        db.log_step(job.id, step_number, "skipped", skip_msg)
                        warnings.append(skip_msg)
                        # Optional: save after skip if any state changed
                        try:
                            doc.save(job.project_file)
                        except Exception as save_err:
                            logger.warning(f"Save after skip failed: {save_err}")
                        step_result = "skipped"
                    else:
                        db.log_step(job.id, step_number, "error", error_msg)
                        db.update_job_status(job.id, "failed", error_msg)
                        # Leave start_step at the failed step so a retry will re-run this step
                        job.start_step = step_number
                        db.set_start_step(job.id, job.start_step)
                        step_result = "failed"
                        raise

                finally:
                    # Only advance progress if step completed successfully or was skipped
                    if step_result in ("completed", "skipped"):
                        job.end_step = max(job.end_step, step_number)
                        db.set_end_step(job.id, job.end_step)
                        job.start_step = step_number + 1
                        db.set_start_step(job.id, job.start_step)
                    else:
                        # Persist current progress without advancing (failed case)
                        db.set_end_step(job.id, job.end_step)
                        db.set_start_step(job.id, job.start_step)

            # Final job status
            final_status = "completed_with_warnings" if warnings else "completed"
            final_msg = "; ".join(warnings) if warnings else None
            db.update_job_status(job.id, final_status, final_msg)
            logger.info(f"Job {job.site_id} processing {final_status}.")
    except Exception:
        # Allow caller/scheduler to handle re-queue or notification as needed
        raise

    finally:
        doc = None  # Release file handles


if __name__ == "__main__":
    logger.info("Starting Metashape batch processing...")
    pending_jobs = db.get_pending_jobs()
    if not pending_jobs:
        logger.info("No pending jobs found.")
    else:
        job = pending_jobs[0]
        try:
            MetashapeProcess(job)
        except Exception as e:
            logger.error(f"Job {job.site_id} failed: {str(e)}")
