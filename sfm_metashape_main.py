import os
import sys
import json
import math
import csv
import Metashape
import shutil
import logging
import html2text
from datetime import datetime
from database.db_manager_mod import MyDatabaseManager, ProcessingJob


import logging
logger = logging.getLogger("reefmapper")  # same named logger

################################################################################################  
# Initialize database manager
db = MyDatabaseManager()
# batch_no = 1  # set number of batch to be processed (1-n)

################################################################################################  

from steps import get_step_function

def MetashapeProcess(job):
    logger.info(f"Processing job: {job.site_id}")

    doc = Metashape.Document()
    if os.path.exists(job.project_path):
        doc.open(job.project_path)
        chunk = doc.chunk
    else:
        doc.clear()
        chunk = doc.addChunk()

    with open(job.output_log_path, 'a') as opf:
        for step_number in range(job.start_step, job.end_step + 1):
            step_status = db.get_step_status(job.id, step_number)
            if step_status == "completed":
                logger.info(f"Step {step_number} already completed. Skipping.")
                continue

            try:
                step_func = get_step_function(step_number)
            except ImportError as e:
                error_msg = str(e)
                logger.error(error_msg)
                db.log_step(job.id, step_number, "error", error_msg)
                db.update_job_status(job.id, "failed", error_msg)
                raise

            try:
                db.log_step(job.id, step_number, "running", f"Step {step_number} in progress")
                step_func(job, doc, chunk, db, logger, opf)
                db.log_step(job.id, step_number, "completed", f"Step {step_number} completed successfully")
                doc.save(job.project_path)
            except Exception as e:
                error_msg = f"Step {step_number} failed: {str(e)}"
                logger.error(error_msg)
                db.log_step(job.id, step_number, "error", error_msg)

                # Special handling for step 4 failure: log and skip it, continue with next steps
                if step_number == 4:
                    logger.warning("Step 4 failed but will be skipped to continue pipeline.")
                    continue  # skip error raise, move on
                else:
                    db.update_job_status(job.id, "failed", error_msg)
                    raise

        db.update_job_status(job.id, "completed")
        logger.info(f"Job {job.site_id} processing completed.")
