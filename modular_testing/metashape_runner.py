import logging
import os
import math
import shutil
import json
import html2text
import Metashape
import concurrent.futures

################################################################################################  

# Import your pipeline runner
from steps import run_metashape_pipeline

################################################################################################  

# Initialize database
from database.db_manager_mod import MyDatabaseManager

################################################################################################  
# Setup logging to both file and console

LOG_FORMAT = '%(asctime)s [%(levelname)s] %(message)s'
formatter = logging.Formatter(LOG_FORMAT)

logger = logging.getLogger("reefmapper")
logger.setLevel(logging.INFO)

# Log to file (overwrites each run)
file_handler = logging.FileHandler(r'T:\DataManagement\DataProjects\reefmapper\reefmapper\reefmapper_batch.log', mode='w')
file_handler.setFormatter(formatter)
logger.addHandler(file_handler)

# Log to console
console_handler = logging.StreamHandler()
console_handler.setFormatter(formatter)
logger.addHandler(console_handler)

################################################################################################  

def process_job_wrapper(job):
    try:
        with MyDatabaseManager() as db:
            run_metashape_pipeline(
                job=job,
                db=db,
                Metashape=Metashape,
                logger=logger,
                os=os,
                shutil=shutil,
                math=math,
                json=json,
                html2text=html2text,
            )
    except Exception as e:
        logger.error(f"Job {job.id} failed: {e}")

def main():
    with MyDatabaseManager() as db:
        jobs = db.get_pending_jobs()

    if not jobs:
        logger.info("No pending jobs found.")
        return

    max_workers = 3  # adjust for your CPU/Metashape limits
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        executor.map(process_job_wrapper, jobs)


if __name__ == "__main__":
    main()
