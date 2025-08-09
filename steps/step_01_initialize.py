import os
import shutil
import Metashape

def run(job, doc, db, logger, os, shutil, Metashape, root_path=None, folder_name=None):
    """
    Step 1: Initialize project: create new or open existing Metashape project,
    prepare folders, open readme file, set logging.
    Returns updated doc, chunk, and opf objects.
    """
    try:
        if folder_name is None:
            folder_name = job.site_id
        if root_path is None:
            root_path = os.path.dirname(job.project_path)

        underline = 50 * '-'
        prod_path = os.path.join(root_path, 'Products_automation')
        psxfile = os.path.join(prod_path, folder_name + '.psx')

        if not os.path.exists(prod_path):
            os.mkdir(prod_path)
            logger.info(f"Created directory {prod_path}")
            db.set_start_step(job.id, 1)
        else:
            logger.info(f"Directory {prod_path} already exists")

        Metashape.app.settings.log_enable = True
        log_file = os.path.join(prod_path, folder_name + '_log.txt')
        Metashape.app.settings.log_path = log_file

        if job.start_step == 1:
            logger.info("Step 1 new project initialization")
            db.log_step(job.id, 1, "info", "Creating new project files")

            opf = open(os.path.join(prod_path, folder_name + '_readme.txt'), 'w')

            # Remove old files if they exist
            if os.path.exists(psxfile):
                os.remove(psxfile)
                logger.info(f"Deleted {psxfile}")

            if os.path.exists(log_file):
                os.remove(log_file)
                logger.info(f"Deleted {log_file}")

            file_path = os.path.join(prod_path, folder_name + '.files')
            if os.path.exists(file_path):
                shutil.rmtree(file_path, ignore_errors=True)
                logger.info(f"Deleted {file_path}")

            # Clear and create new document and chunk
            doc.clear()
            doc.save(psxfile)
            chunk = doc.addChunk()
            logger.info(f"Saved new project to {psxfile}")

        else:
            logger.info("Step 1 open existing project")
            db.log_step(job.id, 1, "info", "Opening existing project files")
            readme_path = os.path.join(prod_path, folder_name + '_readme.txt')
            if not os.path.exists(readme_path):
                opf = open(readme_path, 'w')
            else:
                opf = open(readme_path, 'a')

            if os.path.exists(psxfile):
                doc.open(psxfile)
                chunk = doc.chunk
            else:
                logger.warning(f"Project file {psxfile} not found. Creating new project.")
                doc.clear()
                chunk = doc.addChunk()

        chunk.camera_location_accuracy = Metashape.Vector((0.1, 0.1, 0.15))
        opf.write(f"Readme file for {folder_name}\n{underline}\n")

        db.log_step(job.id, 1, "completed", "Initialization complete")
        db.set_end_step(job.id, 1)
        logger.info("Step 1 complete")
        return doc, chunk, opf

    except Exception as e:
        msg = f"Step 1 failed: {str(e)}"
        logger.error(msg)
        db.log_step(job.id, 1, "error", msg)
        raise
