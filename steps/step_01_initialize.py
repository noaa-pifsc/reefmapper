import os
import Metashape

def run(job, doc, chunk, db, logger, opf, log_path, products_dir, project_file):
    """
    Step 1: Initialize project, create readme, and set up chunk.
    """
    try:
        folder_name = job.site_id
        underline = '-' * 50
        readme_path = os.path.join(products_dir, f"{folder_name}_readme.txt")

        # Ensure products directory exists
        os.makedirs(products_dir, exist_ok=True)
        logger.info(f"Ensured directory exists: {products_dir}")

        # Mark step start
        db.set_start_step(job.id, 1)
        db.log_step(job.id, 1, "running", "Initializing project")
        logger.info("Step 1. Initialize project")
        opf.write(f'\nInitialize project\n{underline}\n')

        # Create readme file if missing
        if not os.path.exists(readme_path):
            with open(readme_path, 'w') as readme_file:
                readme_file.write(f"Readme file for {folder_name}\n{underline}\n")
            logger.info(f"Created new readme file: {readme_path}")
        else:
            logger.info(f"Using existing readme file: {readme_path}")

        # ✅ Do NOT recreate log_path here — main.py already opened it
        logger.info(f"Using log file: {log_path}")

        # Initialize project file
        if not os.path.exists(project_file):
            logger.info(f"Creating new project at {project_file}")
            doc.save(project_file)
        else:
            logger.info(f"Opening existing project at {project_file}")
            doc.open(project_file)

        # Ensure chunk exists
        if not doc.chunks:
            chunk = doc.addChunk()
            logger.info("Created new chunk")
        else:
            chunk = doc.chunks[0]
            logger.info("Using existing chunk")

        # Save project state
        doc.save(project_file)

        # Mark step complete
        db.log_step(job.id, 1, "completed", "Initialization complete")
        db.set_end_step(job.id, 1)
        logger.info("Step 1 complete")

        return chunk, readme_path

    except Exception as e:
        msg = f"Step 1 failed: {str(e)}"
        logger.error(msg)
        db.log_step(job.id, 1, "error", msg)
        opf.write(msg + '\n')
        raise
