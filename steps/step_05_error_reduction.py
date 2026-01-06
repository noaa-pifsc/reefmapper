import Metashape

def run(job, doc, chunk, db, logger, opf, log_path, products_dir, project_file):
    """
    Step 5: Error reduction by filtering points with high reprojection error.
    """
    try:
        underline = '-' * 50
        db.log_step(job.id, 5, "running", "Starting error reduction phase 2")
        db.set_start_step(job.id, 5)
        logger.info("Step 5. Error reduction part 2")
        opf.write(f'\nError Reduction\n{underline}\n')

        points = chunk.tie_points.points
        keep_percent = 90

        # Filter by Reprojection Error
        db.log_step(job.id, 5, "info", "Filtering reprojection error")
        f = Metashape.TiePoints.Filter()
        f.init(chunk, criterion=Metashape.TiePoints.Filter.ReprojectionError)

        valid_values = sorted(val for i, val in enumerate(f.values) if points[i].valid)
        target_index = int(len(valid_values) * keep_percent / 100)
        Rep_Err = max(valid_values[target_index], 0.5)

        msg = f'Reprojection Error threshold to keep {keep_percent}%: {Rep_Err:.2f}'
        logger.info(msg)
        db.log_step(job.id, 5, "info", msg)
        opf.write(msg + '\n')

        f.removePoints(Rep_Err)

        # Final camera optimization
        db.log_step(job.id, 5, "info", "Performing final camera optimization")
        chunk.optimizeCameras()

        # Log remaining valid points
        valid_remaining = sum(1 for p in points if p.valid)
        opf.write(f"Remaining valid tie points: {valid_remaining}\n")
        logger.info(f"Remaining valid tie points: {valid_remaining}")

        doc.save(project_file)

        db.log_step(job.id, 5, "completed", "Error reduction complete")
        logger.info("Step 5. Error reduction part 2 Done")
        db.set_end_step(job.id, 5)

    except Exception as e:
        msg = f"Step 5 failed: {str(e)}"
        logger.error(msg)
        db.log_step(job.id, 5, "error", msg)
        opf.write(msg + '\n')
        raise
