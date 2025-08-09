import Metashape

def run(job, doc, chunk, db, logger, opf):
    """
    Step 5: Error reduction by filtering points with high reprojection error.
    """
    try:
        underline = 50 * '-'
        db.log_step(job.id, 5, "running", "Starting error reduction phase 2")
        db.set_start_step(job.id, 5)
        print("Step 5. Error reduction part 2")
        logger.info("Step 5. Error reduction part 2")
        opf.write(f'\nError Reduction\n{underline}\n')

        points = chunk.tie_points.points
        keep_percent = 90  # Keep top 90% points by reprojection error

        # Filter by Reprojection Error
        db.log_step(job.id, 5, "info", "Filtering reprojection error")
        f = Metashape.TiePoints.Filter()
        f.init(chunk, criterion=Metashape.TiePoints.Filter.ReprojectionError)

        list_values = f.values
        list_values_valid = [val for i, val in enumerate(list_values) if points[i].valid]
        list_values_valid.sort()

        target = int(len(list_values_valid) * keep_percent / 100)
        Rep_Err = list_values_valid[target]

        msg = f'Reprojection Error threshold to keep {keep_percent}%: {Rep_Err}'
        print(msg)
        logger.info(msg)
        db.log_step(job.id, 5, "info", msg)
        opf.write(f'Reprojection Error threshold to keep {keep_percent}%: {Rep_Err:.2f}\n')

        f.removePoints(Rep_Err)

        # Final camera optimization
        db.log_step(job.id, 5, "info", "Performing final camera optimization")
        chunk.optimizeCameras()

        doc.save()

        db.log_step(job.id, 5, "completed", "Error reduction complete")
        print("Step 5. Error reduction part 2 Done")
        logger.info("Step 5. Error reduction part 2 Done")
        db.set_end_step(job.id, 5)

    except Exception as e:
        msg = f"Step 5 failed: {str(e)}"
        logger.error(msg)
        db.log_step(job.id, 5, "error", msg)
        opf.write(msg + '\n')
        raise
