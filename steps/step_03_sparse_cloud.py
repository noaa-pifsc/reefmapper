import Metashape

def run(job, doc, chunk, db, logger, opf, log_path, products_dir=None, project_file=None):
    """
    Step 3: Sparse point cloud filtering in Metashape.
    """
    try:
        underline = '-' * 50
        db.log_step(job.id, 3, "running", "Starting sparse point cloud filtering")
        logger.info("Step 3. Sparse point cloud filtering")
        opf.write(f'\nSparse point cloud filtering\n{underline}\n')
        db.set_start_step(job.id, 3)

        keep_percent = 51
        tie_points = chunk.tie_points
        points = tie_points.points

        # --- Optimize cameras before filtering (1/3) ---
        db.log_step(job.id, 3, "info", "Optimizing cameras (1/3)")
        chunk.optimizeCameras(tiepoint_covariance=True)

        # --- Filter: Reconstruction Uncertainty ---
        db.log_step(job.id, 3, "info", "Filtering reconstruction uncertainty")
        f = Metashape.TiePoints.Filter()
        f.init(chunk, criterion=Metashape.TiePoints.Filter.ReconstructionUncertainty)

        valid_values = sorted(val for i, val in enumerate(f.values) if points[i].valid)
        target_index = int(len(valid_values) * keep_percent / 100)
        RecUncert = max(valid_values[target_index], 10)

        msg = f"Reconstruction Uncertainty threshold to keep {keep_percent}%: {RecUncert:.2f}"
        logger.info(msg)
        db.log_step(job.id, 3, "info", msg)
        opf.write(msg + '\n')

        f.removePoints(RecUncert)

        # --- Optimize cameras again (2/3) ---
        db.log_step(job.id, 3, "info", "Optimizing cameras (2/3)")
        chunk.optimizeCameras(tiepoint_covariance=True)

        # --- Filter: Projection Accuracy ---
        db.log_step(job.id, 3, "info", "Filtering projection accuracy")
        f = Metashape.TiePoints.Filter()
        f.init(chunk, criterion=Metashape.TiePoints.Filter.ProjectionAccuracy)

        valid_values = sorted(val for i, val in enumerate(f.values) if points[i].valid)
        target_index = int(len(valid_values) * keep_percent / 100)
        ProjAcc = max(valid_values[target_index], 2)

        msg = f"Projection Accuracy threshold to keep {keep_percent}%: {ProjAcc:.2f}"
        logger.info(msg)
        db.log_step(job.id, 3, "info", msg)
        opf.write(msg + '\n')

        f.removePoints(ProjAcc)

        # --- Final camera optimization (3/3) ---
        db.log_step(job.id, 3, "info", "Optimizing cameras (3/3)")
        chunk.optimizeCameras(tiepoint_covariance=True)

        # --- Summary of remaining points ---
        valid_points = sum(1 for p in points if p.valid)
        opf.write(f"Remaining valid tie points: {valid_points}\n")
        logger.info(f"Remaining valid tie points: {valid_points}")

        doc.save(project_file)

        db.log_step(job.id, 3, "completed", "Sparse point cloud filtering complete")
        logger.info("Step 3. Sparse point cloud filtering Finished")
        db.set_end_step(job.id, 3)

    except Exception as e:
        msg = f"Step 3 failed: {str(e)}"
        logger.error(msg)
        db.log_step(job.id, 3, "error", msg)
        opf.write(msg + '\n')
        raise
