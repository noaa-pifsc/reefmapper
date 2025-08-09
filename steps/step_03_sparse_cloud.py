import Metashape
import math

def run(job, doc, chunk, db, logger, opf):
    """
    Step 3: Sparse point cloud filtering in Metashape.
    """
    try:
        underline = 50 * '-'
        db.log_step(job.id, 3, "running", "Starting sparse point cloud filtering")
        print("Step 3. Sparse point cloud filtering")
        logger.info("Step 3. Sparse point cloud filtering")
        opf.write(f'\nSparse point cloud filtering\n{underline}\n')
        db.set_start_step(job.id,3)

        keep_percent = 51
        points = chunk.tie_points.points

        # Optimize cameras before filtering (1/3)
        db.log_step(job.id, 3, "info", "Optimizing cameras (1/3)")
        chunk.optimizeCameras(tiepoint_covariance=True)

        # Filter: Reconstruction Uncertainty
        db.log_step(job.id, 3, "info", "Filtering reconstruction uncertainty")
        f = Metashape.TiePoints.Filter()
        f.init(chunk, criterion=Metashape.TiePoints.Filter.ReconstructionUncertainty)

        list_values = f.values
        list_values_valid = [val for i, val in enumerate(list_values) if points[i].valid]
        list_values_valid.sort()

        target = int(len(list_values_valid) * keep_percent / 100)
        RecUncert = list_values_valid[target]

        msg = f'Reconstruction Uncertainty threshold to keep {keep_percent}%: {RecUncert}'
        print(msg)
        logger.info(msg)
        db.log_step(job.id, 3, "info", msg)

        if RecUncert < 10:
            RecUncert = 10
            msg = f'Reconstruction Uncertainty threshold set to {RecUncert}'
            print(msg)
            logger.info(msg)
            db.log_step(job.id, 3, "info", msg)

        opf.write(f'Reconstruction Uncertainty threshold to keep {keep_percent}%: {RecUncert:.2f}\n')
        f.removePoints(RecUncert)

        # Optimize cameras again (2/3)
        db.log_step(job.id, 3, "info", "Optimizing cameras (2/3)")
        chunk.optimizeCameras(tiepoint_covariance=True)

        # Filter: Projection Accuracy
        db.log_step(job.id, 3, "info", "Filtering projection accuracy")
        f = Metashape.TiePoints.Filter()
        f.init(chunk, criterion=Metashape.TiePoints.Filter.ProjectionAccuracy)

        list_values = f.values
        list_values_valid = [val for i, val in enumerate(list_values) if points[i].valid]
        list_values_valid.sort()

        target = int(len(list_values_valid) * keep_percent / 100)
        ProjAcc = list_values_valid[target]

        msg = f'Projection Accuracy threshold to keep {keep_percent}%: {ProjAcc}'
        print(msg)
        logger.info(msg)
        db.log_step(job.id, 3, "info", msg)

        if ProjAcc < 2:
            ProjAcc = 2
            msg = f'Projection Accuracy threshold set to {ProjAcc}'
            print(msg)
            logger.info(msg)
            db.log_step(job.id, 3, "info", msg)

        opf.write(f'Projection Accuracy threshold to keep {keep_percent}%: {ProjAcc:.2f}\n')
        f.removePoints(ProjAcc)

        # Final camera optimization (3/3)
        db.log_step(job.id, 3, "info", "Optimizing cameras (3/3)")
        chunk.optimizeCameras(tiepoint_covariance=True)

        doc.save()

        db.log_step(job.id, 3, "completed", "Sparse point cloud filtering complete")
        print("Step 3. Sparse point cloud filtering Finished")
        logger.info("Step 3. Sparse point cloud filtering Finished")
        db.set_end_step(job.id, 3)

    except Exception as e:
        msg = f"Step 3 failed: {str(e)}"
        logger.error(msg)
        db.log_step(job.id, 3, "error", msg)
        opf.write(msg + '\n')
        raise
