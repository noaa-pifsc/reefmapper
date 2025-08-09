import Metashape
import os

def run(job, doc, chunk, db, logger, opf):
    """
    Step 6: Build Dense Cloud
    """
    try:
        underline = 50 * '-'
        db.log_step(job.id, 6, "running", "Starting dense cloud construction")
        db.set_start_step(job.id, 6)
        print("Step 6. Build Dense Cloud")
        logger.info("Step 6. Build Dense Cloud")
        opf.write(f'\nBuild Dense Cloud\n{underline}\n')

        # Build depth maps
        db.log_step(job.id, 6, "info", "Building depth maps")
        chunk.buildDepthMaps(downscale=4, filter_mode=Metashape.MildFiltering)

        # Build point cloud
        db.log_step(job.id, 6, "info", "Building point cloud")
        chunk.buildPointCloud(point_colors=True, point_confidence=True)

        # Filter low confidence points
        db.log_step(job.id, 6, "info", "Filtering low confidence points")
        chunk.point_cloud.setConfidenceFilter(0, 1)
        chunk.point_cloud.removePoints(list(range(128)))  # Remove points with lowest confidence
        chunk.point_cloud.resetFilters()

        doc.save()
        db.log_step(job.id, 6, "completed", "Dense cloud construction complete")
        print("Step 6. Build Dense Cloud complete")
        logger.info("Step 6. Build Dense Cloud complete")
        db.set_end_step(job.id, 6)

    except Exception as e:
        msg = f"Step 6 failed: {str(e)}"
        logger.error(msg)
        db.log_step(job.id, 6, "error", msg)
        opf.write(msg + '\n')
        raise
