import Metashape
import math

def run(job, doc, chunk, db, logger, opf):
    """
    Step 4: Model scaling using markers and scalebars.
    """
    try:
        underline = 50 * '-'
        db.log_step(job.id, 4, "running", "Starting model scaling")
        db.set_start_step(job.id, 4)
        print("Step 4. Scaling")
        logger.info("Step 4. Scaling")
        opf.write(f'\nScaling\n{underline}\n')

        # Detect markers in the chunk
        db.log_step(job.id, 4, "info", "Detecting markers")
        chunk.detectMarkers()
        error_thresh = 0.4  # error threshold for marker reprojection error

        # Get list of markers
        MarkersList = list(chunk.markers)
        db.log_step(job.id, 4, "info", f"Initial marker count: {len(MarkersList)}")
        print("MarkersList:", MarkersList)
        logger.info(f"MarkersList: {MarkersList}")

        # Filter and validate markers
        for marker in MarkersList[:]:  # Create a copy for iteration
            print(underline, '\n', marker.label)
            if not marker:
                msg = f"Marker {marker.label} not found, skipping..."
                print(msg)
                logger.warning(msg)
                db.log_step(job.id, 4, "warning", msg)
                MarkersList.remove(marker)
                chunk.remove(marker)
                continue
            
            if not marker.position:
                msg = f"Marker {marker.label} not defined in 3D, skipping..."
                print(msg)
                logger.warning(msg)
                db.log_step(job.id, 4, "warning", msg)
                MarkersList.remove(marker)
                chunk.remove(marker)
                continue
            
            if marker.label not in valid_markers:
                msg = f"Marker {marker.label} not valid, skipping..."
                print(msg)
                logger.warning(msg)
                db.log_step(job.id, 4, "warning", msg)
                MarkersList.remove(marker)
                chunk.remove(marker)
                continue
            
            # Process marker errors
            opf.write(f'Marker: {marker.label}\n')
            pix_error = error_thresh
            while pix_error >= error_thresh:
                total = (0, 0)
                cam_err = []
                for camera in marker.projections.keys():
                    if not camera.transform:
                        continue
                    proj = marker.projections[camera].coord
                    reproj = camera.project(marker.position)
                    
                    # Fix: ensure both proj and reproj are 2D vectors
                    if len(proj) == 3:
                        proj = proj[:2]
                    if len(reproj) == 3:
                        reproj = reproj[:2]
                    
                    error = (proj - reproj).norm()
                    total = (total[0] + error**2, total[1] + 1)
                    cam_err.append((camera.label, error))
                
                pix_error = math.sqrt(total[0] / total[1])
                
                if pix_error >= error_thresh:
                    max_err = sorted(cam_err, key=lambda x: x[1], reverse=True)[0]
                    msg = f"Removed {max_err[0]} with error {max_err[1]:.4f}"
                    print(msg)
                    logger.info(msg)
                    db.log_step(job.id, 4, "info", msg)
                    opf.write(f'Removed {max_err[0]} with pix error {max_err[1]:.4f}\n')
                    for camera in marker.projections.keys():
                        if camera.label == max_err[0]:
                            marker.projections[camera] = None
                            
        # Add scalebars between every pair of markers (every other)
        db.log_step(job.id, 4, "info", "Adding scalebars")
        sb_dist = 0.25  # meters, distance for scalebars
        for i, marker in enumerate(MarkersList):
            if i % 2 > 0:
                sb = chunk.addScalebar(chunk.markers[i-1], chunk.markers[i])
                sb.reference.distance = sb_dist

        # Update transform to apply scaling
        db.log_step(job.id, 4, "info", "Updating transform")
        chunk.updateTransform()

        # Check scalebar errors and log
        for scalebar in chunk.scalebars:
            dist_source = scalebar.reference.distance
            if not dist_source:
                continue

            if isinstance(scalebar.point0, Metashape.Camera):
                if not (scalebar.point0.center and scalebar.point1.center):
                    continue
                dist_estimated = (scalebar.point0.center - scalebar.point1.center).norm() * chunk.transform.scale
            else:
                if not (scalebar.point0.position and scalebar.point1.position):
                    continue
                dist_estimated = (scalebar.point0.position - scalebar.point1.position).norm() * chunk.transform.scale

            dist_error = dist_estimated - dist_source
            msg = f'Scalebar {scalebar.label}: source={dist_source:.3f}, estimated={dist_estimated:.3f}, error={dist_error:.6f}'
            print(msg)
            logger.info(msg)
            db.log_step(job.id, 4, "info", msg)
            opf.write(f'{msg}\n')

            if dist_error >= 0.002:
                warn_msg = "Scalebar error too high - consider rerunning with lower error threshold"
                print(warn_msg)
                logger.warning(warn_msg)
                db.log_step(job.id, 4, "warning", warn_msg)
                opf.write(f'{warn_msg}\n')

            scalebar.reference.enabled = False

        doc.save()
        db.log_step(job.id, 4, "completed", "Scaling complete")
        print("Step 4. Scaling Finished")
        logger.info("Step 4. Scaling Finished")
        db.set_end_step(job.id, 4)

    except Exception as e:
        msg = f"Step 4 failed: {str(e)}"
        logger.error(msg)
        db.log_step(job.id, 4, "error", msg)
        opf.write(msg + '\n')
        raise
