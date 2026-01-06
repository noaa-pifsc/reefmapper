import Metashape
import math
from steps.utils import parse_marker_pairs, find_marker

def run(job, doc, chunk, db, logger, opf, log_path, products_dir, project_file):

    """
    Step 4: Model scaling using markers and scalebars.
    """
    try:
        underline = '-' * 50
        db.log_step(job.id, 4, "running", "Starting model scaling")
        db.set_start_step(job.id, 4)
        logger.info("Step 4. Scaling")
        opf.write(f'\nScaling\n{underline}\n')

        # Detect markers
        db.log_step(job.id, 4, "info", "Detecting markers")
        chunk.detectMarkers()
        error_thresh = 0.4

        # Filter markers
        raw_markers = list(chunk.markers)
        valid_markers = []

        for marker in raw_markers:
            if not marker or not marker.position:
                msg = f"Skipping marker {getattr(marker, 'label', 'unknown')} (missing or undefined)"
                logger.warning(msg)
                db.log_step(job.id, 4, "warning", msg)
                chunk.remove(marker)
                continue
            valid_markers.append(marker)

        db.log_step(job.id, 4, "info", f"Valid marker count: {len(valid_markers)}")

        # Reprojection error filtering
        for marker in valid_markers:
            opf.write(f'Marker: {marker.label}\n')
            pix_error = error_thresh

            while pix_error >= error_thresh:
                total_error = 0
                count = 0
                cam_err = []

                for camera, proj_data in marker.projections.items():
                    if not camera.transform:
                        continue

                    proj = proj_data.coord
                    reproj = camera.project(marker.position)

                    # Ensure both are 2D vectors
                    if len(proj) == 3:
                        proj = proj[:2]
                    if len(reproj) == 3:
                        reproj = reproj[:2]

                    error = (proj - reproj).norm()
                    total_error += error**2
                    count += 1
                    cam_err.append((camera.label, error))

                if count == 0:
                    break

                pix_error = math.sqrt(total_error / count)

                if pix_error >= error_thresh:
                    worst = max(cam_err, key=lambda x: x[1])
                    msg = f"Removed {worst[0]} with error {worst[1]:.4f}"
                    logger.info(msg)
                    db.log_step(job.id, 4, "info", msg)
                    opf.write(f'{msg}\n')

                    for camera in list(marker.projections.keys()):
                        if camera.label == worst[0]:
                            marker.projections[camera] = None

        # Add scalebars
        db.log_step(job.id, 4, "info", "Adding scalebars")
        sb_dist = 0.25
        for i in range(1, len(valid_markers), 2):
            sb = chunk.addScalebar(valid_markers[i - 1], valid_markers[i])
            sb.reference.distance = sb_dist

        # Apply scaling
        db.log_step(job.id, 4, "info", "Updating transform")
        chunk.updateTransform()

        # Check scalebar errors
        for scalebar in chunk.scalebars:
            dist_source = scalebar.reference.distance
            if not dist_source:
                continue

            if isinstance(scalebar.point0, Metashape.Camera):
                if not (scalebar.point0.center and scalebar.point1.center):
                    continue
                dist_est = (scalebar.point0.center - scalebar.point1.center).norm() * chunk.transform.scale
            else:
                if not (scalebar.point0.position and scalebar.point1.position):
                    continue
                dist_est = (scalebar.point0.position - scalebar.point1.position).norm() * chunk.transform.scale

            dist_error = dist_est - dist_source
            msg = f'Scalebar {scalebar.label}: source={dist_source:.3f}, estimated={dist_est:.3f}, error={dist_error:.6f}'
            logger.info(msg)
            db.log_step(job.id, 4, "info", msg)
            opf.write(f'{msg}\n')

            if dist_error >= 0.002:
                warn_msg = "Scalebar error too high - consider rerunning with lower error threshold"
                logger.warning(warn_msg)
                db.log_step(job.id, 4, "warning", warn_msg)
                opf.write(f'{warn_msg}\n')

            scalebar.reference.enabled = False

        doc.save(project_file)
        db.log_step(job.id, 4, "completed", "Scaling complete")
        logger.info("Step 4. Scaling Finished")
        db.set_end_step(job.id, 4)

    except Exception as e:
        msg = f"Step 4 failed: {str(e)}"
        logger.error(msg)
        db.log_step(job.id, 4, "error", msg)
        opf.write(msg + '\n')
        raise
