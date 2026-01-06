import Metashape

def run(job, doc, chunk, db, logger, opf, log_path, products_dir, project_file):
    """
    Step 6: Build Dense Cloud
    Hardened version matching single-script behavior.
    """
    try:
        underline = '-' * 50
        db.log_step(job.id, 6, "running", "Starting dense cloud construction")
        db.set_start_step(job.id, 6)
        logger.info("Step 6. Build Dense Cloud")
        opf.write(f'\nBuild Dense Cloud\n{underline}\n')

        # ---------------------------------------------------------
        # PHASE 1 — CLEANUP BLOCK (matches single-script behavior)
        # ---------------------------------------------------------
        try:
            logger.info("Cleaning up old products before dense cloud build...")
            db.log_step(job.id, 6, "info",
                        "Cleaning up old dense clouds, DEMs, and orthomosaics before rebuild")

            # --- Remove depth maps ---
            if getattr(chunk, "depth_maps", None):
                logger.info("Removing existing depth maps")
                if hasattr(chunk, "clearDepthMaps"):
                    chunk.clearDepthMaps()
                else:
                    try:
                        chunk.depth_maps.clear()
                    except Exception:
                        logger.warning("Unable to clear depth maps — continuing")

            # --- Remove dense cloud ---
            if getattr(chunk, "point_cloud", None):
                logger.info("Clearing existing dense cloud")
                try:
                    chunk.point_cloud.clear()
                except Exception:
                    chunk.point_cloud = None

            # --- Remove models ---
            if getattr(chunk, "models", None):
                for model in list(chunk.models):
                    logger.info(f"Removing existing model: {model.label}")
                    chunk.remove(model)

            # --- Remove DEMs ---
            if getattr(chunk, "elevation", None):
                for dem in list(chunk.elevation):
                    logger.info(f"Removing existing DEM: {dem.label}")
                    chunk.remove(dem)

            # --- Remove orthomosaics ---
            if getattr(chunk, "orthomosaics", None):
                for ortho in list(chunk.orthomosaics):
                    logger.info(f"Removing existing orthomosaic: {ortho.label}")
                    chunk.remove(ortho)

            # --- Save → reopen → save (critical memory flush) ---
            doc.save()
            doc.open(doc.path)
            doc.save()

            # Reload chunk reference
            chunk = doc.chunk

            logger.info("Cleanup complete — ready to rebuild dense cloud.")
            db.log_step(job.id, 6, "info", "Cleanup complete — ready to rebuild dense cloud.")

        except Exception as e:
            logger.warning(f"Cleanup before dense cloud failed: {e}")
            logger.info("Cleanup failed.")
            db.log_step(job.id, 6, "warning",
                        f"Cleanup before dense cloud failed: {e}")

        # ---------------------------------------------------------
        # PHASE 2 — BUILD DENSE CLOUD (matches single-script logic)
        # ---------------------------------------------------------
        try:
            db.log_step(job.id, 6, "running", "Building depth maps and dense cloud")
            logger.info("Building depth maps...")

            # Build depth maps
            db.log_step(job.id, 6, "info", "Building depth maps")
            chunk.buildDepthMaps(
                downscale=4,
                filter_mode=Metashape.MildFiltering
            )

            # Build dense cloud
            logger.info("Building dense cloud...")
            db.log_step(job.id, 6, "info", "Building point cloud")
            chunk.buildPointCloud(
                point_colors=True,
                point_confidence=True
            )

            # Filter low-confidence points
            logger.info("Filtering low-confidence points...")
            db.log_step(job.id, 6, "info", "Filtering low confidence points")

            chunk.point_cloud.setConfidenceFilter(0, 1)
            chunk.point_cloud.removePoints(list(range(128)))
            chunk.point_cloud.resetFilters()

            # Log remaining points
            remaining = chunk.point_cloud.point_count
            logger.info(f"Remaining dense cloud points: {remaining}")
            opf.write(f"Remaining dense cloud points: {remaining}\n")

        except Exception as e:
            logger.error(f"Error during dense cloud build: {e}")
            db.log_step(job.id, 6, "failed", f"Dense cloud build failed: {e}")
            raise

        # Save project
        doc.save(project_file)

        # Mark step complete
        db.log_step(job.id, 6, "completed", "Dense cloud construction complete")
        logger.info("Step 6. Build Dense Cloud complete")
        db.set_end_step(job.id, 6)

        # Human review pause (matches single-script)
        logger.info("Exiting after Step 6 for Human Review")
        db.log_step(job.id, 6, "info", "Exiting process after Step 6")
        db.update_job_status(job.id, "ready for review")

        try:
            opf.close()
        except Exception:
            pass

        return  # intentional pause for human review

    except Exception as e:
        msg = f"Step 6 failed: {str(e)}"
        logger.error(msg)
        db.log_step(job.id, 6, "error", msg)
        opf.write(msg + '\n')
        raise
