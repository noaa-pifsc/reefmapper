import os
import json
import csv
import html2text
import Metashape
import shutil
import fnmatch
import re

def run(job, doc, chunk, db, logger, opf, log_path, products_dir, project_file):
    """
    Step 7: Build and export DEM and Orthomosaic
    """
    try:
        underline = '-' * 50
        folder_name = job.site_id
        survey_year = job.survey_year


        db.log_step(job.id, 7, "running", "Starting DEM and orthomosaic generation")
        db.set_start_step(job.id, 7)
        logger.info("Step 7. Build and export DEM and Orthomosaic")
        opf.write(f'\nBuild DEM and Orthomosaic\n{underline}\n')

          # Create ARC directory
        arc_path = os.path.join(products_dir, 'ARC')
        try:
            os.mkdir(arc_path)
        except:
            db.log_step(job.id, 7, "info", "ARC directory already exists")
          
        # Check if DEM and orthomosaic already exist
        db.log_step(job.id, 7, "info", "Checking if DEM and ortho already exist")
        dem_path = os.path.join(arc_path, f'{survey_year}_{folder_name}_dem.tif')
        ortho_path = os.path.join(arc_path, f'{survey_year}_{folder_name}_mos.tif')
        
        if os.path.isfile(dem_path) and os.path.isfile(ortho_path):
            db.log_step(job.id, 7, "info", f"DEM and orthomosaic already exist, skipping build and export")
            logger.info(f"DEM and orthomosaic already exist, skipping build and export")
            # Export thumbnails
            thumbnail_dir = os.path.join(arc_path, "thumbnails")
            os.makedirs(thumbnail_dir, exist_ok=True)
            thumb_path = os.path.join(thumbnail_dir, "thumbnail_ortho.png")
                            
            try:
                if not chunk.orthomosaic:
                    raise RuntimeError("No orthomosaic exists in chunk")

                chunk.exportRaster(
                    path=thumb_path,
                    source_data=Metashape.OrthomosaicData,
                    image_format=Metashape.ImageFormat.ImageFormatPNG,
                    width=512,
                    height=512
                )

                db.log_step(job.id, 7, "info", "Thumbnails exported from existing orthomosaic")
                logger.info("Thumbnails exported from existing orthomosaic")

            except Exception as e:
                logger.warning(f"Failed to export thumbnails from existing orthomosaic: {e}")
                db.log_step(job.id, 7, "warning", f"Thumbnail export failed: {e}")

        else:
            # Clean up previous product files in the project's .files directory to avoid bloat
            try:
                proj_path = doc.path or project_file
                proj_dir, proj_name = os.path.split(proj_path)
                proj_name = proj_name[:-4] if proj_name.lower().endswith('.psx') else proj_name
                files_dir = os.path.join(proj_dir, proj_name + '.files')

                def _cleanup_product_files(files_dir):
                    if not os.path.isdir(files_dir):
                        return []
                    removed = []

                    # Remove directories which typically hold orthomosaic/dem caches
                    for entry in os.listdir(files_dir):
                        lname = entry.lower()
                        entry_path = os.path.join(files_dir, entry)

                        # heuristic: remove directories that match these substrings
                        if os.path.isdir(entry_path) and any(x in lname for x in ("orthomosaic", "dem", "elevation", "raster", "tiles", "mosaic", "dense", "point", "thumbnail")):
                            try:
                                shutil.rmtree(entry_path)
                                removed.append(entry_path)
                            except Exception as e:
                                logger.warning(f"Failed to remove {entry_path}: {e}")

                    # Remove files with common raster/aux extensions inside .files
                    for root, _, files in os.walk(files_dir):
                        for fname in files:
                            lfn = fname.lower()
                            if lfn.endswith(('.tif', '.ovr', '.vrt', '.aux', '.ply', '.jp2', '.png', '.jpg', '.jpeg')):
                                fpath = os.path.join(root, fname)
                                try:
                                    os.remove(fpath)
                                    removed.append(fpath)
                                except Exception as e:
                                    logger.warning(f"Failed to remove {fpath}: {e}")

                    return removed

                cleaned = _cleanup_product_files(files_dir)
                if cleaned:
                    msg = f"Cleaned {len(cleaned)} product file(s) from {files_dir}"
                    logger.info(msg)
                    db.log_step(job.id, 7, "info", msg)
            except Exception as e:
                logger.warning(f"Product cleanup failed or skipped: {e}")

            # Build DEM
            try:
                db.log_step(job.id, 7, "info", "Building DEM")
                #chunk.buildDem(source_data=Metashape.DenseCloudData,
                #                interpolation=Metashape.EnabledInterpolation
                #)
                chunk.buildDem(source_data=Metashape.PointCloudData)
            except Exception as e:
                logger.error(f"DEM build failed: {e}")
                db.log_step(job.id, 7, "error", f"DEM build failed: {e}")
                db.update_job_status(job.id, "failed", f"DEM build failed: {e}")
                raise
        
            # Build Orthomosaic
            try:
                db.log_step(job.id, 7, "info", "Building orthomosaic")
                chunk.buildOrthomosaic(surface_data=Metashape.ElevationData,
                                    fill_holes=True,
                                    ghosting_filter=False,
                                    refine_seamlines=False,
                                    resolution=0.0005)
            except Exception as e:
                logger.error(f"Orthomosaic build failed: {e}")
                db.log_step(job.id, 7, "error", f"Orthomosaic build failed: {e}")
                db.update_job_status(job.id, "failed", f"Orthomosaic build failed: {e}")
                raise
                
            # Export DEM
            try:
                dem_path = os.path.join(arc_path, f'{survey_year}_{folder_name}_dem.tif')
                db.log_step(job.id, 7, "info", f"Exporting DEM to {dem_path}")
                chunk.exportRaster(dem_path,
                                resolution=0.001,
                                save_world=True,
                                source_data=Metashape.ElevationData)
            except Exception as e:
                logger.error(f"DEM export failed: {e}")
                db.log_step(job.id, 7, "error", f"DEM export failed: {e}")
                db.update_job_status(job.id, "failed", f"DEM export failed: {e}")
                raise
            
            # Setup compression for orthomosaic
            compression = Metashape.ImageCompression()
            compression.tiff_compression = Metashape.ImageCompression.TiffCompressionLZW
            compression.jpeg_quality = 99
            compression.tiff_big = False
            compression.tiff_tiled = False
            compression.tiff_overviews = False
            
            # Export Orthomosaic
            try:
                ortho_path = os.path.join(arc_path, f'{survey_year}_{folder_name}_mos.tif')
                db.log_step(job.id, 7, "info", f"Exporting orthomosaic to {ortho_path}")
                chunk.exportRaster(ortho_path,
                                resolution=0.0005,
                                image_compression=compression,
                                save_world=True,
                                save_alpha=False,
                                source_data=Metashape.OrthomosaicData)
            except Exception as e:
                logger.error(f"Orthomosaic export failed: {e}")
                db.log_step(job.id, 7, "error", f"Orthomosaic export failed: {e}")
                db.update_job_status(job.id, "failed", f"Orthomosaic export failed: {e}")
                raise
        
            # Generate thumbnail of ortho
            db.log_step(job.id, 7, "info", "Generating orthomosaic thumbnail")
            # Export a small thumbnail TIFF directly from Metashape
            thumb_path = os.path.join(arc_path, folder_name + "_thumb.jpg")
            try:
                chunk.orthomosaic.exportRaster(
                    thumb_path,
                    image_format=Metashape.ImageFormatJPEG,
                    resolution=5.0  # meters/pixel or coarser for thumbnail
                )
                db.log_step(job.id, 7, "info", f"Thumbnail saved to {thumb_path}")
            except Exception as e:
                logger.warning(f"Failed to generate thumbnail: {e}")
                opf.write(f"Failed to generate thumbnail: {e}\n")       
            
            doc.save()

        # Generate and export report
        report_path = os.path.join(arc_path, f'{survey_year}_{folder_name}_rpt.html')
        db.log_step(job.id, 7, "info", "Generating processing report")
        chunk.exportReport(report_path,
                            title=f'{survey_year}_{folder_name}',
                            description="Processing Report")
        
        # Convert report to CSV
        with open(report_path) as f:
            html = f.read()
        h = html2text.HTML2Text()
        t = h.handle(html)
        lines = t.splitlines()
        
        with open(os.path.join(arc_path, f'{survey_year}_{folder_name}_rpt.csv'), 'w', newline='') as file:
            writer = csv.writer(file)
            for line in lines:
                words = line.strip().split('|')
                writer.writerow(words)
        
        # Extract camera metadata
        db.log_step(job.id, 7, "info", "Extracting camera metadata")
        cams = chunk.cameras
        proj_path = doc.path
        proj_dir, proj_name = os.path.split(proj_path)
        proj_name = proj_name[:-4]
        
        outputs = {}
        cams_filename = proj_dir + '/' + proj_name + '.cams.xml'
        meta_filename = proj_dir + '/' + proj_name + '.meta.json'
        
        chunk.exportCameras(cams_filename)
        
        # export cam x,y,z locations csv - add frp, Dama 12-16-2025
        chunk.exportReference(os.path.join(arc_path, survey_year+'_'+folder_name+'_cams.csv'), format=Metashape.ReferenceFormatCSV,
                            items=Metashape.ReferenceItemsCameras,  # or Metashape.ReferenceItemsMarkers
                            columns='nuvw',                         # Specifies Label, X, Y, Z columns
                            delimiter =',') 
        
        for cam in cams:
            key = cam.key
            path = cam.photo.path
            center = list(cam.center) if cam.center is not None else None
            agi_trans = cam.transform
            trans = [list(agi_trans.row(n)) for n in range(agi_trans.size[1])] if agi_trans else None
            outputs[key] = {'path': path, 'center': center, 'transform': trans}
        
        with open(meta_filename, 'w') as meta_file:
            json.dump({'cameras': outputs}, meta_file, indent=4)
        
        # Export point cloud
        pt_file = os.path.join(arc_path, f'{folder_name}.ply')
        if not os.path.isfile(pt_file):
            db.log_step(job.id, 7, "info", "Exporting point cloud")
            chunk.exportPointCloud(pt_file, source_data=Metashape.PointCloudData)
        
        doc.save()
        db.log_step(job.id, 7, "completed", "DEM and orthomosaic generation complete")
        print("Step 7. Build and export DEM and Orthomosaic Complete")
        logger.info("Step 7. Build and export DEM and Orthomosaic Complete")
        db.set_end_step(job.id, 7)
        db.update_priority(job.id, '')
    
        # On successful completion
        logger.info(f"Completed processing {folder_name}")
        db.update_job_status(job.id, "completed")

    except Exception as e:
        error_msg = f"Error processing {folder_name}: {str(e)}"
        logger.error(error_msg)
        db.log_step(job.id, -1, "error", error_msg)
        # mark job as failed and return so the main loop can continue with other jobs
        try:
            db.update_job_status(job.id, "failed", error_msg)
        except Exception:
            # ensure we don't raise from the error handler
            logger.exception("Failed to update job status after error")
        return