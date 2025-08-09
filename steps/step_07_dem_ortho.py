import os
import json
import csv
import html2text
import Metashape

def run(job, doc, chunk, db, logger, opf):
    """
    Step 7: Build and export DEM and Orthomosaic
    """
    try:
        underline = 50 * '-'
        folder_name = job.site_id
        survey_year = job.survey_year
        prod_path = os.path.join(str(job.project_path).replace('\\\\', '\\'), 'Products_automation')
        
        db.log_step(job.id, 7, "running", "Starting DEM and orthomosaic generation")
        db.set_start_step(job.id, 7)
        print("Step 7. Build and export DEM and Orthomosaic")
        logger.info("Step 7. Build and export DEM and Orthomosaic")
        opf.write(f'\nBuild DEM and Orthomosaic\n{underline}\n')

        # Create ARC directory
        arc_path = os.path.join(prod_path, 'ARC')
        try:
            os.mkdir(arc_path)
        except FileExistsError:
            db.log_step(job.id, 7, "info", "ARC directory already exists")

        # Build DEM
        db.log_step(job.id, 7, "info", "Building DEM")
        chunk.buildDem(source_data=Metashape.PointCloudData)

        # Build Orthomosaic
        db.log_step(job.id, 7, "info", "Building orthomosaic")
        chunk.buildOrthomosaic(
            surface_data=Metashape.ElevationData,
            fill_holes=True,
            ghosting_filter=False,
            refine_seamlines=False,
            resolution=0.0005
        )

        # Export DEM
        dem_path = os.path.join(arc_path, f'{survey_year}_{folder_name}_dem.tif')
        db.log_step(job.id, 7, "info", f"Exporting DEM to {dem_path}")
        chunk.exportRaster(
            dem_path,
            resolution=0.001,
            save_world=True,
            source_data=Metashape.ElevationData
        )

        # Setup compression for orthomosaic
        compression = Metashape.ImageCompression()
        compression.tiff_compression = Metashape.ImageCompression.TiffCompressionLZW
        compression.jpeg_quality = 99
        compression.tiff_big = False
        compression.tiff_tiled = False
        compression.tiff_overviews = False

        # Export Orthomosaic
        ortho_path = os.path.join(arc_path, f'{survey_year}_{folder_name}_mos.tif')
        db.log_step(job.id, 7, "info", f"Exporting orthomosaic to {ortho_path}")
        chunk.exportRaster(
            ortho_path,
            resolution=0.0005,
            image_compression=compression,
            save_world=True,
            save_alpha=False,
            source_data=Metashape.OrthomosaicData
        )

        # Generate and export report (html)
        report_path = os.path.join(arc_path, f'{survey_year}_{folder_name}_rpt.html')
        db.log_step(job.id, 7, "info", "Generating processing report")
        chunk.exportReport(
            report_path,
            title=f'{survey_year}_{folder_name}',
            description="Processing Report"
        )

        # Convert report html to CSV
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
        cams_filename = os.path.join(proj_dir, proj_name + '.cams.xml')
        meta_filename = os.path.join(proj_dir, proj_name + '.meta.json')

        chunk.exportCameras(cams_filename)

        for cam in cams:
            key = cam.key
            path = cam.photo.path
            center = list(cam.center) if cam.center is not None else None
            agi_trans = cam.transform
            trans = [list(agi_trans.row(n)) for n in range(agi_trans.size[1])] if agi_trans else None
            outputs[key] = {'path': path, 'center': center, 'transform': trans}

        with open(meta_filename, 'w') as meta_file:
            json.dump({'cameras': outputs}, meta_file, indent=4)

        # Export point cloud if not already present
        pt_file = os.path.join(prod_path, f'{folder_name}.ply')
        if not os.path.isfile(pt_file):
            db.log_step(job.id, 7, "info", "Exporting point cloud")
            chunk.exportPointCloud(pt_file, source_data=Metashape.PointCloudData)

        doc.save()
        db.log_step(job.id, 7, "completed", "DEM and orthomosaic generation complete")
        print("Step 7. Build and export DEM and Orthomosaic Complete")
        logger.info("Step 7. Build and export DEM and Orthomosaic Complete")
        db.set_end_step(job.id, 7)

    except Exception as e:
        msg = f"Step 7 failed: {str(e)}"
        logger.error(msg)
        db.log_step(job.id, 7, "error", msg)
        opf.write(msg + '\n')
        raise
