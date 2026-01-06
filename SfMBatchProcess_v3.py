# Auto batch process for Agisoft Metashape
# Following Structure-from-Motion workflow
# Database-driven version
#
# Author(s): F. Lichowski, M. Akrdige, D. Torres-Pulliza
#
# Changelog - Updates to match Metashape Python API updates
# Changes as follows:  
#     chunk.analyzePhotos to chunk.analyzeImages
#     PointCloud to TiePoints
#     point_cloud to tie_points
#     buildDenseCloud to buildPointCloud
#     Chunk.exportPoints() to exportPointCloud()
#     DataSource.DenseCloudData to DataSource.PointCloudData
#
# reference: https://www.agisoft.com/forum/index.php?topic=9578.15
#            https://github.com/gisportsmouth/Agisoft-Metashape-Automation-Script
#            https://www.agisoft.com/pdf/metashape_python_api_2_1_0.pdf
################################################################################################
import os
import sys
import json
import math
import csv
import Metashape
import re
import shutil
import logging
import html2text
from datetime import datetime
from database.db_manager import MyDatabaseManager, ProcessingJob
from pathlib import Path

################################################################################################  
# Setup logging to both file and console
LOG_FORMAT = '%(asctime)s [%(levelname)s] %(message)s'
logger = logging.getLogger("reefmapper")
logger.setLevel(logging.INFO)
formatter = logging.Formatter(LOG_FORMAT)

# --- Prepare unique log file ---
timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
log_dir = Path(r"C:\Users\PICHLMRUser\Desktop\reefmapper\logs")
log_dir.mkdir(parents=True, exist_ok=True)  # ensure directory exists
log_file = log_dir / f"reefmapper_{timestamp}.log"

# --- Set up logging ---
formatter = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s")
file_handler = logging.FileHandler(log_file, mode="w", encoding="utf-8")
file_handler.setFormatter(formatter)

logger = logging.getLogger()  # or use your existing logger
logger.addHandler(file_handler)
logger.setLevel(logging.INFO)

logger.info(f"Logging to {log_file}")

# Log to console
console_handler = logging.StreamHandler()
console_handler.setFormatter(formatter)
logger.addHandler(console_handler)


################################################################################################
# Define function to rename images based on sequence and datetime

def rename_images_by_datetime(folder, start_index=1, logger=None):
    """
    Rename images in a folder based on datetime (including milliseconds) and a sequential index.
    Ensures proper sorting for Metashape.
    
    Parameters:
        folder (str): Path to folder containing images.
        start_index (int): Starting number for sequence.
        logger (logging.Logger): Optional logger.
    """
    # Collect image files (common image extensions)
    image_exts = (".jpg", ".jpeg",".JPG",".JPEG")

    # Collect image files
    images = [f for f in os.listdir(folder) if f.lower().endswith(image_exts)]
    if not images:
        if logger:
            logger.warning(f"No images found in folder: {folder}")
        return

    # Gather datetime info
    image_info = []
    for fname in images:
        path = os.path.join(folder, fname)
        try:
            dt = datetime.fromtimestamp(os.path.getmtime(path))
            ms = dt.microsecond // 1000  # milliseconds
            image_info.append((fname, dt, ms))
        except Exception as e:
            if logger:
                logger.warning(f"Failed to read timestamp for {fname}: {e}")

    # Sort by datetime and milliseconds
    image_info.sort(key=lambda x: (x[1], x[2]))

    # Rename sequentially
    for idx, (fname, dt, ms) in enumerate(image_info, start=start_index):
        ext = os.path.splitext(fname)[1].lower()
        new_name = f"{dt.strftime('%Y%m%d_%H%M%S')}_{ms:03d}_{idx:04d}{ext}"
        old_path = os.path.join(folder, fname)
        new_path = os.path.join(folder, new_name)

        # Avoid overwriting or renaming to same name
        if old_path != new_path:
            try:
                shutil.move(old_path, new_path)
                if logger:
                    logger.info(f"Renamed: {fname} -> {new_name}")
            except Exception as e:
                if logger:
                    logger.error(f"Failed to rename {fname}: {e}")

################################################################################################  
# Initialize database manager
db = MyDatabaseManager()
# batch_no = 1  # set number of batch to be processed (1-n)

################################################################################################  
# function to make list of all photos in root_path
PHOTO_EXTENSIONS = {'.jpg', '.jpeg', '.JPG', '.JPEG'}
def getPhotoList(root_path, photoList):
    files = [f for f in os.listdir(root_path) if os.path.isfile(os.path.join(root_path, f))]
    for photo in files:
        if os.path.splitext(photo)[1] in PHOTO_EXTENSIONS:
            add_path = os.path.join(root_path, photo)
            photoList.append(add_path)           

def find_marker(label, chunk):
    for marker in chunk.markers:
        if label == marker.label:
            return marker
    return None

def parse_marker_pairs(job: ProcessingJob) -> list:
    """Parse marker pairs from job and return valid markers list."""
    valid_markers = []
    for pair in [job.marker_pair1, job.marker_pair2, job.marker_pair3, job.marker_pair4]:
        if pair and pair.lower() != 'na':
            t1, t2 = pair.split(',')
            valid_markers.append('target ' + t1.strip())
            valid_markers.append('target ' + t2.strip())
    return valid_markers


def orient_model_by_marker_depth(chunk, job, logger):
    """
    Rotates the chunk so that:
    - lowest depth marker pair is on the LEFT (-X)
    - highest depth marker pair is on the RIGHT (+X)
    """
    logger.info("Rotating and orienting model based on lowest/highest scalebars")

    # ---- Build marker depth lookup ----
    marker_depths = {}
    for idx, pair in enumerate(
        (job.marker_pair1, job.marker_pair2, job.marker_pair3, job.marker_pair4),
        start=1
    ):
        depth = getattr(job, f"marker_pair{idx}_depth", None)
        if not pair or pair.lower() == "na" or depth is None:
            continue

        try:
            t1, t2 = pair.split(",")
            marker_depths[f"target {t1.strip()}"] = depth
            marker_depths[f"target {t2.strip()}"] = depth
        except Exception:
            continue

    # ---- Collect marker centers ----
    marker_centers = []
    for marker in chunk.markers:
        if marker.label in marker_depths and marker.position:
            marker_centers.append((marker_depths[marker.label], marker.position))

    if len(marker_centers) < 2:
        logger.warning("Not enough valid markers to orient model, skipping rotation")
        return

    # ---- Find lowest & highest depth groups ----
    marker_centers.sort(key=lambda x: x[0])  # shallow → deep
    lowest_depth = marker_centers[0][0]
    highest_depth = marker_centers[-1][0]

    low_pts = [p for d, p in marker_centers if d == lowest_depth]
    high_pts = [p for d, p in marker_centers if d == highest_depth]

    # ---- Compute midpoints ----
    low_center = sum(low_pts, Metashape.Vector((0, 0, 0))) / len(low_pts)
    high_center = sum(high_pts, Metashape.Vector((0, 0, 0))) / len(high_pts)

    # ---- Direction vector (low → high) ----
    direction = high_center - low_center
    if direction.norm() < 1e-6:
        logger.info("Markers are coincident or too close; skipping rotation")
        return
    direction = Metashape.Vector(direction)
    direction.normalize()

    # Target direction: +X axis
    target = Metashape.Vector((1, 0, 0))

    # ---- Compute rotation safely ----
    dot = max(-1.0, min(1.0, direction * target))
    angle = math.acos(dot)

    try:
        axis = Metashape.cross(direction, target)
        if axis.norm() < 1e-6:
            logger.info("Model already oriented correctly (axis too small)")
            return
        axis.normalize()
    except Exception as e:
        logger.error(f"Cannot compute cross product: direction={direction}, target={target}, error={e}")
        return

    R = Metashape.Matrix.Rotation(axis, angle)

    # ---- Apply rotation while preserving scale and translation ----
    T = chunk.transform.matrix
    chunk.transform.matrix = Metashape.Matrix.Translation(T.translation()) * R * Metashape.Matrix.Diagonal(T.scale())
    chunk.updateTransform()

    logger.info("Model rotation/orientation updated: lowest depth left, highest depth right")



# helper script to prevent .files bloat before rebuild
def cleanup_rebuild_products(doc, chunk, logger):
    """
    Safely remove all derived products that cause .files bloat
    without touching alignment, markers, or cameras.
    """

    logger.info("Cleaning derived products before rebuild")

    # Depth maps (largest bloat source)
    try:
        if hasattr(chunk, "clearDepthMaps"):
            chunk.clearDepthMaps()
            logger.info("Depth maps cleared")
        elif hasattr(chunk, "depth_maps") and chunk.depth_maps:
            chunk.depth_maps.clear()
            logger.info("Depth maps cleared (legacy)")
    except Exception as e:
        logger.warning(f"Depth map cleanup failed: {e}")

    # Dense cloud
    if chunk.point_cloud:
        chunk.point_cloud = None
        logger.info("Dense cloud removed")

    # Mesh models
    if getattr(chunk, "models", None):
        for model in list(chunk.models):
            chunk.remove(model)
            logger.info(f"Removed model: {model.label}")

    # DEMs
    if getattr(chunk, "elevations", None):
        for dem in list(chunk.elevations):
            chunk.remove(dem)
            logger.info(f"Removed DEM: {dem.label}")

    # Orthomosaics
    if getattr(chunk, "orthomosaics", None):
        for ortho in list(chunk.orthomosaics):
            chunk.remove(ortho)
            logger.info(f"Removed orthomosaic: {ortho.label}")

    doc.save()

def MetashapeProcess(job: ProcessingJob):

    """Main function for automatic batch processing using job information from database."""
    try:
        folder_name = job.site_id
        project_path = str(job.project_path)
        root_path = project_path.replace('\\\\', '\\')  # Handle double backslashes
        start_step = job.start_step
        end_step = job.end_step
        valid_markers = parse_marker_pairs(job)
        quality = job.quality
        survey_year = job.survey_year
        
        # logger.info(f"Renaming images for: {folder_name}")
        # Rename images based on sequence and datetime
       #  rename_images_by_datetime(root_path, start_index=1, logger=logger)

        
        logger.info(f"Starting batch: {folder_name} (steps {start_step}-{end_step})")
        db.log_step(job.id, 0, "started", f"Starting batch processing for {folder_name}")
        underline = 50*'-'
         # create new Products folder for project files
        prod_path = os.path.join(root_path, 'Products_automation')
        psxfile = os.path.join(prod_path, folder_name+'.psx')
        print(root_path)
        print(prod_path)

        try:
            os.mkdir(prod_path)
            db.set_start_step(job.id, 1)

        except:
            print(prod_path, 'directory already exists')    
            logger.info(f'{prod_path} directory already exists')

        #==============================================================================================
        # Step 1: Initialize
        db.log_step(job.id, 1, "running", "Initializing project")
        print("Step 1. Initialize running")
        logger.info("Step 1. Initialize running")
        db.update_job_status(job.id,'running')

        # construct the document class (with no chunks)
        doc = Metashape.app.document

        if start_step > 1:
                    # Continuing an existing job
                    if os.path.exists(psxfile):
                        try:
                            doc.open(psxfile)
                            chunk = doc.chunk
                            logger.info(f"Reopened existing project: {psxfile}")
                        except Exception as e:
                            logger.error(f"Error reopening {psxfile}: {e}")
                            doc.clear()
                            chunk = doc.addChunk()
                            logger.warning("Created new project due to open failure")
                    else:
                        logger.warning(f"No project found at {psxfile}, creating new.")
                        doc.clear()
                        chunk = doc.addChunk()
                        chunk.label = f"{job.site_id}_chunk"
                        #may wan tto remove the chunk.label text above and below

                        # Save the empty project file
                        try:
                            doc.save(psxfile)
                            logger.info(f"Created new .psx project: {psxfile}")
                            db.log_step(job.id, 1, "info", f"Created new project at {psxfile}")
                        except Exception as e:
                            logger.error(f"Failed to create project file: {e}")
                            db.log_step(job.id, 1, "error", f"Failed to save project file: {e}")
                            db.update_status(job.id, "failed")
                            raise
        else:
            # Starting fresh
            doc.clear()
            chunk = doc.addChunk()
            doc.save(psxfile)
            logger.info(f"Created new project: {psxfile}")
        
        
        # write log and readme text file with processing information
        Metashape.app.settings.log_enable = True
        log_file = os.path.join(prod_path, folder_name+'_log.txt')
        Metashape.app.settings.log_path = log_file
        
        if start_step==1 and end_step ==0:
            print('Step 1 new')
            logger.info('Step 1 new')
            db.log_step(job.id, 1, "info", "Creating new project files")

            # start new readme file
            opf = open(os.path.join(prod_path, folder_name+'_readme.txt'), 'w')       

            # remove .psx file, if it exists
            if os.path.exists(psxfile):
                os.remove(psxfile)
                print('Deleted', psxfile)
                logger.info(f'Deleted {psxfile}')
                
            # remove log file, if it exists
            if os.path.exists(log_file):          
                os.remove(log_file)
                print('Deleted', log_file)
                logger.info(f'Deleted {log_file}')
            
            # remove .files folder, if it exists
            file_path = os.path.join(prod_path, folder_name+'.files')
            if os.path.exists(file_path):
                shutil.rmtree(file_path, ignore_errors=True)
                print('Deleted', file_path)
                logger.info(f'Deleted {file_path}')
                           
            # save .psx file and add new chunk     
            doc.save(psxfile)
            chunk = doc.addChunk()
            print('Saved project to: ' + psxfile)
            logger.info(f'Saved project to: {psxfile}')
            db.set_end_step(job.id,1)
        
        else:
            print('Step 1 existing')
            logger.info('Step 1 existing')
            db.log_step(job.id, 1, "info", "Opening existing project files")
             # for start_step>1: open existing readme and .psx file and load chunk

            readme_path = os.path.join(prod_path, folder_name+'_readme.txt')
            if not os.path.exists(readme_path):
                opf = open(readme_path, 'w')  # create new readme if missing
            else:
                opf = open(readme_path, 'a')
            
            if os.path.exists(psxfile):
                doc.open(psxfile)
                chunk = doc.chunk
            else:
                logger.warning(f"Project file {psxfile} not found. Creating new project.")
                doc.clear()
                chunk = doc.addChunk()
        
        chunk.camera_location_accuracy = Metashape.Vector((0.1, 0.1, 0.15))
        opf.write('Readme file for {0}\n{1}\n'.format(folder_name, underline))
        print("Step 1. Initialize Finished")          
        logger.info("Step 1. Initialize Finished")
        db.log_step(job.id, 1, "completed", "Initialization complete")


        #==============================================================================================
        # Step 2: Add and align photos
        if (start_step ==1 and end_step ==1) or (start_step ==2 and job.status == 'failed'):
            db.log_step(job.id, 2, "running", "Adding and aligning photos")
            print("Step 2. Add and align photos")
            logger.info("Step 2. Add and align photos")
            opf.write('\nAlign photos\n{0}\n'.format(underline))

            try:
                ## get photo list
                photoList = []
                extra = ''
                getPhotoList(os.path.join(root_path, extra), photoList)
                n = len(photoList)
                db.set_start_step(job.id, 2)

                if n == 0:
                    raise Exception(f'No photos found in {os.path.join(root_path, extra)}')
                if len(chunk.cameras) == 0:
                ## add photos
                    db.log_step(job.id, 2, "info", f"Adding {n} photos")
                    print(f'Adding {n} photos')
                    logger.info(f'Adding {n} photos')
                    opf.write(f'Adding {n} photos\n')
                    chunk.addPhotos(photoList)
                    doc.save(psxfile)
                else:
                    logger.info(f"{len(chunk.cameras)} photos already in chunk, skipping addPhotos()")  
                    opf.write(f"{len(chunk.cameras)} photos already in chunk, skipping addPhotos()\n")  
                ## estimate image quality and disable poor quality images
                db.log_step(job.id, 2, "info", "Analyzing image quality")
                chunk.analyzeImages(chunk.cameras)
                doc.save(psxfile)

                qualities = []
                bad_quality = 0
                quality_log = [0] * n  # 0 = good, 1 = poor

                for qc, camera in enumerate(chunk.cameras):
                    img_quality = float(camera.meta['Image/Quality'])
                    qualities.append(img_quality)
                    print('Image Quality:', img_quality)
                    quality = 0.5

                    if img_quality < quality and quality_log[qc - 1] == 0:
                        bad_quality += 1
                        quality_log[qc] = 1
                        camera.enabled = False
                        print(qc, camera.label, img_quality)
                        logger.info(f'Image {camera.label} with quality {img_quality:.3f} disabled')
                        opf.write(f'Image {camera.label} with quality {img_quality:.3f} disabled\n')

                # Calculate and update average image quality
                if qualities:
                    avg_quality = sum(qualities) / len(qualities)
                    db.update_quality(job.id, avg_quality)
                    logger.info(f"Average image quality: {avg_quality:.3f}")
                    opf.write(f"Average image quality: {avg_quality:.3f}\n")
                else:
                    logger.warning("No image qualities found to compute average.")
                    avg_quality = None

                # Logging poor image results
                n_enabled = n - bad_quality
                msg = f'{bad_quality} photos ({round(bad_quality / n * 100, 1)}%) below quality threshold {quality}'
                db.log_step(job.id, 2, "info", msg)
                print(msg)
                logger.info(msg)
                opf.write(msg + '\n')

                ## perform image matching and alignment
                db.log_step(job.id, 2, "info", "Matching and aligning photos")
                chunk.matchPhotos(generic_preselection=True, reference_preselection=False,
                                filter_mask=False, keypoint_limit=40000, tiepoint_limit=0)
                chunk.alignCameras(adaptive_fitting=False)

                ## verify alignment
                thresh_align = 15
                counter = sum(1 for camera in chunk.cameras if camera.transform)
                msg = f'Enabled: {n_enabled}, Aligned: {counter}'
                db.log_step(job.id, 2, "info", msg)
                print(msg)
                logger.info(msg)
                opf.write(msg + '\n')

                if n_enabled - counter == 0:
                    db.log_step(job.id, 2, "info", "All enabled images aligned successfully")
                elif n_enabled - counter > 0 and n_enabled - counter <= thresh_align:
                    db.log_step(job.id, 2, "warning", f"{n_enabled - counter} images failed to align")
                else:
                    raise Exception("Too many images failed to align")

                ## backup point cloud
                doc.save()
                shutil.copyfile(psxfile, os.path.join(prod_path, folder_name + '_bkup.psx'))
                if os.path.exists(os.path.join(prod_path, folder_name + '_bkup.files')):
                    shutil.rmtree(os.path.join(prod_path, folder_name + '_bkup.files'), ignore_errors=True)
                shutil.copytree(os.path.join(prod_path, folder_name + '.files'),
                                os.path.join(prod_path, folder_name + '_bkup.files'))

                db.log_step(job.id, 2, "completed", "Photo alignment complete")
                print("Step 2. Add and align photos finished")
                logger.info("Step 2. Add and align photos finished")

                db.set_end_step(job.id, 2)

            except Exception as e:
                msg = f"Step 2 failed: {str(e)}"
                logger.error(msg)
                db.update_job_status(job.id,'failed')
                db.log_step(job.id, 2, "error", msg)
                opf.write(msg + '\n')
                try:
                    opf.close()
                except Exception:
                    pass
                raise Exception(f"Step 2 failed: {str(e)}")  # Raise to be caught by main()
        #==============================================================================================
        # Step 3: Sparse point cloud filtering
        if (start_step ==2 and end_step ==2) or (start_step ==3 and job.status == 'failed'):
            db.log_step(job.id, 3, "running", "Starting sparse point cloud filtering")
            print("Step 3. Sparse point cloud filtering")
            logger.info("Step 3. Sparse point cloud filtering")
            opf.write('\nSparse point cloud filtering\n{0}\n'.format(underline))
            db.set_start_step(job.id, 3)

            keep_percent = 51
            points = chunk.tie_points.points
            
            # Optimize cameras
            db.log_step(job.id, 3, "info", "Optimizing cameras (1/3)")
            chunk.optimizeCameras(tiepoint_covariance=True)
            
            # Filter: Reconstruction Uncertainty
            db.log_step(job.id, 3, "info", "Filtering reconstruction uncertainty")
            f = Metashape.TiePoints.Filter()
            f.init(chunk, criterion = Metashape.TiePoints.Filter.ReconstructionUncertainty)
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
            
            # Optimize cameras again
            db.log_step(job.id, 3, "info", "Optimizing cameras (2/3)")
            chunk.optimizeCameras(tiepoint_covariance=True)
            
            # Filter: Projection Accuracy
            db.log_step(job.id, 3, "info", "Filtering projection accuracy")
            f = Metashape.TiePoints.Filter()
            f.init(chunk, criterion = Metashape.TiePoints.Filter.ProjectionAccuracy)
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
            
            # Final camera optimization
            db.log_step(job.id, 3, "info", "Optimizing cameras (3/3)")
            chunk.optimizeCameras(tiepoint_covariance=True)
            doc.save()
            

            db.log_step(job.id, 3, "completed", "Sparse point cloud filtering complete")
            print("Step 3. Sparse point cloud filtering Finished")
            logger.info("Step 3. Sparse point cloud filtering Finished")
            
            db.set_end_step(job.id, 3)

        #==============================================================================================
        # Step 4: Scaling
        if (start_step ==3 and end_step ==3) or (start_step ==4 and job.status == 'failed'):
            db.log_step(job.id, 4, "running", "Starting model scaling")
            print("Step 4. Scaling")
            print("Metashape version:", Metashape.app.version)
            db.log_step(job.id, 4, "running", f"Metashape version: {Metashape.app.version}")
            logger.info("Step 4. Scaling")
            opf.write('\nScaling\n{0}\n'.format(underline))
            db.set_start_step(job.id,4)
            

            # Detect markers
            db.log_step(job.id, 4, "info", "Detecting markers")
            chunk.detectMarkers()
            error_thresh = 0.4  # marker error threshold
            
            
            # Process markers
            MarkersList = list(chunk.markers)
            db.log_step(job.id, 4, "info", f"Initial marker count: {len(MarkersList)}")
            print("MarkersList:", MarkersList)
            logger.info(f"MarkersList: {MarkersList}")
            
            # Validate and filter markers
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
           
            # Detect existing scalebars
            existing_scalebars = list(chunk.scalebars)

            if existing_scalebars:
                logger.info(f"{len(existing_scalebars)} scalebars already exist — reusing them")
                db.log_step(job.id, 4, "info", "Reusing existing scalebars")
            else:
                logger.info("No existing scalebars found — creating new ones")

            if not existing_scalebars:
                db.log_step(job.id, 4, "info", "Adding scalebars")
                sb_dist = 0.25
                for i, marker in enumerate(MarkersList):
                    if i % 2 > 0:
                        sb = chunk.addScalebar(chunk.markers[i-1], chunk.markers[i])
                        sb.reference.distance = sb_dist
            
            # Enable markers
            for marker in MarkersList:
                marker.enabled = True
                db.log_step(job.id, 4, "info", f"Enabled markers")

            # Enable all scalebars explicitly
            enabled_count = 0
            for sb in chunk.scalebars:
                if sb.reference:
                    sb.reference.enabled = True
                    enabled_count += 1
                    logger.info(f"Enabled {enabled_count} scalebars for scaling")
                    db.log_step(job.id, 4, "info", f"Enabled {enabled_count} scalebars for scaling")

            # Update transform
            db.log_step(job.id, 4, "info", "Updating transform")

            chunk.updateTransform()
            
            # Check scalebar errors
            for scalebar in chunk.scalebars:
                dist_source = scalebar.reference.distance
                if not dist_source:
                    continue
                
                if type(scalebar.point0) == Metashape.Camera:
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
                    msg = "Scalebar error too high - consider rerunning with lower error threshold"
                    print(msg)
                    logger.warning(msg)
                    db.log_step(job.id, 4, "warning", msg)
                    opf.write(f'{msg}\n')
                
               #  scalebar.reference.enabled = False
            
            # --- Log marker depths and scalebar endpoint depths ---
            # Prefer depths provided by the database (marker_pair?_depth fields). If not available, fall back to model Z.
            try:
                # Build mapping from marker label -> depth (from DB)
                marker_depth_map = {}
                for idx, pair in enumerate((job.marker_pair1, job.marker_pair2, job.marker_pair3, job.marker_pair4), start=1):
                    depth_val = getattr(job, f'marker_pair{idx}_depth', None)
                    if pair and pair.lower() != 'na':
                        try:
                            t1, t2 = pair.split(',')
                            lbl1 = 'target ' + t1.strip()
                            lbl2 = 'target ' + t2.strip()
                            # store depth if present (could be None)
                            marker_depth_map[lbl1] = depth_val
                            marker_depth_map[lbl2] = depth_val
                        except Exception:
                            # ignore malformed pair strings
                            continue

                db.log_step(job.id, 4, "info", "Recording marker depths (DB if available, else model)")
                opf.write('\nMarker depths:\n')
                for marker in MarkersList:
                    depth_to_report = None
                    if marker and marker.label:
                        # Prefer DB depth
                        if marker.label in marker_depth_map and marker_depth_map[marker.label] is not None:
                            depth_to_report = marker_depth_map[marker.label]
                        else:
                            # Fallback: model z
                            try:
                                depth_to_report = float(marker.position[2]) if marker.position is not None else None
                            except Exception:
                                depth_to_report = None

                    msg = f"Marker {marker.label} depth: {depth_to_report if depth_to_report is not None else 'N/A'}"
                    print(msg)
                    logger.info(msg)
                    db.log_step(job.id, 4, "info", msg)
                    opf.write(msg + '\n')

                opf.write('\nScalebar endpoint depths:\n')
                for scalebar in chunk.scalebars:
                    # Try to infer scalebar endpoint labels; if endpoints are markers, use DB depths where available.
                    def _depth_for_endpoint(pt):
                        # If endpoint is a Camera, try camera.center z (fallback to N/A)
                        try:
                            if isinstance(pt, Metashape.Camera):
                                return float(pt.center[2]) if pt.center is not None else None
                            # If endpoint is a Marker-like object with label and position
                            label = getattr(pt, 'label', None)
                            if label and label in marker_depth_map and marker_depth_map[label] is not None:
                                return marker_depth_map[label]
                            if getattr(pt, 'position', None) is not None:
                                return float(pt.position[2])
                        except Exception:
                            return None
                        return None

                    d0 = _depth_for_endpoint(scalebar.point0)
                    d1 = _depth_for_endpoint(scalebar.point1)
                    msg = f"Scalebar {scalebar.label} depths: point0={d0 if d0 is not None else 'N/A'}, point1={d1 if d1 is not None else 'N/A'}"
                    print(msg)
                    logger.info(msg)
                    db.log_step(job.id, 4, "info", msg)
                    opf.write(msg + '\n')
            except Exception as e:
                logger.warning(f"Failed to record marker/scalebar depths: {e}")
                opf.write(f"Failed to record marker/scalebar depths: {e}\n")

            doc.save()
            db.log_step(job.id, 4, "completed", "Scaling complete")
            print("Step 4. Scaling Finished")
            logger.info("Step 4. Scaling Finished")
            db.set_end_step(job.id, 4)

            # Exit after completing Step 4 as requested
            logger.info("Exiting after Step 4 for Human Review")
            db.log_step(job.id, 4, "info", "Exiting process after Step 4")
            db.update_job_status(job.id, "ready for review")
            try:
                opf.close()
            except Exception:
                pass
            return  # This return is intentional - waiting for human review

        #==============================================================================================
        # Step 5: Error reduction part 2
        if (start_step ==4 and end_step ==4 and job.status =='reviewed') or (start_step ==5 and job.status =='failed'):
            db.log_step(job.id, 5, "running", "Starting error reduction phase 2")
            print("Step 5. Error reduction part 2")
            logger.info("Step 5. Error reduction part 2")
            opf.write('\nError Reduction\n{0}\n'.format(underline))
            db.set_start_step(job.id,5)
            
            points = chunk.tie_points.points
            keep_percent = 90
            
            # Filter Reprojection Error
            db.log_step(job.id, 5, "info", "Filtering reprojection error")
            f = Metashape.TiePoints.Filter()
            f.init(chunk, criterion = Metashape.TiePoints.Filter.ReprojectionError)
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
            
            # Final optimization
            db.log_step(job.id, 5, "info", "Performing final camera optimization")
            chunk.optimizeCameras()
            doc.save()
            
            db.log_step(job.id, 5, "completed", "Error reduction complete")
            print("Step 5. Error reduction part 2 Done")
            logger.info("Step 5. Error reduction part 2 Done")
            db.set_end_step(job.id, 5)

        #==============================================================================================
        # Step 6: Build Dense Cloud
        if (start_step == 5 and end_step ==5):
                # --- CLEANUP BLOCK --- #
            # Prevent reprocessing from accumulating previous models, DEMs, or orthos
            try:
                logger.info("Cleaning up old products before dense cloud build...")
                db.log_step(job.id, 6, "info", "Cleaning up old dense clouds, DEMs, and orthomosaics before rebuild")

                # Remove depth maps
                if chunk.depth_maps:
                    logger.info("Removing existing depth maps")
                    # Clear existing depth maps (method depends on Metashape version)
                    if hasattr(chunk, 'clearDepthMaps'):
                        chunk.clearDepthMaps()
                    elif hasattr(chunk, 'depth_maps'):
                        chunk.depth_maps.clear()
                    else:
                        logger.warning("No method found to clear depth maps for this Metashape version.")
                # Remove dense cloud
                if getattr(chunk, 'point_cloud', None):
                    logger.info("Clearing existing dense cloud")
                    try:
                        chunk.point_cloud.clear()  # preferred
                    except Exception:
                        chunk.point_cloud = None  # fallback

                # Remove models (meshes)
                if getattr(chunk, "models", None):
                    for model in list(chunk.models):
                        logger.info(f"Removing existing model: {model.label}")
                        chunk.remove(model)
                else:
                    logger.debug("No models to remove")

                # Remove DEMs
                if getattr(chunk, "elevation", None):
                    for dem in list(chunk.elevation):
                        logger.info(f"Removing existing DEM: {dem.label}")
                        chunk.remove(dem)
                else:
                    logger.debug("No DEMs to remove")

                # Remove orthomosaics
                if getattr(chunk, "orthomosaics", None):
                    for ortho in list(chunk.orthomosaics):
                        logger.info(f"Removing existing orthomosaic: {ortho.label}")
                        chunk.remove(ortho)
                else:
                    logger.debug("No orthomosaics to remove")

        # Save and reload project to release .files memory
                doc.save()
                doc.open(doc.path)
                doc.save()
                chunk = doc.chunk
                logger.info("Cleanup complete - ready to rebuild dense cloud.")
                db.log_step(job.id, 6, "info", "Cleanup complete - ready to rebuild dense cloud.")
            except Exception as e:
                logger.warning(f"Cleanup before dense cloud failed: {e}")
                db.log_step(job.id, 6, "warning", f"Cleanup before dense cloud failed: {e}")

            # --- Continue with your existing dense cloud build process --- #
            try:
                db.log_step(job.id, 6, "running", "Starting dense cloud construction")
                print("Step 6. Build Dense Cloud")
                logger.info("Step 6. Build Dense Cloud")
                db.set_start_step(job.id,6)
        
                # Build depth maps
                db.log_step(job.id, 6, "info", "Building depth maps")
                chunk.buildDepthMaps(downscale=4, filter_mode=Metashape.MildFiltering)

                # Build point cloud
                db.log_step(job.id, 6, "info", "Building point cloud")
                chunk.buildPointCloud(point_colors=True, point_confidence=True)
                
                # Filter low confidence points
                db.log_step(job.id, 6, "info", "Filtering low confidence points")
                chunk.point_cloud.setConfidenceFilter(0, 1)
                chunk.point_cloud.removePoints(list(range(128)))
                chunk.point_cloud.resetFilters()
            except Exception as e:
                logger.error(f"Error during dense cloud build: {e}")
                db.log_step(job.id, 6, "failed", f"Dense cloud build failed: {e}")
                raise
            doc.save()

            db.log_step(job.id, 6, "completed", "Dense cloud construction complete")
            print("Step 6. Build Dense Cloud complete")
            logger.info("Step 6. Build Dense Cloud complete")
            db.set_end_step(job.id, 6)
            
            # Exit after completing Step 6 as requested
            logger.info("Exiting after Step 6 for Human Review")
            db.log_step(job.id, 6, "info", "Exiting process after Step 6")
            db.update_job_status(job.id, "ready for review")
            try:
                opf.close()
            except Exception:
                pass
            return  # This return is intentional - waiting for human review

        #==============================================================================================
        # Step 7: Build and export DEM and Orthomosaic
        if (start_step ==6 and end_step ==6 and job.status == 'reviewed') or (start_step == 7 and job.status =='failed'):

            # --- Pre-check dense cloud existence ---
         #   if not getattr(chunk, 'point_cloud', None) or len(chunk.point_cloud.points) == 0:
          #      msg = "Step 7 skipped: Dense cloud missing. Cannot build DEM/orthomosaic."
          #      logger.warning(msg)
           #     db.log_step(job.id, 7, "warning", msg)
          #      db.update_job_status(job.id, "failed", msg)
            #    return  # exit Step 7 gracefully

            # db.log_step(job.id, 7, "info", "Cleaning DEMs and orthomosaics before rebuild")
          #   cleanup_rebuild_products(doc, chunk, logger)
                    
            db.log_step(job.id, 7, "running", "Starting DEM and orthomosaic generation")
            print("Step 7. Build and export DEM and Orthomosaic")
            logger.info("Step 7. Build and export DEM and Orthomosaic")
            db.set_start_step(job.id, 7)

            
            # Create ARC directory
            arc_path = os.path.join(prod_path, 'ARC')
            try:
                os.mkdir(arc_path)
            except:
                db.log_step(job.id, 7, "info", "ARC directory already exists")
                
            # Check and log region size before processing
       #     region_size = chunk.region.size
        #    logger.info(f"Region size before processing: {region_size}")
        #    db.log_step(job.id, 7, "info", f"Initial region size: {region_size}")
            
            # Reset region if it seems too large (optional)
#max_expected_size = 100  # meters - adjust this based on your typical project size
        #     if any(size > max_expected_size for size in region_size):
         #       logger.warning(f"Region size unusually large: {region_size}")
          #      db.log_step(job.id, 7, "warning", f"Unusually large region detected: {region_size}")
                
            # Re-verify orientation if not a spiral survey
         #   if not job.site_id.startswith("OCC"):
          #      logger.info("Re-checking orientation before export")
          #      db.log_step(job.id, 7, "info", "Re-verifying orientation")
         #       orient_model(chunk, logger)
            
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
                thumb_path = os.path.join(prod_path, folder_name + "_thumb.jpg")
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
            pt_file = os.path.join(prod_path, f'{folder_name}.ply')
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

# def main():
 #   """Main function to process jobs from database in priority order."""
 #   try:
  #      # Get pending jobs ordered by priority
 #       jobs = db.get_pending_jobs()
  #      
 #       if not jobs:
  #          logger.info("No pending jobs found.")
 #           return
            
   #     for job in jobs:
 #           logger.info(f"Processing job {job.id} (priority: {getattr(job, 'priority', 'N/A')}): {getattr(job, 'site_id', '')}")
 #           MetashapeProcess(job)
            
 #   except Exception as e:
 #       logger.error(f"Error in main processing loop: {str(e)}")
 #       raise


def main():
    db = MyDatabaseManager()
    db.connect()
    try:
        while True:
            jobs = db.get_pending_jobs()
            if not jobs:
                logger.info("No more pending jobs found. Exiting.")
                break
                
            for job in jobs:
                try:
                    logger.info(f"Starting processing for job {job.id} ({job.site_id})")
                    MetashapeProcess(job)
                except Exception as job_error:
                    error_msg = f"Error processing job {job.id} ({job.site_id}): {str(job_error)}"
                    logger.error(error_msg)
                    
                    # Update job status in database
                    try:
                        db.update_job_status(job.id, "failed", error_msg)
                        db.log_step(job.id, -1, "error", error_msg)
                        logger.info(f"Job {job.id} marked as failed, continuing with next job")
                    except Exception as db_error:
                        logger.error(f"Failed to update database for failed job {job.id}: {str(db_error)}")
                    
                    # Clear Metashape document to free memory
                    try:
                        doc = Metashape.app.document
                        doc.clear()
                    except:
                        pass
                    
                    continue  # Move on to next job
                    
            # Optional pause between job batches
            import time
            time.sleep(5)  # 5 second pause between batches
            
    except Exception as e:
        logger.exception("Fatal error in main processing loop")
        raise
    finally:
        try:
            db.close()
        except:
            pass


if __name__ == "__main__":
    main()
