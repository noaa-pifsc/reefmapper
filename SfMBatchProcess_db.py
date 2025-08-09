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
import shutil
import logging
import html2text
from datetime import datetime
from database.db_manager import MyDatabaseManager, ProcessingJob

################################################################################################  
# Setup logging to both file and console
LOG_FORMAT = '%(asctime)s [%(levelname)s] %(message)s'
logger = logging.getLogger("reefmapper")
logger.setLevel(logging.INFO)
formatter = logging.Formatter(LOG_FORMAT)

# Log to file (overwrites each run)
# file_handler = logging.FileHandler('reefmapper_batch.log', mode='w') #issues with permission
file_handler = logging.FileHandler(r'C:\Users\PICHLMRUser\Desktop\reefmapper\reefmapper_batch.log', mode='w') 
file_handler.setFormatter(formatter)
logger.addHandler(file_handler)

# Log to console
console_handler = logging.StreamHandler()
console_handler.setFormatter(formatter)
logger.addHandler(console_handler)

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
        
        logger.info(f"Starting batch: {folder_name} (steps {start_step}-{end_step})")
        db.log_step(job.id, 0, "started", f"Starting batch processing for {folder_name}")
        underline = 50*'-'

        #==============================================================================================
        # Step 1: Initialize
        db.log_step(job.id, 1, "running", "Initializing project")
        print("Step 1. Initialize running")
        logger.info("Step 1. Initialize running")
        db.update_job_status(job.id,'running')


        # construct the document class (with no chunks)
        doc = Metashape.app.document
        doc.clear()
        
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
        
        # write log and readme text file with processing information
        Metashape.app.settings.log_enable = True
        log_file = os.path.join(prod_path, folder_name+'_log.txt')
        Metashape.app.settings.log_path = log_file
        
        if start_step==1:
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
        if start_step ==1 and end_step ==1:
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

                ## add photos
                db.log_step(job.id, 2, "info", f"Adding {n} photos")
                print(f'Adding {n} photos')
                logger.info(f'Adding {n} photos')
                opf.write(f'Adding {n} photos\n')
                chunk.addPhotos(photoList)
                doc.save(psxfile)

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
                db.log_step(job.id, 2, "error", msg)
                opf.write(msg + '\n')
                return  # ✅ Prevent continuing to later steps


        #==============================================================================================
        # Step 3: Sparse point cloud filtering
        if start_step==2 and end_step==2:
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
        if start_step==3 and end_step==3:
            db.log_step(job.id, 4, "running", "Starting model scaling")
            print("Step 4. Scaling")
            print("Metashape version:", Metashape.app.version)
            db.log_step(job.id, 4, "running", f"Metashape version: {Metashape.app.version}")
            logger.info("Step 4. Scaling")
            opf.write('\nScaling\n{0}\n'.format(underline))
            db.set_start_step(job_id,4)
            
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
                                
            # Add scalebars
            db.log_step(job.id, 4, "info", "Adding scalebars")
            sb_dist = 0.25  # distance between markers in meters
            for i, marker in enumerate(MarkersList):
                if i % 2 > 0:
                    sb = chunk.addScalebar(chunk.markers[i-1], chunk.markers[i])
                    sb.reference.distance = sb_dist
            
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
                
                scalebar.reference.enabled = False
            
            doc.save()
            db.log_step(job.id, 4, "completed", "Scaling complete")
            print("Step 4. Scaling Finished")
            logger.info("Step 4. Scaling Finished")
            db.set_end_step(job.id, 4)

        #==============================================================================================
        # Step 5: Error reduction part 2
        if (start_step==4 and end_step==4) or (start_step ==4 and job.status =='failed'):
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
        if start_step==5 and end_step==5:
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
            
            doc.save()
            db.log_step(job.id, 6, "completed", "Dense cloud construction complete")
            print("Step 6. Build Dense Cloud complete")
            logger.info("Step 6. Build Dense Cloud complete")
            db.set_end_step(job.id, 6)

        #==============================================================================================
        # Step 7: Build and export DEM and Orthomosaic
        if start_step==6 and end_step ==6:
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
            
            # Build DEM
            db.log_step(job.id, 7, "info", "Building DEM")
            chunk.buildDem(source_data=Metashape.PointCloudData)
            
            # Build Orthomosaic
            db.log_step(job.id, 7, "info", "Building orthomosaic")
            chunk.buildOrthomosaic(surface_data=Metashape.ElevationData,
                                 fill_holes=True,
                                 ghosting_filter=False,
                                 refine_seamlines=False,
                                 resolution=0.0005)
            
            # Export DEM
            dem_path = os.path.join(arc_path, f'{survey_year}_{folder_name}_dem.tif')
            db.log_step(job.id, 7, "info", f"Exporting DEM to {dem_path}")
            chunk.exportRaster(dem_path,
                             resolution=0.001,
                             save_world=True,
                             source_data=Metashape.ElevationData)
            
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
            chunk.exportRaster(ortho_path,
                             resolution=0.0005,
                             image_compression=compression,
                             save_world=True,
                             save_alpha=False,
                             source_data=Metashape.OrthomosaicData)
            
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

        
        # On successful completion
        logger.info(f"Completed processing {folder_name}")
        db.update_job_status(job.id, "completed")

        
    except Exception as e:
        error_msg = f"Error processing {folder_name}: {str(e)}"
        logger.error(error_msg)
        db.log_step(job.id, -1, "error", error_msg)
        db.update_job_status(job.id, "failed", error_msg)
        raise

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
                MetashapeProcess(job)
    except Exception as e:
        logger.exception("Fatal error during Metashape processing.")
    finally:
        db.close()


if __name__ == "__main__":
    main()
