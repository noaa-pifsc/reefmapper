# Auto batch process for Agisoft Metashape
# Following Structure-from-Motion workflow
# M. Akridge, 2024/07/11
# F. Lichowski, 2024/07/10
#
################################################################################################
import os
import sys
import csv
import json
import math
import Metashape
import shutil
# path where html2text folder resides
# sys.path.insert(0, 'M:\\FL_test\\Scripts\\') 
# pip install html2text
import html2text
              
################################################################################################  
# set path and filename of processing log (as text string)
process_log = 'N:\StRS_Sites\2024\MP2404_MHI\OAH\StRS_Sites_2024_ProcessingLog_V1.csv'
batch_no = 1 # set number of batch to be processed (1-n)

################################################################################################  
# function to make list of all photos in root_path
def getPhotoList(root_path, photoList):
    files = [f for f in os.listdir(root_path) if os.path.isfile(os.path.join(root_path, f))]
    for photo in files:
        if ("jpg" or "jpeg" or "JPG" or "JPEG") in photo.lower(): 
            add_path = os.path.join(root_path, photo)
            photoList.append(add_path)           

# function to detect marker
def find_marker(label, chunk):
	for marker in chunk.markers:
		if label == marker.label:
			return marker
	return None

# main function for automatic batch processing
def MetashapeProcess(root_path, folder_name, start_step=1, end_step=7, 
                     valid_markers=[], quality=0.5, survey_year='2024'):
                    
                    # root_path = folder name containing of subfolders with different sites
                    # folder_name = subfolder with individual site
                    # start_step = number of step at which the script starts: 
                    #              1=create new .psx file and .files folder
                    #              2=adds photos, match photos, and align cameras
                    #              3=sparse point cloud filtering
                    #              4=scaling, detect markers and add scalebar
                    #              5=error reduction, part 2
                    #              6=build dense cloud
                    #              7=build and export DEM and Orthomosaic, generate report 
                    # end_step = step to which the model runs, including that step
                    # valid_markers= list of marker pair numbers
                    # quality = image quality threshold
                    # sur_year = survey year for export
                    
    #==============================================================================================
    ### 1) Initialize			
    # Metashape.app.ConsolePane.clear() # not woking in 1.7.4
    underline = 50*'-'
    
    # construct the document class (with no chunks)
    doc = Metashape.app.document
    doc.clear()
    
     # create new Products folder for project files
    prod_path = os.path.join(root_path, 'Products')
    psxfile = os.path.join(prod_path, folder_name+'.psx')
    try:
        os.mkdir(prod_path)
    except:
        print(prod_path, 'directory already exists')    
    
    # write log and readme text file with processing information
    Metashape.app.settings.log_enable = True
    log_file = os.path.join(prod_path, folder_name+'_log.txt')
    Metashape.app.settings.log_path = log_file
    
    if start_step==1:
        print('Step 1 new')
        # start new readme file
        opf = open(os.path.join(prod_path, folder_name+'_readme.txt'), 'w')       

        # remove .psx file, if it exists
        if os.path.exists(psxfile):
            os.remove(psxfile)
            print('Deleted', psxfile)
            
        # remove log file, if it exists
        if os.path.exists(log_file):          
            os.remove(log_file)
            print('Deleted', log_file)
        
        # remove .files folder, if it exists
        file_path = os.path.join(prod_path, folder_name+'.files')
        if os.path.exists(file_path):
            shutil.rmtree(file_path, ignore_errors=True)
            print('Deleted', file_path)
                       
        # save .psx file and add new chunk     
        doc.save(psxfile)
        chunk = doc.addChunk()
        print('Saved project to: ' + psxfile)
    
    else:
        print('Step 1 existing')
        # for start_step>1: open existing readme and .psx file and load chunk
        opf = open(os.path.join(prod_path, folder_name+'_readme.txt'), 'a')
        doc.open(psxfile)
        chunk = doc.chunk
    
    chunk.camera_location_accuracy = Metashape.Vector((0.1, 0.1, 0.15))
    opf.write('Readme file for {0}\n{1}\n'.format(folder_name, underline)) 
        
    #==============================================================================================
    ### 2) Add and align photos   
    if start_step<=2 and end_step>=2:
        print('Step 2')
        opf.write('\nAlign photos\n{0}\n'.format(underline))        
        
        ## get photo list
        photoList = []
        extra = ''
        getPhotoList(os.path.join(root_path, extra), photoList)
        n = len(photoList)
       
        ## add photos
        print('Adding '+str(n)+' photos')
        opf.write('Adding '+str(n)+' photos\n')
        chunk.addPhotos(photoList) 

        ## estimate image quality and disable those with poor quality
        chunk.analyzePhotos(chunk.cameras)
        bad_quality = 0
        quality_log = [0]*n # 0=good, 1=poor
        for qc, camera in enumerate(chunk.cameras):
            if float(camera.meta['Image/Quality'])<quality and quality_log[qc-1]==0:
                bad_quality += 1
                quality_log[qc] = 1
                camera.enabled = False
                print(qc, camera.label, float(camera.meta['Image/Quality'])) 
                opf.write('Image {0} with quality {1:.3f} disabled\n'.format(
                                camera.label, float(camera.meta['Image/Quality'])))
        n_enabled = n-bad_quality
        print(str(bad_quality)+' photo(s) ('+str(round(bad_quality/n*100,1))+ \
              '%) of quality below '+str(quality))                 
        opf.write(str(bad_quality)+' photo(s) ('+str(round(bad_quality/n*100,1))+ \
              '%) of quality below '+str(quality)+'\n')         
        print(quality_log)

        ## perform image matching for the chunk frame
        chunk.matchPhotos(generic_preselection=True, reference_preselection=False, 
                          filter_mask=False, keypoint_limit=40000, tiepoint_limit=0)
    
        ## perform photo alignment for the chunk
        chunk.alignCameras(adaptive_fitting=False) # disable adaptive fitting of distortion coefficients
         
        ## ensure that all photos are aligned
        thresh_align = 15
        counter = 0 
        for camera in chunk.cameras:
            if camera.transform: 
                counter += 1
        print('Enabled images: {0} '.format(n_enabled))                
        print('Aligned images: {0} '.format(counter))
        opf.write('Enabled images: {0}\nAligned images: {1}\n'.format(n_enabled, counter)) 
        if n_enabled-counter==0:
            print('All enabled images aligned!')
        elif n_enabled-counter>0 and n_enabled-counter<=thresh_align:
            print(n_enabled-counter, 'enabled images did not align!')       
        elif n_enabled-counter>thresh_align:
            print('Try again or reevaluate the quality of corresponding photos') 
            quit()
             
        # save an original copy of the SPC in case you need to undo any point filtering 
        doc.save()          
        shutil.copyfile(psxfile, os.path.join(prod_path, folder_name+'_bkup.psx'))
        if os.path.exists(os.path.join(prod_path, folder_name+'_bkup.files')):
            shutil.rmtree(os.path.join(prod_path, folder_name+'_bkup.files'), ignore_errors=True)
        shutil.copytree(os.path.join(prod_path, folder_name+'.files'), 
                        os.path.join(prod_path, folder_name+'_bkup.files'))
             
    #==============================================================================================  
    ### 3) Sparse point cloud filtering
    if start_step<=3 and end_step>=3:
        print('Step 3')
        opf.write('\nSparse point cloud filtering\n{0}\n'.format(underline))
        
        keep_percent = 51
        points = chunk.point_cloud.points
        
        # Optimize 1
        chunk.optimizeCameras(tiepoint_covariance=True)
        
        # Filter: Reconstruction Uncertainty
        f = Metashape.PointCloud.Filter()
        f.init(chunk, criterion = Metashape.PointCloud.Filter.ReconstructionUncertainty)
        list_values = f.values
        list_values_valid = list()
        for i in range(len(list_values)):
            if points[i].valid:
                list_values_valid.append(list_values[i])
        list_values_valid.sort()
        target = int(len(list_values_valid) * keep_percent / 100)
        RecUncert = list_values_valid[target]
        print('Reconstruction Uncertainty threshold to keep 51%:', RecUncert)
        if RecUncert < 10:
            RecUncert = 10
            print('Reconstruction Uncertainty threshold set to', RecUncert)
        opf.write('Reconstruction Uncertainty threshold to keep 51%: {0:.2f}\n'.format(RecUncert))
        f.removePoints(RecUncert) 
        
        # Optimize 2
        chunk.optimizeCameras(tiepoint_covariance=True)

        # Filter Projection Accuracy     
        f = Metashape.PointCloud.Filter()
        f.init(chunk, criterion = Metashape.PointCloud.Filter.ProjectionAccuracy)
        list_values = f.values
        list_values_valid = list()
        for i in range(len(list_values)):
            if points[i].valid:
                list_values_valid.append(list_values[i])
        list_values_valid.sort()
        target = int(len(list_values_valid) * keep_percent / 100)
        ProjAcc = list_values_valid[target]
        print('Projection Accuracy threshold to keep 51%:', ProjAcc)
        if ProjAcc < 2:
            ProjAcc = 2
            print('Projection Accuracy threshold set to', ProjAcc)
        opf.write('Projection Accuracy threshold to keep 51%: {0:.2f}\n'.format(ProjAcc))
        f.removePoints(ProjAcc)     
        
        # Optimize 3
        chunk.optimizeCameras(tiepoint_covariance=True)
        doc.save()
    
    #==============================================================================================
    ### 4) Scaling
    if start_step<=4 and end_step>=4:
        
        print('Step 4')
        opf.write('\nScaling\n{0}\n'.format(underline))
                  
        # detect markers   
        chunk.detectMarkers()
        error_thresh = 0.4  # marker error threshold      
        
        # reduce marker list to valid markers     
        MarkersList = list(chunk.markers)
        print(MarkersList)
        for marker in MarkersList:
            print(underline, '\n', marker.label)
            if not marker:
                print(marker.label + " not found, skipping...")
                MarkersList.remove(marker)
                chunk.remove(marker)
                continue
            if not marker.position:
                print(marker.label + " is not defined in 3D, skipping...")
                MarkersList.remove(marker)
                chunk.remove(marker)
                continue            
            if marker.label not in valid_markers:
                print(marker.label + " not a valid marker, skipping...")
                MarkersList.remove(marker)
                chunk.remove(marker)
                continue
            
            # make list of error for each image
            opf.write('Marker: {0}\n'.format(marker.label))
            pix_error = error_thresh
            while pix_error >= error_thresh:
                total = (0, 0)
                cam_err = []                
                for camera in marker.projections.keys():
                    if not camera.transform:
                        continue
                    proj = marker.projections[camera].coord
                    reproj = camera.project(marker.position)
                    error = (proj - reproj).norm()
                    total = (total[0] + error**2, total[1] + 1)
                    # print(total)
                    cam_err.append((camera.label, error))
                
                # calculate total error for marker
                pix_error = math.sqrt(total[0] / total [1]) 
                # print('intermediate error', pix_error)
                
                # remove images with highest error as long as total error is >= threshold                 
                if pix_error >= error_thresh:
                    max_err = sorted(cam_err, key=lambda x: x[1], reverse=True)[0]
                    # print (sorted(cam_err, key=lambda x: x[1], reverse=True))
                    print('Removed:', max_err[0], max_err[1])
                    opf.write('Removed {0} with pix error {1:.4}\n'.format(
                              max_err[0], max_err[1]))
                    for camera in marker.projections.keys():
                        if camera.label == max_err[0]:
                            marker.projections[camera] = None

            # write final error to text file
            opf.write('Marker: {1}, projections: {2}, total error {3:.4f} pix\n{0}\n\n'.format(
                              underline, marker.label, len(cam_err), pix_error))
            print('Marker:', marker.label, len(cam_err), pix_error)
        
        # check validity of marker (thanks GUA-2662)
        print(MarkersList)
        for marker in MarkersList:
            print(marker.label)
            if marker.label not in valid_markers:
                print(marker.label + " not a valid marker, skipping 2...")
                MarkersList.remove(marker)
                chunk.remove(marker)
                continue
        print(MarkersList)

        # double check validity of marker (thanks GUA-2662)
        for marker in MarkersList:
            print(marker.label)
            if marker.label not in valid_markers:
                print(marker.label + " not a valid marker, skipping 3...")
                MarkersList.remove(marker)
                chunk.remove(marker)
                continue
        print(MarkersList)

        # add scalebar
        sb_dist = 0.25 # distance between the two markers in meters
        for i, marker in enumerate(MarkersList):
            if i % 2 > 0:
                sb = chunk.addScalebar(chunk.markers[i-1], chunk.markers[i])
                sb.reference.distance = sb_dist
                
        # update transform to scale the model based on known ground control distances        
        chunk.updateTransform() 
        
        # Check the error of scalebars
        for scalebar in chunk.scalebars:
            print(scalebar)
            dist_source = scalebar.reference.distance
            if not dist_source:
                continue # skipping scalebars without source values
            if type(scalebar.point0) == Metashape.Camera:
                if not (scalebar.point0.center and scalebar.point1.center):
                    continue # skipping scalebars with undefined ends
                dist_estimated = (scalebar.point0.center - scalebar.point1.center).norm() * chunk.transform.scale
            else:
                if not (scalebar.point0.position and scalebar.point1.position):
                    continue # skipping scalebars with undefined ends
                dist_estimated = (scalebar.point0.position - scalebar.point1.position).norm() * chunk.transform.scale            
            dist_error = dist_estimated - dist_source
            print(scalebar.label, str(dist_source), str(dist_estimated), str(dist_error))
            opf.write('Scalebar: {0}, distance: {1:.2f}, esimated distance: {2:.5f}, error: {3:.6f}\n'.format(
                              scalebar.label, dist_source, dist_estimated, dist_error))
            if dist_error>=0.002:
                print('Scalebar error too high: lower error_thresh and run again starting at step 4!')
                opf.write('Scalebar error too high: lower error_thresh and run again starting at step 4!\n')
            scalebar.reference.enabled = False
        doc.save()
        
    #==============================================================================================
    ### 5) Error reduction, part 2
    if start_step<=5 and end_step>=5:
        print('Step 5')
        opf.write('\nError Reduction\n{0}\n'.format(underline))
        points = chunk.point_cloud.points
    
        ## Filter Reprojection Error
        keep_percent = 90
        f = Metashape.PointCloud.Filter()
        f.init(chunk, criterion = Metashape.PointCloud.Filter.ReprojectionError)
        list_values = f.values
        list_values_valid = list()
        for i in range(len(list_values)):
            if points[i].valid:
                list_values_valid.append(list_values[i])
        list_values_valid.sort()
        target = int(len(list_values_valid) * keep_percent / 100)
        Rep_Err = list_values_valid[target]
        print('Reprojection Error threshold to keep 90%:', Rep_Err)
        opf.write('Reprojection Error threshold to keep 90%: {0:.2f}\n'.format(Rep_Err))
        f.removePoints(Rep_Err)        
        
        ## Optimize
        chunk.optimizeCameras() # Optimization will be based on estimated camera poses
        doc.save()

    #==============================================================================================
    ### 6) Build Dense Cloud
    if start_step<=6 and end_step>=6:
        print('Step 6')
        # buildDepthMaps(downscale=4, filter_mode=MildFiltering[, cameras ], reuse_depth=False,
        # max_neighbors=40, subdivide_task=True, workitem_size_cameras=20,
        # max_workgroup_size=100[, progress ])
        chunk.buildDepthMaps(downscale=4, filter_mode=Metashape.MildFiltering)
        
        # Generate dense cloud for the chunk.
        # buildDenseCloud(point_colors=True, point_confidence=False, keep_depth=True,
        # max_neighbors=100, subdivide_task=True, workitem_size_cameras=20, 
        # max_workgroup_size=100[, progress ])
        chunk.buildDenseCloud(point_colors=True, point_confidence=True)
        
        # select points to remove (0-1 confidence)
        chunk.dense_cloud.setConfidenceFilter(0, 1)
        chunk.dense_cloud.removePoints(list(range(128))) # removes all "visible" points of dense cloud
        chunk.dense_cloud.resetFilters()
        
        # Add depth and rectify model
        # DONE MANUALLY (for now, could easily be read from the processing log)
        
        doc.save()
    #==============================================================================================
    ### 7) Build and export DEM and Orthomosaic, generate report
    if start_step==7:
        print('Step 7')
        
        # make ARC directory
        arc_path = os.path.join(prod_path, 'ARC')
        try:     
            os.mkdir(arc_path)
        except:
            pass 
                                 
        # build DEM using all defaults
        chunk.buildDem(source_data=Metashape.DenseCloudData)
        
        # build Orthomosaic  
        chunk.buildOrthomosaic(surface_data=Metashape.ElevationData, fill_holes=True, 
                                ghosting_filter=False, refine_seamlines=False, resolution=0.0005)
        
        # set export projection
        out_projection = Metashape.OrthoProjection()
        out_projection.crs = Metashape.CoordinateSystem("EPSG::32604")
        
        # export DEM
        chunk.exportRaster(os.path.join(arc_path, survey_year+'_'+folder_name+'_dem.tif'), 
                           projection = out_projection, resolution=0.001, save_world=True, 
                           source_data=Metashape.ElevationData)
        
        # export DEM for Hannah
        chunk.exportRaster(os.path.join(arc_path, survey_year+'_'+folder_name+'_dem_1cm.tif'), 
                           projection = out_projection, resolution=0.01,
                           save_world=True, source_data=Metashape.ElevationData)    
            
        # compression parameters
        compression = Metashape.ImageCompression()
        compression.tiff_compression = Metashape.ImageCompression.TiffCompressionLZW
        compression.jpeg_quality = 99
        compression.tiff_big = False
        compression.tiff_tiled = False
        compression.tiff_overviews = False
        
        # export Orthomosaic         
        chunk.exportRaster(os.path.join(arc_path, survey_year+'_'+folder_name+'_mos.tif'), 
                            projection = out_projection, resolution=0.0005, 
                            image_compression=compression, save_world=True,
                            save_alpha=False, source_data=Metashape.OrthomosaicData)
        
        # export report
        url = os.path.join(arc_path, survey_year+'_'+folder_name+'_rpt.html')
        chunk.exportReport(url, title=survey_year+'_'+folder_name, 
                            description="Processing Report")
        
        # read text from html report
        with open(url) as f:
            html = f.read()
        h = html2text.HTML2Text()
        t = h.handle(html)
        lines = t.splitlines()
        
        # save text from html report as csv file
        with open(os.path.join(arc_path, url[:-5]+'.csv'), 'w', newline='') as file:
            writer = csv.writer(file)
            for line in lines:
                words = line.strip().split('|')
                writer.writerow(words)

        # extract meta and cams
        cams = chunk.cameras
        proj_path = doc.path
        proj_dir, proj_name = os.path.split(proj_path)
        proj_name = proj_name[:-4]
        outputs = {}
        cams_filename = proj_dir + '/' + proj_name + '.cams.xml'
        meta_filename = proj_dir + '/' + proj_name + '.meta.json'
        meta_file = open(meta_filename, 'w')
        chunk.exportCameras(cams_filename)
        
        for cam in cams:
            key = cam.key
            path = cam.photo.path
            center = cam.center
            if center is not None:
                center = list(center)
            agi_trans = cam.transform
            trans = None
            if agi_trans is not None:
                trans = [list(agi_trans.row(n)) for n in range(agi_trans.size[1])]
            outputs[key] = {'path' : path, 'center' : center, 'transform' : trans}
            
        print(outputs)
        meta_file.write(json.dumps({'cameras' : outputs}, indent=4))
        meta_file.close()
        
        # export point cloud
        pt_file = os.path.join(prod_path, folder_name+'.ply')
        if os.path.isfile(pt_file)==False:
            chunk.exportPoints(pt_file, source_data = Metashape.DenseCloudData)
       
    #==============================================================================================
    # Save the project 
    doc.save()
    opf.close()
    print('Finished processing: ', folder_name)

###################################################################################################   
# open processing log and execute MetashapeProcess function
with open(os.path.join(r''+process_log), newline='') as csvfile:  
    ml = [] # empty marker list
    spamreader = csv.reader(csvfile, delimiter=',')
    for col in spamreader:
        if col[2]==str(batch_no): # number of batch to be processed
            fpath = r''+col[0] # read project filepath
            print(fpath)
            if os.path.exists(fpath):
                quality_thres = float(col[9]) # image quality threshold 
                sur_year = str(col[10]) # survey year for export
                for r in range(5, 9): # make list of valid markers from column F-I
                    if col[r] != 'NA':
                        t1, t2 = col[r].split(',')
                        ml.append('target '+t1.strip())
                        ml.append('target '+t2.strip())
                # run metashape function with values from processing log
                MetashapeProcess(fpath, col[1], start_step=int(col[3]), end_step=int(col[4]), 
                                 valid_markers=ml, quality=quality_thres, survey_year=sur_year)
            else:
                print(fpath, 'does not exist, move to next project!')
