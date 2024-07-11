# Auto batch process for Agisoft Metashape
# Following Structure-from-Motion workflow
# M. Akridge, 2024/07/11
# F. Lichowski, 2024/07/10
################################################################################################
import os
import sys
import csv
import json
import logging
import math
import Metashape
import shutil
import html2text
# path where html2text folder resides
# sys.path.insert(0, 'M:\\FL_test\\Scripts\\') 
# pip install html2text
################################################################################################
# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)
################################################################################################
     
def get_photo_list(path, photo_list):
    """
    Recursively add photo paths from the given directory to the photo list.
    
    Args:
        path (str): Directory path to search for photos.
        photo_list (list): List to append photo paths to.
    """
    for root, _, files in os.walk(path):
        for file in files:
            if file.endswith(('jpg', 'jpeg', 'tif', 'tiff')):
                photo_list.append(os.path.join(root, file))

def add_photos(chunk, photo_list, quality):
    """
    Add photos to the chunk and disable poor quality images.
    
    Args:
        chunk: Metashape chunk object.
        photo_list (list): List of photo paths.
        quality (float): Quality threshold for disabling photos.
    """
    n = len(photo_list)
    logger.info('Adding %d photos', n)
    chunk.addPhotos(photo_list)
    chunk.analyzePhotos(chunk.cameras)

    bad_quality = 0
    quality_log = [0] * n

    for qc, camera in enumerate(chunk.cameras):
        if float(camera.meta['Image/Quality']) < quality and quality_log[qc-1] == 0:
            bad_quality += 1
            quality_log[qc] = 1
            camera.enabled = False
            logger.info('Image %s with quality %.3f disabled', camera.label, float(camera.meta['Image/Quality']))

    n_enabled = n - bad_quality
    logger.info('%d photo(s) (%.1f%%) of quality below %.2f', bad_quality, (bad_quality / n) * 100, quality)

def match_and_align_photos(chunk):
    """
    Perform image matching and alignment for the chunk.
    
    Args:
        chunk: Metashape chunk object.
    """
    chunk.matchPhotos(generic_preselection=True, reference_preselection=False, 
                      filter_mask=False, keypoint_limit=40000, tiepoint_limit=0)
    chunk.alignCameras(adaptive_fitting=False)

    thresh_align = 15
    counter = 0
    for camera in chunk.cameras:
        if camera.transform:
            counter += 1

    logger.info('Aligned %d of %d images', counter, len(chunk.cameras))

def save_project(doc, path):
    """
    Save the Metashape project.
    
    Args:
        doc: Metashape document object.
        path (str): Path to save the project.
    """
    doc.save(path)
    logger.info('Project saved at %s', path)

def process_chunk(chunk, prod_path, folder_name, survey_year):
    """
    Process the chunk and export necessary files.
    
    Args:
        chunk: Metashape chunk object.
        prod_path (str): Path to the production directory.
        folder_name (str): Name of the folder.
        survey_year (str): Survey year for export.
    """
    meta_file_path = os.path.join(prod_path, folder_name + '_meta.json')
    with open(meta_file_path, 'w') as meta_file:
        outputs = {}
        for key, value in chunk.cameras.items():
            path = value.photo.path
            center = value.center
            trans = value.transform
            outputs[key] = {'path': path, 'center': center, 'transform': trans}
        
        json.dump({'cameras': outputs}, meta_file, indent=4)
        logger.info('Metadata saved to %s', meta_file_path)

    pt_file = os.path.join(prod_path, folder_name + '.ply')
    if not os.path.isfile(pt_file):
        chunk.exportPoints(pt_file, source_data=Metashape.DenseCloudData)
        logger.info('Point cloud exported to %s', pt_file)

def MetashapeProcess(root_path, folder_name, start_step, end_step, valid_markers, quality, survey_year):
    """
    Main processing function for Metashape.
    
    Args:
        root_path (str): Root directory path.
        folder_name (str): Name of the folder.
        start_step (int): Starting step for the process.
        end_step (int): Ending step for the process.
        valid_markers (list): List of valid markers.
        quality (float): Quality threshold for images.
        survey_year (str): Survey year for export.
    """
    doc = Metashape.Document()
    doc.open(root_path)
    chunk = doc.chunk

    photo_list = []
    get_photo_list(os.path.join(root_path, ''), photo_list)
    add_photos(chunk, photo_list, quality)
    match_and_align_photos(chunk)
    process_chunk(chunk, root_path, folder_name, survey_year)
    save_project(doc, root_path)
    logger.info('Finished processing: %s', folder_name)

def main(process_log, batch_no):
    """
    Main function to read the processing log and execute Metashape processing.
    
    Args:
        process_log (str): Path to the processing log file.
        batch_no (int): Batch number to be processed.
    """
    with open(process_log, newline='') as csvfile:
        spamreader = csv.reader(csvfile, delimiter=',')
        for col in spamreader:
            if col[2] == str(batch_no):
                fpath = col[0]
                logger.info('Processing file: %s', fpath)
                if os.path.exists(fpath):
                    quality_thres = float(col[9])
                    survey_year = col[10]
                    valid_markers = [f'target {marker.strip()}' for marker in col[5:9] if marker != 'NA']

                    MetashapeProcess(fpath, col[1], start_step=int(col[3]), end_step=int(col[4]), 
                                     valid_markers=valid_markers, quality=quality_thres, survey_year=survey_year)
                else:
                    logger.error('%s does not exist, skipping to next project!', fpath)

if __name__ == '__main__':
    # set path and filename of processing log (as text string)
    process_log = 'N:\StRS_Sites\2024\MP2404_MHI\OAH\StRS_Sites_2024_ProcessingLog_V1.csv'  # Update with the actual path
    batch_no = 1 # Update with the actual batch number | set number of batch to be processed (1-n)
    main(process_log, batch_no)
