# https://agisoft.freshdesk.com/support/solutions/articles/31000148930-how-to-install-metashape-stand-alone-python-module
#import Metashape
import os
import csv
import json
import math
import shutil
from gooey import Gooey, GooeyParser
import signal
import html2text

# Mock Metashape module for testing
class MockMetashape:
    class app:
        class settings:
            log_enable = False
            log_path = ""
        document = None

    class Document:
        def __init__(self):
            self.chunk = MockMetashape.Chunk()

        def clear(self):
            print("Document cleared")

        def save(self, path):
            print(f"Document saved at {path}")

        def open(self, path):
            print(f"Document opened from {path}")

        def addChunk(self):
            return MockMetashape.Chunk()

    class Chunk:
        def __init__(self):
            self.cameras = []
            self.point_cloud = MockMetashape.PointCloud()
            self.transform = MockMetashape.Transform()
            self.markers = []
            self.scalebars = []

        def addPhotos(self, photo_list):
            print(f"Added {len(photo_list)} photos")

        def analyzePhotos(self, cameras):
            print("Analyzed photos")

        def matchPhotos(self, **kwargs):
            print("Matched photos")

        def alignCameras(self, **kwargs):
            print("Aligned cameras")

        def buildDepthMaps(self, **kwargs):
            print("Built depth maps")

        def buildDenseCloud(self, **kwargs):
            print("Built dense cloud")

        def buildDem(self, **kwargs):
            print("Built DEM")

        def buildOrthomosaic(self, **kwargs):
            print("Built orthomosaic")

        def exportPoints(self, path, **kwargs):
            print(f"Exported points to {path}")

        def exportRaster(self, path, **kwargs):
            print(f"Exported raster to {path}")

        def exportReport(self, path, **kwargs):
            print(f"Exported report to {path}")

        def detectMarkers(self):
            print("Detected markers")

        def optimizeCameras(self, **kwargs):
            print("Optimized cameras")

        def addScalebar(self, marker1, marker2):
            print("Added scalebar")
            return MockMetashape.Scalebar()

        def updateTransform(self):
            print("Updated transform")

        def dense_cloud(self):
            return MockMetashape.DenseCloud()

    class PointCloud:
        def __init__(self):
            self.points = [MockMetashape.Point() for _ in range(100)]

        class Filter:
            ReconstructionUncertainty = "ReconstructionUncertainty"
            ProjectionAccuracy = "ProjectionAccuracy"
            ReprojectionError = "ReprojectionError"

            def init(self, chunk, criterion):
                print(f"Initialized filter with criterion: {criterion}")

            def values(self):
                return [1] * 100

            def removePoints(self, threshold):
                print(f"Removed points with threshold: {threshold}")

    class DenseCloud:
        def setConfidenceFilter(self, min_val, max_val):
            print("Set confidence filter")

        def removePoints(self, points):
            print("Removed points")

        def resetFilters(self):
            print("Reset filters")

    class Transform:
        scale = 1.0

    class Point:
        valid = True

    class Scalebar:
        reference = type("Reference", (object,), {"distance": 0, "enabled": False})

Metashape = MockMetashape()  # Use mock Metashape for testing

@Gooey(
    dump_build_config=True,
    program_name="ReefMapper - SfM Batch Processing Toolbox",
    default_size=(900, 600),
    console=True,
    richtext_controls=True,
    shutdown_signal=signal.CTRL_C_EVENT,
    progress_regex=r"^progress: (?P<current>\d+)/(?P<total>\d+)$",
    progress_expr="current / total * 100",
    hide_progress_msg=True,
    timing_options={
        'show_time_remaining': True,
        'hide_time_remaining_on_complete': True,
    }
)
def main():
    desc = 'Batch process Structure from Motion (SfM) data using Agisoft Metashape.'
    parser = GooeyParser(description=desc)

    subs = parser.add_subparsers(help='commands', dest='command')

    # ------------------------------------------------------------------------------------------------------------------
    # Initialize
    # ------------------------------------------------------------------------------------------------------------------
    initialize_parser = subs.add_parser('Initialize')
    initialize_panel = initialize_parser.add_argument_group('Initialize Project', 
                                                            'Initialize the Metashape project', 
                                                            gooey_options={'show_border': True})

    initialize_panel.add_argument('root_path', type=str, widget='DirChooser', help='Root path for the project')
    initialize_panel.add_argument('folder_name', type=str, help='Folder name for the project')
    initialize_panel.add_argument('start_step', type=int, choices=range(1, 8), help='Starting step')
    
    # ------------------------------------------------------------------------------------------------------------------
    # Add and Align Photos
    # ------------------------------------------------------------------------------------------------------------------
    add_align_parser = subs.add_parser('Add and Align Photos')
    add_align_panel = add_align_parser.add_argument_group('Add and Align Photos', 
                                                          'Add photos and align them in Metashape', 
                                                          gooey_options={'show_border': True})

    add_align_panel.add_argument('root_path', type=str, widget='DirChooser', help='Root path for the project')
    add_align_panel.add_argument('folder_name', type=str, help='Folder name for the project')
    add_align_panel.add_argument('quality', type=float, help='Image quality threshold')
    
    # ------------------------------------------------------------------------------------------------------------------
    # Sparse Point Cloud Filtering
    # ------------------------------------------------------------------------------------------------------------------
    spc_filter_parser = subs.add_parser('Sparse Point Cloud Filtering')
    spc_filter_panel = spc_filter_parser.add_argument_group('Sparse Point Cloud Filtering', 
                                                            'Filter the sparse point cloud in Metashape', 
                                                            gooey_options={'show_border': True})

    spc_filter_panel.add_argument('chunk', type=str, help='Chunk object')
    
    # ------------------------------------------------------------------------------------------------------------------
    # Scaling
    # ------------------------------------------------------------------------------------------------------------------
    scaling_parser = subs.add_parser('Scaling')
    scaling_panel = scaling_parser.add_argument_group('Scaling', 
                                                      'Scale the project in Metashape', 
                                                      gooey_options={'show_border': True})

    scaling_panel.add_argument('chunk', type=str, help='Chunk object')
    scaling_panel.add_argument('valid_markers', type=str, help='Comma-separated list of valid markers')

    # ------------------------------------------------------------------------------------------------------------------
    # Error Reduction, Part 2
    # ------------------------------------------------------------------------------------------------------------------
    error_reduction_parser = subs.add_parser('Error Reduction')
    error_reduction_panel = error_reduction_parser.add_argument_group('Error Reduction', 
                                                                      'Reduce errors in Metashape', 
                                                                      gooey_options={'show_border': True})

    error_reduction_panel.add_argument('chunk', type=str, help='Chunk object')
    
    # ------------------------------------------------------------------------------------------------------------------
    # Build Dense Cloud
    # ------------------------------------------------------------------------------------------------------------------
    dense_cloud_parser = subs.add_parser('Build Dense Cloud')
    dense_cloud_panel = dense_cloud_parser.add_argument_group('Build Dense Cloud', 
                                                              'Build the dense cloud in Metashape', 
                                                              gooey_options={'show_border': True})

    dense_cloud_panel.add_argument('chunk', type=str, help='Chunk object')
    
    # ------------------------------------------------------------------------------------------------------------------
    # Build and Export DEM and Orthomosaic
    # ------------------------------------------------------------------------------------------------------------------
    dem_ortho_parser = subs.add_parser('Build and Export DEM and Orthomosaic')
    dem_ortho_panel = dem_ortho_parser.add_argument_group('Build and Export DEM and Orthomosaic', 
                                                          'Build and export DEM and orthomosaic in Metashape', 
                                                          gooey_options={'show_border': True})

    dem_ortho_panel.add_argument('chunk', type=str, help='Chunk object')
    dem_ortho_panel.add_argument('folder_name', type=str, help='Folder name for the project')
    dem_ortho_panel.add_argument('prod_path', type=str, widget='DirChooser', help='Production path for the project')
    dem_ortho_panel.add_argument('survey_year', type=str, help='Survey year for export')

    # ------------------------------------------------------------------------------------------------------------------
    # All Steps
    # ------------------------------------------------------------------------------------------------------------------
    all_steps_parser = subs.add_parser('All Steps')
    all_steps_panel = all_steps_parser.add_argument_group('All Steps', 
                                                          'Run all processing steps from start to finish', 
                                                          gooey_options={'show_border': True})

    all_steps_panel.add_argument('process_log', type=str, widget='FileChooser', help='Path to the process log CSV file')
    all_steps_panel.add_argument('batch_no', type=int, help='Batch number to process')
    all_steps_panel.add_argument('--quality_threshold', type=float, default=0.6, help='Image quality threshold')
    all_steps_panel.add_argument('--survey_year', type=str, default='2024', help='Survey year for export')
    all_steps_panel.add_argument('--valid_markers', type=str, help='Comma-separated list of valid markers', default='')

    # ------------------------------------------------------------------------------------------------------------------
    # Execute
    # ------------------------------------------------------------------------------------------------------------------
    args = parser.parse_args()

    if args.command == 'Initialize':
        doc, chunk, prod_path, opf = initialize_project(args.root_path, args.folder_name, args.start_step)
    
    elif args.command == 'Add and Align Photos':
        root_path = args.root_path
        folder_name = args.folder_name
        quality = args.quality
        chunk = Metashape.app.document.chunk  # Assuming chunk is the currently active chunk
        opf = open(os.path.join(root_path, 'Products', f'{folder_name}_readme.txt'), 'a')
        add_and_align_photos(root_path, folder_name, chunk, opf, quality)
    
    elif args.command == 'Sparse Point Cloud Filtering':
        chunk = Metashape.app.document.chunk  # Assuming chunk is the currently active chunk
        opf = open(os.path.join(chunk.path, 'Products', f'{chunk.label}_readme.txt'), 'a')
        sparse_point_cloud_filtering(chunk, opf)
    
    elif args.command == 'Scaling':
        chunk = Metashape.app.document.chunk  # Assuming chunk is the currently active chunk
        opf = open(os.path.join(chunk.path, 'Products', f'{chunk.label}_readme.txt'), 'a')
        valid_markers = args.valid_markers.split(',')
        scaling(chunk, opf, valid_markers)
    
    elif args.command == 'Error Reduction':
        chunk = Metashape.app.document.chunk  # Assuming chunk is the currently active chunk
        opf = open(os.path.join(chunk.path, 'Products', f'{chunk.label}_readme.txt'), 'a')
        error_reduction_part2(chunk, opf)
    
    elif args.command == 'Build Dense Cloud':
        chunk = Metashape.app.document.chunk  # Assuming chunk is the currently active chunk
        build_dense_cloud(chunk)
    
    elif args.command == 'Build and Export DEM and Orthomosaic':
        chunk = Metashape.app.document.chunk  # Assuming chunk is the currently active chunk
        build_and_export_dem_orthomosaic(chunk, args.folder_name, args.prod_path, args.survey_year)
    
    elif args.command == 'All Steps':
        valid_markers = [f'target {marker.strip()}' for marker in args.valid_markers.split(',')] if args.valid_markers else []
        process_log = args.process_log
        batch_no = args.batch_no
        quality_thres = args.quality_threshold
        survey_year = args.survey_year

        with open(process_log, newline='') as csvfile:
            spamreader = csv.reader(csvfile, delimiter=',')
            for col in spamreader:
                if col[2] == str(batch_no):
                    fpath = col[0]
                    if os.path.exists(fpath):
                        metashape_process(fpath, col[1], int(col[3]), int(col[4]), valid_markers, quality_thres, survey_year)
                    else:
                        print(f'{fpath} does not exist, move to next project!')

def metashape_process(fpath, folder_name, start_step, end_step, valid_markers, quality, survey_year):
    root_path = os.path.dirname(fpath)
    doc, chunk, prod_path, opf = initialize_project(root_path, folder_name, start_step)
    
    if start_step <= 2 <= end_step:
        add_and_align_photos(root_path, folder_name, chunk, opf, quality)
    
    if start_step <= 3 <= end_step:
        sparse_point_cloud_filtering(chunk, opf)
    
    if start_step <= 4 <= end_step:
        scaling(chunk, opf, valid_markers)
    
    if start_step <= 5 <= end_step:
        error_reduction_part2(chunk, opf)
    
    if start_step <= 6 <= end_step:
        build_dense_cloud(chunk)
    
    if start_step <= 7 <= end_step:
        build_and_export_dem_orthomosaic(chunk, folder_name, prod_path, survey_year)
    
    doc.save()
    opf.close()
    print('Finished processing:', folder_name)

def initialize_project(root_path, folder_name, start_step):
    underline = 50 * '-'
    doc = Metashape.app.document
    doc.clear()
    
    prod_path = os.path.join(root_path, 'Products')
    psxfile = os.path.join(prod_path, f'{folder_name}.psx')
    
    if not os.path.exists(prod_path):
        os.mkdir(prod_path)
    
    Metashape.app.settings.log_enable = True
    log_file = os.path.join(prod_path, f'{folder_name}_log.txt')
    Metashape.app.settings.log_path = log_file
    
    if start_step == 1:
        opf = open(os.path.join(prod_path, f'{folder_name}_readme.txt'), 'w')
        if os.path.exists(psxfile):
            os.remove(psxfile)
        if os.path.exists(log_file):
            os.remove(log_file)
        file_path = os.path.join(prod_path, f'{folder_name}.files')
        if os.path.exists(file_path):
            shutil.rmtree(file_path, ignore_errors=True)
        
        doc.save(psxfile)
        chunk = doc.addChunk()
    else:
        opf = open(os.path.join(prod_path, f'{folder_name}_readme.txt'), 'a')
        doc.open(psxfile)
        chunk = doc.chunk
    
    chunk.camera_location_accuracy = Metashape.Vector((0.1, 0.1, 0.15))
    opf.write(f'Readme file for {folder_name}\n{underline}\n')
    return doc, chunk, prod_path, opf

def get_photo_list(root_path, extra=''):
    photo_list = []
    for dirpath, _, filenames in os.walk(root_path):
        for filename in filenames:
            if filename.endswith((".jpg", ".jpeg", ".tif", ".tiff", ".png")):
                photo_list.append(os.path.join(dirpath, filename))
    return photo_list

def add_and_align_photos(root_path, folder_name, chunk, opf, quality):
    underline = 50 * '-'
    opf.write(f'\nAlign photos\n{underline}\n')
    
    photo_list = get_photo_list(os.path.join(root_path, ''))
    n = len(photo_list)
    chunk.addPhotos(photo_list)
    
    chunk.analyzePhotos(chunk.cameras)
    bad_quality = 0
    quality_log = [0] * n
    
    for qc, camera in enumerate(chunk.cameras):
        if float(camera.meta['Image/Quality']) < quality and quality_log[qc - 1] == 0:
            bad_quality += 1
            quality_log[qc] = 1
            camera.enabled = False
            opf.write(f'Image {camera.label} with quality {float(camera.meta["Image/Quality"]):.3f} disabled\n')
    
    n_enabled = n - bad_quality
    opf.write(f'{bad_quality} photo(s) ({round(bad_quality / n * 100, 1)}%) of quality below {quality}\n')
    
    chunk.matchPhotos(generic_preselection=True, reference_preselection=False, filter_mask=False, keypoint_limit=40000, tiepoint_limit=0)
    chunk.alignCameras(adaptive_fitting=False)
    
    counter = sum(1 for camera in chunk.cameras if camera.transform)
    opf.write(f'Enabled images: {n_enabled}\nAligned images: {counter}\n')
    
    if n_enabled - counter > 15:
        print('Try again or reevaluate the quality of corresponding photos')
        quit()
    
    shutil.copyfile(os.path.join(prod_path, f'{folder_name}.psx'), os.path.join(prod_path, f'{folder_name}_bkup.psx'))
    shutil.copytree(os.path.join(prod_path, f'{folder_name}.files'), os.path.join(prod_path, f'{folder_name}_bkup.files'))

def sparse_point_cloud_filtering(chunk, opf):
    underline = 50 * '-'
    opf.write(f'\nSparse point cloud filtering\n{underline}\n')
    keep_percent = 51
    points = chunk.point_cloud.points
    
    chunk.optimizeCameras(tiepoint_covariance=True)
    f = Metashape.PointCloud.Filter()
    f.init(chunk, criterion=Metashape.PointCloud.Filter.ReconstructionUncertainty)
    list_values = sorted([f.values[i] for i in range(len(f.values)) if points[i].valid])
    RecUncert = list_values[int(len(list_values) * keep_percent / 100)]
    RecUncert = max(RecUncert, 10)
    f.removePoints(RecUncert)
    opf.write(f'Reconstruction Uncertainty threshold to keep 51%: {RecUncert:.2f}\n')
    
    chunk.optimizeCameras(tiepoint_covariance=True)
    f.init(chunk, criterion=Metashape.PointCloud.Filter.ProjectionAccuracy)
    list_values = sorted([f.values[i] for i in range(len(f.values)) if points[i].valid])
    ProjAcc = list_values[int(len(list_values) * keep_percent / 100)]
    ProjAcc = max(ProjAcc, 2)
    f.removePoints(ProjAcc)
    opf.write(f'Projection Accuracy threshold to keep 51%: {ProjAcc:.2f}\n')
    
    chunk.optimizeCameras(tiepoint_covariance=True)

def scaling(chunk, opf, valid_markers):
    underline = 50 * '-'
    opf.write(f'\nScaling\n{underline}\n')
    chunk.detectMarkers()
    error_thresh = 0.4
    MarkersList = [marker for marker in chunk.markers if marker.label in valid_markers and marker.position]
    
    for marker in MarkersList:
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
                total = (total[0] + error ** 2, total[1] + 1)
                cam_err.append((camera.label, error))
            pix_error = math.sqrt(total[0] / total[1])
            if pix_error >= error_thresh:
                max_err = sorted(cam_err, key=lambda x: x[1], reverse=True)[0]
                for camera in marker.projections.keys():
                    if camera.label == max_err[0]:
                        marker.projections[camera] = None
                opf.write(f'Removed {max_err[0]} with pix error {max_err[1]:.4}\n')
        opf.write(f'Marker: {marker.label}, projections: {len(cam_err)}, total error {pix_error:.4f} pix\n{underline}\n\n')
    
    sb_dist = 0.25
    for i, marker in enumerate(MarkersList):
        if i % 2 > 0:
            sb = chunk.addScalebar(chunk.markers[i - 1], chunk.markers[i])
            sb.reference.distance = sb_dist
    chunk.updateTransform()
    
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
        opf.write(f'Scalebar: {scalebar.label}, distance: {dist_source:.2f}, estimated distance: {dist_estimated:.5f}, error: {dist_error:.6f}\n')
        if dist_error >= 0.002:
            opf.write('Scalebar error too high: lower error_thresh and run again starting at step 4!\n')
        scalebar.reference.enabled = False

def error_reduction_part2(chunk, opf):
    underline = 50 * '-'
    opf.write(f'\nError Reduction\n{underline}\n')
    keep_percent = 90
    points = chunk.point_cloud.points
    f = Metashape.PointCloud.Filter()
    f.init(chunk, criterion=Metashape.PointCloud.Filter.ReprojectionError)
    list_values = sorted([f.values[i] for i in range(len(f.values)) if points[i].valid])
    Rep_Err = list_values[int(len(list_values) * keep_percent / 100)]
    f.removePoints(Rep_Err)
    opf.write(f'Reprojection Error threshold to keep 90%: {Rep_Err:.2f}\n')
    chunk.optimizeCameras()

def build_dense_cloud(chunk):
    chunk.buildDepthMaps(downscale=4, filter_mode=Metashape.MildFiltering)
    chunk.buildDenseCloud(point_colors=True, point_confidence=True)
    chunk.dense_cloud.setConfidenceFilter(0, 1)
    chunk.dense_cloud.removePoints(list(range(128)))
    chunk.dense_cloud.resetFilters()

def build_and_export_dem_orthomosaic(chunk, folder_name, prod_path, survey_year):
    arc_path = os.path.join(prod_path, 'ARC')
    if not os.path.exists(arc_path):
        os.mkdir(arc_path)
    
    chunk.buildDem(source_data=Metashape.DenseCloudData)
    chunk.buildOrthomosaic(surface_data=Metashape.ElevationData, fill_holes=True, ghosting_filter=False, refine_seamlines=False, resolution=0.0005)
    
    out_projection = Metashape.OrthoProjection()
    out_projection.crs = Metashape.CoordinateSystem("EPSG::32604")
    
    chunk.exportRaster(os.path.join(arc_path, f'{survey_year}_{folder_name}_dem.tif'), projection=out_projection, resolution=0.001, save_world=True, source_data=Metashape.ElevationData)
    chunk.exportRaster(os.path.join(arc_path, f'{survey_year}_{folder_name}_dem_1cm.tif'), projection=out_projection, resolution=0.01, save_world=True, source_data=Metashape.ElevationData)
    
    compression = Metashape.ImageCompression()
    compression.tiff_compression = Metashape.ImageCompression.TiffCompressionLZW
    compression.jpeg_quality = 99
    chunk.exportRaster(os.path.join(arc_path, f'{survey_year}_{folder_name}_mos.tif'), projection=out_projection, resolution=0.0005, image_compression=compression, save_world=True, save_alpha=False, source_data=Metashape.OrthomosaicData)
    
    url = os.path.join(arc_path, f'{survey_year}_{folder_name}_rpt.html')
    chunk.exportReport(url, title=f'{survey_year}_{folder_name}', description="Processing Report")
    
    with open(url) as f:
        html = f.read()
    h = html2text.HTML2Text()
    t = h.handle(html)
    lines = t.splitlines()
    
    with open(os.path.join(arc_path, f'{url[:-5]}.csv'), 'w', newline='') as file:
        writer = csv.writer(file)
        for line in lines:
            words = line.strip().split('|')
            writer.writerow(words)
    
    cams = chunk.cameras
    proj_path = chunk.document.path
    proj_dir, proj_name = os.path.split(proj_path)
    proj_name = proj_name[:-4]
    outputs = {}
    cams_filename = os.path.join(proj_dir, f'{proj_name}.cams.xml')
    meta_filename = os.path.join(proj_dir, f'{proj_name}.meta.json')
    chunk.exportCameras(cams_filename)
    
    with open(meta_filename, 'w') as meta_file:
        for cam in cams:
            key = cam.key
            path = cam.photo.path
            center = list(cam.center) if cam.center else None
            agi_trans = cam.transform
            trans = [list(agi_trans.row(n)) for n in range(agi_trans.size[1])] if agi_trans else None
            outputs[key] = {'path': path, 'center': center, 'transform': trans}
        meta_file.write(json.dumps({'cameras': outputs}, indent=4))
    
    pt_file = os.path.join(prod_path, f'{folder_name}.ply')
    if not os.path.isfile(pt_file):
        chunk.exportPoints(pt_file, source_data=Metashape.DenseCloudData)

if __name__ == '__main__':
    main()
