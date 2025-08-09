import os
import Metashape
import shutil

def run(job, doc, chunk, db, logger, opf, root_path=None, folder_name=None):
    """
    Step 2: Add and align photos in Metashape.
    """
    try:
        underline = 50 * '-'
        db.log_step(job.id, 2, "running", "Adding and aligning photos")
        print("Step 2. Add and align photos")
        logger.info("Step 2. Add and align photos")
        opf.write(f'\nAlign photos\n{underline}\n')
        db.set_start_step(job.id, 2)


        # Helper function to recursively get all photo paths from a folder
        def getPhotoList(path, photoList):
            for root, _, files in os.walk(path):
                for f in files:
                    if f.lower().endswith(('.jpg', '.jpeg', '.tif', '.tiff', '.png')):
                        photoList.append(os.path.join(root, f))

        # Get photo list from root_path
        photoList = []
        getPhotoList(os.path.join(root_path), photoList)
        n = len(photoList)

        if n == 0:
            raise Exception(f'No photos found in {root_path}')

        db.log_step(job.id, 2, "info", f"Adding {n} photos")
        print(f'Adding {n} photos')
        logger.info(f'Adding {n} photos')
        opf.write(f'Adding {n} photos\n')

        # Add photos to chunk
        chunk.addPhotos(photoList)
        doc.save(os.path.join(root_path, 'Products_automation', folder_name + '.psx'))

        # Analyze image quality
        db.log_step(job.id, 2, "info", "Analyzing image quality")
        chunk.analyzeImages(chunk.cameras)

        qualities = []
        bad_quality = 0
        quality_log = [0] * n  # 0 = good, 1 = poor

        for qc, camera in enumerate(chunk.cameras):
            img_quality = float(camera.meta.get('Image/Quality', 1.0))  # default 1.0 if missing
            qualities.append(img_quality)
            quality = 0.5

            if img_quality < quality and (qc == 0 or quality_log[qc - 1] == 0):
                bad_quality += 1
                quality_log[qc] = 1
                camera.enabled = False
                print(qc, camera.label, img_quality)
                logger.info(f'Image {camera.label} with quality {img_quality:.3f} disabled')
                opf.write(f'Image {camera.label} with quality {img_quality:.3f} disabled\n')

        if qualities:
            avg_quality = sum(qualities) / len(qualities)
            db.update_quality(job.id, "quality", avg_quality)
            logger.info(f"Average image quality: {avg_quality:.3f}")
            opf.write(f"Average image quality: {avg_quality:.3f}\n")
        else:
            logger.warning("No image qualities found to compute average.")

        n_enabled = n - bad_quality
        msg = f'{bad_quality} photos ({round(bad_quality / n * 100, 1)}%) below quality threshold {quality}'
        db.log_step(job.id, 2, "info", msg)
        print(msg)
        logger.info(msg)
        opf.write(msg + '\n')

        # Match and align photos
        db.log_step(job.id, 2, "info", "Matching and aligning photos")
        chunk.matchPhotos(generic_preselection=True, reference_preselection=False,
                          filter_mask=False, keypoint_limit=40000, tiepoint_limit=0)
        chunk.alignCameras(adaptive_fitting=False)

        # Verify alignment
        thresh_align = 15
        counter = sum(1 for camera in chunk.cameras if camera.transform)
        msg = f'Enabled: {n_enabled}, Aligned: {counter}'
        db.log_step(job.id, 2, "info", msg)
        print(msg)
        logger.info(msg)
        opf.write(msg + '\n')

        if n_enabled - counter == 0:
            db.log_step(job.id, 2, "info", "All enabled images aligned successfully")
        elif 0 < n_enabled - counter <= thresh_align:
            db.log_step(job.id, 2, "warning", f"{n_enabled - counter} images failed to align")
        else:
            raise Exception("Too many images failed to align")

        # Backup project file and associated .files folder
        prod_path = os.path.join(root_path, 'Products_automation')
        psxfile = os.path.join(prod_path, folder_name + '.psx')

        doc.save()
        shutil.copyfile(psxfile, os.path.join(prod_path, folder_name + '_bkup.psx'))
        bkup_files_path = os.path.join(prod_path, folder_name + '_bkup.files')
        orig_files_path = os.path.join(prod_path, folder_name + '.files')

        if os.path.exists(bkup_files_path):
            shutil.rmtree(bkup_files_path, ignore_errors=True)
        shutil.copytree(orig_files_path, bkup_files_path)

        db.log_step(job.id, 2, "completed", "Photo alignment complete")
        print("Step 2. Add and align photos finished")
        logger.info("Step 2. Add and align photos finished")
        db.set_end_step(job.id, 2)

    except Exception as e:
        msg = f"Step 2 failed: {str(e)}"
        logger.error(msg)
        db.log_step(job.id, 2, "error", msg)
        opf.write(msg + '\n')
        raise
