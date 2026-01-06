
# step_02_add_align.py
import os
import shutil
from pathlib import Path
import Metashape  # Ensure Metashape Python API is available

# ---- Tunable parameters ----
VALID_EXTS = {".jpg", ".jpeg", ".tif", ".tiff", ".png"}
MIN_SIZE_BYTES = 1024          # Skip tiny files (helps avoid AppleDouble stubs)
QUALITY_THRESHOLD = 0.5        # Disable images with quality below this
KEYPOINT_LIMIT = 40000
TIEPOINT_LIMIT = 0
DOWNSCALE = 2
THRESH_ALIGN = 15              # Max number of enabled images allowed to fail alignment

def _get_photo_list(root_path: str):
    """
    Collect a list of image file paths from a folder.
    Skips hidden files and AppleDouble sidecars (names starting with '.' or '._').
    Restricts to VALID_EXTS and min file size.
    """
    photos = []
    root = Path(root_path)

    # If you need recursive, replace 'iterdir()' with 'rglob("*")'
    for p in root.iterdir():
        if not p.is_file():
            continue
        name = p.name
        if name.startswith(".") or name.startswith("._"):
            # Skip hidden and AppleDouble sidecar files
            continue
        if p.suffix.lower() in VALID_EXTS and p.stat().st_size >= MIN_SIZE_BYTES:
            photos.append(str(p))
    return photos


def run(job, doc, chunk, db, logger, opf, log_path, products_dir, psx_path):
    """
    Step 2: Add and align photos
    Compatible with your pipeline: logs progress, handles failures, and makes a backup.
    """
    # --- Step header ---
    logger.info("Step 2. Add and align photos")
    print("Step 2. Add and align photos")
    underline = "=" * 50
    opf.write(f"\nAlign photos\n{underline}\n")

    # --- Gather images ---
    # Prefer job.project_path if you’re staging images there; adjust as needed.
    image_root = getattr(job, "project_path", None) or getattr(job, "images_path", None) or os.path.dirname(psx_path)
    photo_list = _get_photo_list(image_root)
    n_total = len(photo_list)

    if n_total == 0:
        msg = f"No valid images found in {image_root}"
        logger.error(msg)
        opf.write(msg + "\n")
        raise RuntimeError(msg)

    # --- Add photos ---
    logger.info(f"Adding {n_total} photos from {image_root}")
    print(f"Adding {n_total} photos")
    opf.write(f"Adding {n_total} photos\n")
    chunk.addPhotos(photo_list)

    # --- Estimate image quality & disable poor ones ---
    # This writes a quality metric to camera.meta['Image/Quality'] for each camera
    chunk.analyzeImages(chunk.cameras)

    bad_quality = 0
    # Logs: 0 = enabled (good), 1 = disabled (poor quality)
    quality_log = []

    for camera in chunk.cameras:
        # Some cameras may not have quality for various reasons; guard it

        meta = camera.meta
        q_val = None
        if meta is not None and "Image/Quality" in meta:
            try:
                q_val = float(meta["Image/Quality"])
            except Exception:
                q_val = None

        # Disable if below threshold and quality metric is available
        if q_val is not None and q_val < QUALITY_THRESHOLD:
            camera.enabled = False
            bad_quality += 1
            quality_log.append(1)
            logger.debug(f"Disabled {camera.label} (quality {q_val:.3f})")
            opf.write(f"Image {camera.label} with quality {q_val:.3f} disabled\n")
        else:
            quality_log.append(0)

    n_enabled = len(chunk.cameras) - bad_quality
    pct_disabled = (bad_quality / len(chunk.cameras) * 100.0) if chunk.cameras else 0.0

    print(f"{bad_quality} photo(s) ({pct_disabled:.1f}%) below quality {QUALITY_THRESHOLD}")
    opf.write(f"{bad_quality} photo(s) ({pct_disabled:.1f}%) of quality below {QUALITY_THRESHOLD}\n")
    logger.info(f"Enabled images after quality filter: {n_enabled}")
    # Optional: write the quality log vector for diagnostics
    opf.write(f"Quality log (0=good,1=poor): {quality_log}\n")

    # --- Match photos ---
    # Keep generic preselection on; reference preselection off unless you have EXIF/GPS you trust.
    logger.info("Matching photos...")
    chunk.matchPhotos(
        generic_preselection=True,
        reference_preselection=False,
        filter_mask=False,
        keypoint_limit=KEYPOINT_LIMIT,
        tiepoint_limit=TIEPOINT_LIMIT,
        downscale=DOWNSCALE,
    )

    # --- Align cameras ---
    logger.info("Aligning cameras...")
    # Disable adaptive fitting for deterministic distortion behavior
    chunk.alignCameras(adaptive_fitting=False)

    # --- Count aligned cameras ---
    aligned_count = 0
    for camera in chunk.cameras:
        # Count only enabled cameras that have a transform (aligned)
        if camera.enabled and camera.transform:
            aligned_count += 1

    print(f"Enabled images: {n_enabled}")
    print(f"Aligned images: {aligned_count}")
    opf.write(f"Enabled images: {n_enabled}\nAligned images: {aligned_count}\n")

    not_aligned = n_enabled - aligned_count
    if not_aligned == 0:
        print("All enabled images aligned!")
        logger.info("All enabled images aligned.")
    elif 0 < not_aligned <= THRESH_ALIGN:
        msg = f"{not_aligned} enabled images did not align."
        print(msg)
        logger.warning(msg)
        opf.write(msg + "\n")
    else:
        # Too many failures—signal to pipeline controller
        msg = (
            f"{not_aligned} enabled images failed to align (threshold {THRESH_ALIGN}). "
            f"Consider re-running with lower accuracy/downscale, revisiting quality, or "
            f"checking overlap."
        )
        print(msg)
        logger.error(msg)
        opf.write(msg + "\n")
        # Persist state before failing
        try:
            doc.save(psx_path)
        except Exception as e:
            logger.warning(f"Save before abort failed: {e}")
        raise RuntimeError(msg)

    # --- Save & make a backup of the project ---
    try:
        doc.save(psx_path)
    except Exception as e:
        logger.warning(f"Primary save failed: {e}")

    try:
        psx = Path(psx_path)
        prod_path = Path(products_dir)
        folder_name = psx.stem

        backup_psx = prod_path / f"{folder_name}_bkup.psx"
        src_files_dir = psx.with_suffix(".files")
        dst_files_dir = prod_path / f"{folder_name}_bkup.files"

        shutil.copyfile(psx, backup_psx)
        if dst_files_dir.exists():
            shutil.rmtree(dst_files_dir, ignore_errors=True)
        shutil.copytree(src_files_dir, dst_files_dir)

        logger.info(f"Backup written to {backup_psx} and {dst_files_dir}")
        opf.write(f"Backup written to {backup_psx} and {dst_files_dir}\n")
    except Exception as e:
        # Backup failure should not crash the pipeline; log and continue
        logger.warning(f"Backup failed: {e}")
        opf.write(f"Backup failed: {e}\n")

    print("Step 2. Add and align photos finished")
