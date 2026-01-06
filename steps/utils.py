import os

PHOTO_EXTENSIONS = {'.jpg', '.jpeg', '.JPG', '.JPEG'}


import os

def get_photo_list(root_path):
    """Return list of photo paths from root directory, skipping hidden and AppleDouble files."""
    photo_list = []
    for root, _, files in os.walk(root_path):
        for photo in files:
            # Skip hidden files and AppleDouble sidecars
            if photo.startswith('.') or photo.startswith('._'):
                continue

            ext = os.path.splitext(photo)[1].lower()
            if ext in {'.jpg', '.jpeg'}:
                full_path = os.path.join(root, photo)
                if os.path.isfile(full_path):
                    photo_list.append(full_path)

def find_marker(label, chunk):
    """Find marker by label in chunk."""
    for marker in chunk.markers:
        if label == marker.label:
            return marker
    return None

def parse_marker_pairs(job):
    """Parse marker pairs from job and return valid marker labels."""
    valid_markers = []
    for pair in [job.marker_pair1, job.marker_pair2, job.marker_pair3, job.marker_pair4]:
        if pair and pair.lower() != 'na':
            try:
                t1, t2 = pair.split(',')
                valid_markers.append('target ' + t1.strip())
                valid_markers.append('target ' + t2.strip())
            except ValueError:
                continue
    return valid_markers
