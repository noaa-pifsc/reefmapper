import Metashape
import logging
from SfMBatchProcess_20251002 import orient_model   # wherever your function lives

logger = logging.getLogger("reefmapper")
logger.setLevel(logging.DEBUG)

doc = Metashape.Document()
doc.open("N:/StRS_Sites/2022/RA2201_MARI/ASC/ASC-637/Products_automation/ASC-637.psx")
chunk = doc.chunk

orient_model(chunk, logger)
doc.save()