"""
blaze.py: Module for making blaze calibration files for xwavecal spectrographs.

Authors
    G. Mirek Brandt (gmbrandt@ucsb.edu)
"""
import numpy as np
import logging
logger = logging.getLogger(__name__)
from copy import deepcopy

from xwavecal.stages import Stage, ApplyCalibration
from xwavecal.images import Image
from xwavecal.utils.blaze_utils import normalize_orders


class ApplyBlaze(ApplyCalibration):
    """
    Divides input images by the per-order normalized blaze spectrum.
    """
    def __init__(self, runtime_context):
        super(ApplyBlaze, self).__init__(runtime_context)

    @property
    def calibration_type(self):
        return 'BLAZE'

    def apply_master_calibration(self, image, master_calibration_path):
        blaze = Image.load(master_calibration_path, extension_name=self.runtime_context.blaze_name,
                           translator=image.translator)
        if not np.allclose(image.data.shape, blaze.data.shape):
            logger.error('Shape of blaze data and image data do not agree. Aborting blaze correction. Wavelength'
                         'solution may suffer.')  # pragma: no cover
        else:
            logger.info('Dividing by blaze from file {0}'.format(master_calibration_path))
            image.ivar = None if image.ivar is None else image.ivar * np.power(blaze.data, 2)
            # TODO full error propagation. The above does not hold for low (< 5) signal-to-noise pixels.
            image.data = np.divide(image.data, blaze.data, out=image.data.astype(float), where=~np.isclose(blaze.data, 0))
            image.set_header_val('IDBLAZE', (master_calibration_path, 'ID of blaze file used'))
        return image


class BlazeMaker(Stage):
    def __init__(self, runtime_context):
        super(BlazeMaker, self).__init__(runtime_context)

    @property
    def calibration_type(self):
        return 'BLAZE'

    def do_stage(self, image):
        logger.info('Making blaze file.')
        blaze = Image(data=deepcopy(image.data), translator=image.translator,
                      trace=image.trace, header=deepcopy(image.header), data_name=self.runtime_context.blaze_name)
        blaze.data = blaze.data / normalize_orders(blaze.data, blaze.trace, half_window=self.runtime_context.max_extraction_half_window)
        if image.ivar is not None:
            blaze.data[image.data * np.sqrt(image.ivar) < self.runtime_context.min_blaze_sn] = 0
        blaze.set_header_val('type', self.calibration_type.lower())
        return image, blaze
