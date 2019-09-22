"""
blaze.py: Driver script for making blaze calibration files for echelle spectrographs.

Authors
    G. Mirek Brandt (gmbrandt@ucsb.edu)
"""
import numpy as np
import sep
import logging as logger
from copy import deepcopy

from echelle.stages import Stage, ApplyCalibration
from echelle.images import Image
from echelle.utils.blaze_utils import normalize_orders


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
        logger.info('Dividing by blaze.')
        blaze = Image.load(master_calibration_path, extension_name=self.runtime_context.blaze_name)
        if not np.allclose(image.data.shape, blaze.data.shape):
            logger.error('Shape of blaze data and image data do not agree. Aborting blaze correction. Wavelength'
                         'solution may suffer.')
        else:
            image.data = image.data / blaze.data
            image.ivar = None if image.ivar is None else image.ivar * np.power(blaze.data, 2)
            import pdb
            pdb.set_trace()
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
        blaze.data = normalize_orders(blaze.data, blaze.trace, minval=300 * blaze.get_header_val('read_noise'),
                                      half_window=self.runtime_context.max_extraction_half_window)
        blaze.set_header_val('type', self.calibration_type.lower())
        return image, blaze


class BackgroundSubtractSpectrum(Stage):
    def __init__(self, runtime_context=None):
        super(BackgroundSubtractSpectrum, self).__init__(runtime_context=runtime_context)

    def do_stage(self, image):
        #TODO refactor storing information about which spectra exist.
        logger.info('Background subtracting the extracted 1d spectra')
        for key in [self.runtime_context.box_spectrum_name, self.runtime_context.blaze_corrected_spectrum_name]:
            if image.data_tables.get(key) is not None:
                spectrum = image.data_tables[key]
                spectrum['flux'] -= sep.Background(spectrum['flux'].data).back()
                image.data_tables[key] = spectrum
        return image
