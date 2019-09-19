"""
blaze.py: Driver script for making blaze calibration files for echelle spectrographs.

Authors
    G. Mirek Brandt (gmbrandt@ucsb.edu)
"""
import numpy as np
import os
import copy

from echelle.stages import Stage, ApplyCalibration
from echelle.images import DataProduct
from echelle.fibers import lit_wavecal_fibers, lit_fibers

import echelle.settings as nres_settings
import sep

import logging as logger


# TODO divide the standard errors in the frame by the calibration as well, so that reduced_chi_sq in overlap
#  is accurate when S/N is low.
class ApplyBlaze(ApplyCalibration):
    #TODO do blaze correction with the 2d frame properly normalized, dont do it after extraction.
    """
    Loads the blaze spectrum from file and divides any arc fibers by the blaze.
    """
    def __init__(self, runtime_context):
        super(ApplyBlaze, self).__init__(runtime_context)

    @property
    def calibration_type(self):
        return 'BLAZE'

    def apply_master_calibration(self, image, master_calibration_path):
        logger.info('Image fiber information',
                    extra={'Lit fibers': str(lit_fibers(image)),
                                'Wavecal fibers': str(lit_wavecal_fibers(image))})
        blaze = DataProduct.load(master_calibration_path, extension_name=nres_settings.BLAZE_TABLE_NAME)
        master_filename = os.path.basename(master_calibration_path)
        image.set_header_val('L1IDBLAZ', (master_filename, 'ID of blaze file'))
        logger.info('Loaded blaze file',   extra={'L1IDBLAZ': image.get_header_val('L1IDBLAZ')})
        spectrum = copy.deepcopy(image.data_tables[nres_settings.BOX_SPECTRUM_NAME])
        blaze_spectrum = blaze.data
        if len(spectrum['id']) != len(blaze_spectrum['id']) or (not np.allclose(spectrum['id'].data,
                                                                                blaze_spectrum['id'].data)):
            logger.error('Trace IDs from blaze spectrum and spectrum do not agree. Aborting '
                         'blaze correction.', )
        else:
            # TODO only divide the arc lamp fibers. Right now this is only be valid for DOUBLE exposures.
            logger.info('Dividing by Blaze', )
            spectrum['flux'] = spectrum['flux'].data / blaze_spectrum['flux'].data
            spectrum.name = nres_settings.BLAZE_CORRECTED_BOX_SPECTRUM_NAME
            image.data_tables[nres_settings.BLAZE_CORRECTED_BOX_SPECTRUM_NAME] = spectrum
        return image

    def do_stage(self, image):
        master_calibration_path = self.get_calibration_filename(image)
        if master_calibration_path is None:
            self.on_missing_master_calibration(image)
            return image
        return self.apply_master_calibration(image, master_calibration_path)


class BackgroundSubtractSpectrum(Stage):
    def __init__(self, runtime_context=None):
        super(BackgroundSubtractSpectrum, self).__init__(runtime_context=runtime_context)

    def do_stage(self, image):
        if len(image.data_tables[nres_settings.BOX_SPECTRUM_NAME]['flux']) > 0:
            spectrum = image.data_tables[nres_settings.BOX_SPECTRUM_NAME]
            spectrum['flux'] -= sep.Background(spectrum['flux'].data).back()
        return image
