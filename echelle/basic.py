"""
Module for basic reduction steps like overscan subtraction,
overscan trimming and gain normalization.

Author:
        G. Mirek Brandt

"""

import numpy as np
import sep

from echelle.stages import Stage
from echelle.utils.runtime_utils import import_obj
from echelle.utils.basic_utils import median_subtract_channels_y

import logging
logger = logging.getLogger(__name__)


class OverscanSubtractor(Stage):
    def __init__(self, runtime_context):
        super(OverscanSubtractor, self).__init__(runtime_context)

    def do_stage(self, image):
        parse_region_keyword = import_obj(self.runtime_context.parse_region_keyword)

        logger.info('Subtracting by median overscan')
        image.data = np.ascontiguousarray(image.data.astype(float))
        data_section = parse_region_keyword(image.get_header_val('data_section'))
        overscan_section = parse_region_keyword(image.get_header_val('overscan_section'))

        logger.info('data section (Y, X) is {0} and overscan section (Y, X) is {1}'.format(data_section, overscan_section))
        image.data[data_section] -= np.median(image.data[overscan_section])
        return image


class GainNormalizer(Stage):
    def __init__(self, runtime_context):
        super(GainNormalizer, self).__init__(runtime_context)

    def do_stage(self, image):
        logger.info('Multiplying by gain')
        image.data = image.data.astype(float) * image.get_header_val('gain')
        image.set_header_val('gain', 1.0)
        return image


class Trimmer(Stage):
    def __init__(self, runtime_context):
        super(Trimmer, self).__init__(runtime_context)

    def do_stage(self, image):
        parse_region_keyword = import_obj(self.runtime_context.parse_region_keyword)

        data_section = parse_region_keyword(image.get_header_val('data_section'))
        logger.info('Trimming image to (Y, X) {0}'.format(data_section))
        image.data = image.data[data_section]
        image.data = np.ascontiguousarray(image.data)
        return image


class MedianSubtractReadoutsAlongY(Stage):
    def __init__(self, runtime_context):
        super(MedianSubtractReadoutsAlongY, self).__init__(runtime_context)

    def do_stage(self, image):
        logger.info('Subtracting the median from each readout channel')
        num = image.get_header_val('num_rd_channels')
        image.data = median_subtract_channels_y(image.data, num)
        return image


class BackgroundSubtract(Stage):
    def __init__(self, runtime_context):
        super(BackgroundSubtract, self).__init__(runtime_context)

    def do_stage(self, image):
        logger.info('Background subtracting the 2d frame')
        image.data = image.data.copy(order='C')
        background = sep.Background(image.data).back()
        image.data = image.data - background
        image.ivar = None if image.ivar is None else (image.ivar ** (-1) + np.abs(background)) ** (-1)
        del background
        return image


class BackgroundSubtractSpectrum(Stage):
    def __init__(self, runtime_context=None):
        super(BackgroundSubtractSpectrum, self).__init__(runtime_context=runtime_context)

    def do_stage(self, image):
        logger.info('Background subtracting the extracted 1d spectra')
        for key in [self.runtime_context.box_spectrum_name, self.runtime_context.blaze_corrected_spectrum_name]:
            if image.data_tables.get(key) is not None:
                spectrum = image.data_tables[key]
                if len(spectrum) > 0:
                    background = sep.Background(spectrum['flux'].data).back()
                    spectrum['flux'] -= background
                    spectrum['stderr'] = np.sqrt(spectrum['stderr']**2 + np.abs(background))
                    image.data_tables[key] = spectrum
        return image
