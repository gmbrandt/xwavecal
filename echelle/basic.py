"""
Module for basic reduction procedures such as overscan subtraction,
overscan trimming and gain normalization.

Author:
        G. Mirek Brandt

"""

import numpy as np

from echelle.stages import Stage
from echelle.utils.fits_utils import parse_region_keyword
import sep
import logging as logger


class OverscanSubtractor(Stage):
    def __init__(self, runtime_context):
        super(OverscanSubtractor, self).__init__(runtime_context)

    def do_stage(self, image):
        logger.info('Subtracting by median overscan')
        image.data = np.ascontiguousarray(image.data.astype(float))
        data_section = parse_region_keyword(image.get_header_val('data_section'))
        overscan_section = parse_region_keyword(image.get_header_val('overscan_section'))
        image.data[data_section] -= np.median(image.data[overscan_section])
        return image


class GainNormalizer(Stage):
    def __init__(self, runtime_context):
        super(GainNormalizer, self).__init__(runtime_context)

    def do_stage(self, image):
        logger.info('Multiplying by gain')
        image.data *= image.get_header_val('gain')
        image.set_header_val('gain', 1.0)
        return image


class Trimmer(Stage):
    def __init__(self, runtime_context):
        super(Trimmer, self).__init__(runtime_context)

    def do_stage(self, image):
        logger.info('Trimming image to {0}'.format(image.get_header_val('data_section')))
        data_section = parse_region_keyword(image.get_header_val('data_section'))
        image.data = image.data[data_section]
        image.data = np.ascontiguousarray(image.data)
        return image


class BackgroundSubtract(Stage):
    def __init__(self, runtime_context):
        super(BackgroundSubtract, self).__init__(runtime_context)

    def do_stage(self, image):
        logger.info('Background subtracting the 2d frame')
        image.data = image.data - sep.Background(image.data).back()
        return image
