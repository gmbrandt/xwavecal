"""
Module for basic reduction procedures such as overscan subtraction, overscan trimming and gain normalization.

"""

import numpy as np

from echelle.stages import Stage
from echelle.utils.fits_utils import parse_region_keyword


class OverscanSubtractor(Stage):
    def __init__(self, runtime_context):
        super(OverscanSubtractor, self).__init__(runtime_context)

    def do_stage(self, image):
        image.data = np.array(image.data.copy(order='C'), dtype=float)
        data_section = parse_region_keyword(image.get_header_val('data_section'))
        overscan_section = parse_region_keyword(image.get_header_val('overscan_section'))
        image.data[data_section] -= np.median(image.data[overscan_section])
        return image


class GainNormalizer(Stage):
    def __init__(self, runtime_context):
        super(GainNormalizer, self).__init__(runtime_context)

    def do_stage(self, image):
        image.data *= image.get_header_val('gain')
        # set gain in image to 1.
        return image


class Trimmer(Stage):
    def __init__(self, runtime_context):
        super(Trimmer, self).__init__(runtime_context)

    def do_stage(self, image):
        data_section = parse_region_keyword(image.get_header_val('data_section'))
        image.data = image.data[data_section]
        return image
