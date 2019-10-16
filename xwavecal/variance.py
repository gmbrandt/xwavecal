"""
Module for calculating the inverse variances of the image.
"""

from xwavecal.stages import Stage
import numpy as np


class CalcInverseVariance(Stage):
    def __init__(self, runtime_context=None):
        super(CalcInverseVariance, self).__init__(runtime_context=runtime_context)

    def do_stage(self, image):
        var = estimate_var(image.data, image.get_header_val('read_noise'))
        image.ivar = var ** (-1)
        del var
        return image


def estimate_var(data, rdnoise):
    var = np.copy(data).astype(float)
    var[var < rdnoise ** 2] = 0
    return var + rdnoise ** 2
