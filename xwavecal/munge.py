"""
Module for rotating and flipping frames

This is to make frames such that the mean wavelengths of the diffraction orders increase
from left to right (e.g bluest horizontal coordinate is x=0 and reddest is x=4096)
and increase from top to bottom (e.g. bluest vertical coordinate is y=4096 and reddest is y=0).
"""

import numpy as np

from xwavecal.stages import Stage


class Rot90(Stage):
    """
    rotate image.data counter-clockwise by 90 degrees.
    """
    def __init__(self, runtime_context):
        super(Rot90, self).__init__(runtime_context)

    def do_stage(self, image):
        image.data = np.rot90(image.data)
        image.ivar = None if image.ivar is None else np.rot90(image.ivar)
        return image


class FlipHoriz(Stage):
    def __init__(self, runtime_context):
        super(FlipHoriz, self).__init__(runtime_context)

    def do_stage(self, image):
        image.data = np.flip(image.data, axis=1)
        image.ivar = None if image.ivar is None else np.flip(image.ivar, axis=1)
        return image


class FlipVert(Stage):
    def __init__(self, runtime_context):
        super(FlipVert, self).__init__(runtime_context)

    def do_stage(self, image):
        image.data = np.flip(image.data, axis=0)
        image.ivar = None if image.ivar is None else np.flip(image.ivar, axis=0)
        return image
