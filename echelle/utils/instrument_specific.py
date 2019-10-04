"""
Module for any functions which are instrument specific and specified in the config files.
E.g. parse_region_keyword will parse the overscan and data region values in the .fits files
and return a tuple of splice objects specifying either the overscan or data region respectively.
"""

import numpy as np


def parse_nres_region_keyword(key, index_from_one=True):
    """
    :param key: string of the form '[y1:y2],[x1:x2]'
    :return: tuple
             (slice(y1,y2,None), slice(x1,x2,None))
             or if index_from_one:
             (slice(y1-1,y2,None), slice(x1-1,x2,None))
    """
    boundaries = [np.array(i.split(':')).astype(int) for i in key.replace('[', '').replace(']', '').split(',')]
    if index_from_one:
        boundaries[0][0], boundaries[1][0] = boundaries[0][0] - 1, boundaries[1][0] - 1
    return tuple(slice(*boundary) for boundary in boundaries)