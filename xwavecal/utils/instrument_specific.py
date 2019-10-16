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


def parse_harps_region_keyword(key):
    # WARN this function assumes all HARPS frames have overscans in x only (and zero overscan in y).
    """
    :param key: tuple of the form (rx, ry, osx, osy, psx, psy, Nx, Ny)
    where each entry is an integer.
    osx and osy are the number of overscan pixels in x and y (not the extent, the number of pixels)
    osx and osy are the number of prescan pixels in x and y
    Nx, Ny are the full size of the image (i.e. NAXIS1 and NAXIS2)
    rx and ry are the number of region pixels in x and y. E.g. for the overscan region, rx, ry will equal osx, osy.
    :return: tuple
             (slice(y1,y2,None), slice(x1,x2,None))

             the rows (y1:y2) and columns (x1:x2) corresponding to the input region Nx, Ny.
             Note x2!= Nx because of how HARPS frames store the overscan information.
    """
    nx, ny, osx, osy, psx, psy, Nx, Ny = key
    if np.isclose(nx, Nx, atol=300):
        # we are processing the data region
        return slice(psy, psy + ny), slice(psx, psx + nx)
    else:
        # we are processing the overscan region (because we do not subtract by the prescan anywhere in this pipeline.
        return slice(0, Ny), slice(Nx - osx, Nx)
