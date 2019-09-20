import os
import tempfile
import shutil
import numpy as np

import logging as logger


class Translator(dict):
    def __init__(self, header_keys, type_keys):
        self.header_keys = header_keys
        self.type_keys = type_keys

    def __getitem__(self, key):
        return self.header_keys.get(key, key)


def writeto(hdu_list, filepath, fpack=False, overwrite=True, output_verify='fix+warn', quant=64):
    if not fpack:
        hdu_list.writeto(filepath.split('.fz')[0],
                         overwrite=overwrite, output_verify=output_verify)
    else:
        if not filepath.endswith('.fz'):
            filepath += '.fz'
        logger.info('Writing file to {filepath}'.format(filepath=filepath))
        base_filename = os.path.basename(filepath).split('.fz')[0]
        with tempfile.TemporaryDirectory() as temp_directory:
            hdu_list.writeto(os.path.join(temp_directory, base_filename),
                             overwrite=overwrite, output_verify=output_verify)
            if os.path.exists(filepath):
                os.remove(filepath)
            command = 'fpack -q {quantization} {temp_directory}/{basename}'
            os.system(command.format(quantization=int(quant), temp_directory=temp_directory, basename=base_filename))
            base_filename += '.fz'
            shutil.move(os.path.join(temp_directory, base_filename), filepath)


def parse_region_keyword(key, index_from_one=True):
    """
    :param key: string of the form '[x:y],[x:y]'
    :return: tuple
             (slice(x,y,None), slice(x,y,None))
             or if index_from_one:
             (slice(x-1,y,None), slice(x-1,y,None))
    """
    boundaries = [np.array(i.split(':')).astype(int) for i in key.replace('[', '').replace(']', '').split(',')]
    if index_from_one:
        boundaries[0][0], boundaries[1][0] = boundaries[0][0] - 1, boundaries[1][0] - 1
    return tuple(slice(*boundary) for boundary in boundaries)
