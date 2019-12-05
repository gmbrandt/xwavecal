import argparse
import glob
import os
import numpy as np
import importlib
from ast import literal_eval

from xwavecal.utils.fits_utils import Translator


def parse_args(args=None):
    parser = argparse.ArgumentParser(description='Reduce an xwavecal spectrograph frame.')
    parser.add_argument("--output-dir", required=True,
                        help="Directory within which to save the processed data files.")
    parser.add_argument("--input-dir", required=False, default=None,
                        help="Directory which contains the raw data files.")
    parser.add_argument('--data-paths', nargs='+', required=False, default=None,
                        help="path(s) to data, usage: '--data_paths path/to/first.fits path/to/second.fits'")
    parser.add_argument("--fpack", required=False, action='store_true',
                        help="fpack output files with the default quantization.")
    parser.add_argument("--config-file", required=True,
                        help="Path to the instrument specific configuration file.")
    parser.add_argument("--frame-type", required=False, default='any',
                        help="Frame type to either fit traces to or wavelength calibrate."
                             "Make sure frame type settings are appropriately set in the config file."
                             "lampflat files are used for tracing, wavecals are wavelength calibration"
                             "frames such as ThAr exposures. Must agree with the frame names in [stages],"
                             "e.g. lampflat, wavecal etc. Ignore to reduce all valid files",
                        type=str.lower)
    args = parser.parse_args(args)
    if args.data_paths is None and args.input_dir is None:
        raise ValueError('both input_dir and data_paths are None. Must specify raw data or a directory of raw data to process.')
    if not os.path.exists(args.config_file):
        raise FileNotFoundError('{0} not found.'.format(args.config_file))
    return args


def get_data_paths(dir_path, files_contain=None):
    all_files = glob.glob(os.path.join(dir_path, '*'))
    return [file for file in all_files if all([item in file for item in files_contain])]


def order_data(data_paths, data_class, primary_ext, header_keys, type_keys):
    translator = Translator(header_keys, type_keys)
    is_not_lampflat = lambda path: 0 if data_class.load(path, primary_ext, translator).get_header_val('type') == 'lampflat' else 1
    data_paths = list(data_paths)
    if len(data_paths) > 0:
        data_paths.sort(key=is_not_lampflat)
    return data_paths


def select_data_of_type(data_paths, data_class, primary_ext, header_keys, type_keys, frame_type='any'):
    is_type = lambda x: x == frame_type
    if frame_type == 'any':
        is_type = lambda x: type(x) is str
    translator = Translator(header_keys, type_keys)
    correct = lambda path: 1 if is_type(data_class.load(path, primary_ext, translator).get_header_val('type')) else 0
    return np.array(data_paths)[np.where([correct(path) for path in data_paths])]


def import_obj(full_class_string):
    """
    dynamically import a class or function from a string
    """

    class_data = full_class_string.split(".")
    module_path = ".".join(class_data[:-1])
    class_str = class_data[-1]

    module = importlib.import_module(module_path)
    return getattr(module, class_str)


def safe_eval(item):
    """
    :param item: str, int or dict
    :return: Any strings have erroneous leading " or ' removed.
    """
    out = literal_eval(item)
    if isinstance(out, str):
        return out.replace("'", '').replace('"', '')
    return out
