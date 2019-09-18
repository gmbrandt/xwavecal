import datetime
import logging as logger
import argparse
import glob
import os


def parse_args():
    parser = argparse.ArgumentParser(description='Reduce an echelle spectrograph frame.')
    parser.add_argument("--output-dir", required=True,
                        help="Directory within which to save the processed data files.")
    parser.add_argument("--input-dir", required=True,
                        help="Directory which contains the raw data files.")
    parser.add_argument("--fpack", required=False, action='store_true',
                        help="fpack output files with default quantization of 64")
    parser.add_argument("--config-file", required=True,
                        help="Path to the instrument specific configuration file.")
    parser.add_argument("--frame-type", required=False, default='all',
                        help="Frame type to either fit traces to or wavelength calibrate."
                             "Make sure frame type settings are appropriately set in the config file."
                             "lampflat files are used for tracing, wavecals are wavelength calibration"
                             "frames such as ThAr exposures.",
                        choices=['lampflat', 'wavecal', 'all'], type=str.lower)
    args = parser.parse_args()
    return args


def get_data_paths(dir_path, files_contain=['.fits']):
    all_files = glob.glob(os.path.join(dir_path, '*'))
    return [file for file in all_files if all([item in file for item in files_contain])]
