"""
main.py: Main driver script for the minimal pipeline included with this package.

Author: G. Mirek Brandt

Note: This is built closely according to the structure of Banzai (https://github.com/lcogt/banzai) by Curtis McCully and
Las Cumbres Observatory.

"""
import logging as logger
from configparser import ConfigParser
from datetime import datetime
from ast import literal_eval
import os

from echelle.utils.runtime_utils import parse_args, get_data_paths, order_data, select_data_of_type, import_class, safe_eval
from echelle.utils.fits_utils import Translator
from echelle.database import format_db_info, add_data_to_db


class RuntimeContext(object):
    def __init__(self, dictionary):
        for attribute, value in dictionary.items():
            setattr(self, attribute, value)


def reduce_data(data_paths=None, args=None, config=None):
    logger.basicConfig(level=logger.DEBUG)
    if args is None:
        args = parse_args()
    if config is None:
        config = ConfigParser()
        config.read(args.config_file)
    if data_paths is None:
        data_paths = args.data_paths

    runtime_context, data_class, extension, header_keys, type_keys = organize_config(config)
    translator = Translator(header_keys, type_keys)
    DataClass = import_class(data_class)

    for data_path in data_paths:
        logger.info('Reducing {path} assuming a data class of {data_class} and raw data in extension {extension}'
                    ''.format(path=data_path, data_class=data_class, extension=extension))

        data = DataClass.load(data_path, extension, translator)
        stages_todo = [import_class(stage) for stage in literal_eval(config.get('stages', data.get_header_val('type')))]
        auxiliary_products = []
        for stage in stages_todo:
            data = stage(runtime_context).do_stage(data)
            if isinstance(data, list):
                auxiliary_products.extend(data[1:])
                data = data[0]

        to_write = [data] + auxiliary_products

        for data in to_write:
            data.filepath = make_output_path(args.output_dir, data)
            data.write(fpack=args.fpack)
            logger.info('Adding file to processed image database at {path}'.format(path=runtime_context.database_path))
            db_info = format_db_info(data, runtime_context.time_format)
            add_data_to_db(runtime_context.database_path, db_info)


def run():
    logger.basicConfig(level=logger.DEBUG)
    # parse command line arguments and the configuration file.
    args = parse_args()
    config = ConfigParser()
    config.read(args.config_file)

    # get the data paths of the data to reduce.
    runtime_context, data_class, extension, header_keys, type_keys = organize_config(config)
    DataClass = import_class(data_class)
    data_paths = select_data(args.input_dir, args.frame_type, literal_eval(config.get('data', 'files_contain')),
                             DataClass, extension, header_keys, type_keys)

    logger.info('Found {0} files of {1} type'.format(len(data_paths), args.frame_type))
    reduce_data(data_paths, args, config)


def organize_config(config):
    # pull information from the configuration file.
    runtime_context = RuntimeContext({key: safe_eval(item) for key, item in config.items('reduction')})
    data_class = config.get('data', 'data_class', fallback='echelle.images.Image')
    extension = config.getint('data', 'primary_data_extension')
    header_keys = literal_eval(config.get('data', 'header_keys'))
    type_keys = literal_eval(config.get('data', 'type_keys'))
    return runtime_context, data_class, extension, header_keys, type_keys


def select_data(input_dir, frame_type, files_contain, data_class, extension, header_keys, type_keys):
    data_paths = get_data_paths(input_dir, files_contain)
    data_paths = order_data(data_paths, data_class, extension, header_keys, type_keys)
    data_paths = select_data_of_type(data_paths, data_class, extension, header_keys, type_keys, frame_type)
    return data_paths


def make_output_path(output_dir, data, time_fmt='%Y-%m-%dT%H:%M:%S.%f'):
    """
    :param data: Image
    :return: string
    """
    id = str(data.get_header_val('unique_id')).zfill(4)
    dayobs = datetime.strptime(data.get_header_val('observation_date'), time_fmt).strftime('%Y%m%d')
    filename = '{site}_{inst}_{dayobs}_{id}_{type}_{f0}{f1}{f2}.fits'.format(inst=data.get_header_val('instrument'),
                                                                         site=data.get_header_val('site_name'),
                                                                         dayobs=dayobs, id=id,
                                                                         type=data.get_header_val('type'),
                                                                         f0=data.fiber0_lit, f1=data.fiber1_lit,
                                                                         f2=data.fiber2_lit)
    return os.path.join(output_dir, filename)
