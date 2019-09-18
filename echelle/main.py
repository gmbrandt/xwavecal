"""
main.py: Main driver script for the minimal pipeline included with this package.

This is inspired by the Template Method Pattern and the Framework of Banzai by Curtis McCully.

"""
import logging as logger
from ast import literal_eval as to_list
from configparser import ConfigParser
import importlib

from echelle.utils.runtime_utils import parse_args, get_data_paths


def run():
    args = parse_args()
    config = ConfigParser()
    config.read(args.config_file)

    data_class = config.get('data', 'data_class', fallback='echelle.images.Image')
    extension = config.get('data', 'primary_data_extension')
    DataClass = importlib.import_module(data_class)

    data_paths = get_data_paths(args.input_dir, ['.fits'])
    stages_todo = [importlib.import_module(stage) for stage in to_list(config.get('reduction', 'stages'))]
    for data_path in data_paths:
        data = DataClass.load(data_path, )
        logger.info('Reducing {path} assuming a data class of {data_class} and raw data in extension {extension}'
                    ''.format(path=data_path, data_class=data_class, extension=extension))
        # need a function to get a specific attribute from a header based on that specified in config.ini
        # e.g. data.type = 'lampflat' or 'wavecal'
        for stage in stages_todo:
            if should_do(stage, data):
                data = stage.do_stage(data)
                # consider feeding in the logger as do_stage(data, logger)
        # data.update_filepath(base_dir=args.output_dir)
        logger.info('Writing output to {path}'.format(path=data.filepath))
        data.write(fpack=args.fpack)
