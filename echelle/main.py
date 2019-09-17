"""
main.py: Main driver script for the banzai-NRES pipeline.
    reduce_night() is the console entry point.
Authors
    Curtis McCully (cmccully@lcogt.net)

    G. Mirek Brandt (gmbrandt@ucsb.edu)
"""
import datetime

import banzai_nres.settings as nres_settings  # import to override banzai settings.
import banzai.settings as settings
from banzai.main import process_directory, parse_directory_args
from banzai.utils import date_utils
from banzai import dbs
from banzai import logs
from banzai.main import run_master_maker

from banzai_nres.utils.runtime_utils import validate_raw_path, get_frame_types, get_reduction_date_window

import logging as logger




def process_master_maker(runtime_context, instrument, frame_type_to_stack, min_date, max_date,
                         master_frame_type=None, use_masters=False):
    if master_frame_type is None:
        master_frame_type = frame_type_to_stack
    extra_tags = {'instrument': instrument.camera, 'master_frame_type': master_frame_type,
                  'min_date': min_date.strftime(date_utils.TIMESTAMP_FORMAT),
                  'max_date': max_date.strftime(date_utils.TIMESTAMP_FORMAT)}
    logger.info("Making master frames", extra_tags=extra_tags)
    image_path_list = dbs.get_individual_calibration_images(instrument, frame_type_to_stack, min_date, max_date,
                                                            use_masters=use_masters,
                                                            db_address=runtime_context.db_address)
    if len(image_path_list) == 0:
        logger.info("No calibration frames found to stack", extra_tags=extra_tags)

    try:
        run_master_maker(image_path_list, runtime_context, master_frame_type)
    except Exception:
        logger.error(logs.format_exception())


def reduce_night(runtime_context=None, raw_path=None):
    all_frame_types = list(settings.LAST_STAGE.keys())
    extra_console_arguments = [{'args': ['--site'],
                                'kwargs': {'dest': 'site', 'help': 'Site code (e.g. ogg)', 'required': True}},
                               {'args': ['--camera'],
                                'kwargs': {'dest': 'camera', 'help': 'Camera (e.g. kb95)', 'required': True}},
                               {'args': ['--instrument-name'],
                                'kwargs': {'dest': 'instrument_name', 'help': 'Instrument (e.g. nres04)', 'required': True}},
                               {'args': ['--enclosure'],
                                'kwargs': {'dest': 'enclosure', 'help': 'Enclosure code (e.g. clma)',
                                           'required': False}},
                               {'args': ['--telescope'],
                                'kwargs': {'dest': 'telescope', 'help': 'Telescope code (e.g. 0m4a)',
                                           'required': False}},
                               {'args': ['--frame-type'],
                                'kwargs': {'dest': 'frame_type', 'help': 'Type of frames to process',
                                           'choices': all_frame_types, 'required': False}},
                               {'args': ['--min-date'],
                                'kwargs': {'dest': 'min_date', 'required': False, 'type': date_utils.valid_date,
                                           'help': 'Earliest observation time of the individual calibration frames. '
                                                   'Must be in the format "YYYY-MM-DDThh:mm:ss".'}},
                               {'args': ['--max-date'],
                                'kwargs': {'dest': 'max_date', 'required': False, 'type': date_utils.valid_date,
                                           'help': 'Latest observation time of the individual calibration frames. '
                                                   'Must be in the format "YYYY-MM-DDThh:mm:ss".'}}]

    runtime_context, raw_path = parse_directory_args(runtime_context, raw_path=raw_path,
                                                     extra_console_arguments=extra_console_arguments)
    instrument = dbs.query_for_instrument(runtime_context.db_address, runtime_context.site,
                                          camera=runtime_context.camera, name=runtime_context.instrument_name,
                                          enclosure=None, telescope=None)
    frame_types = get_frame_types(runtime_context, default_frames_to_reduce=all_frame_types)
    min_date, max_date = get_reduction_date_window(runtime_context)
    raw_path = validate_raw_path(runtime_context, raw_path)

    for frame_type in frame_types:
        if frame_type == 'TRACE':
            frame_type_to_stack = 'LAMPFLAT'
            use_masters = True
            master_frame_type = 'TRACE'
        else:
            frame_type_to_stack = frame_type
            use_masters = False
            master_frame_type = None
            # must reduce frames before making the master calibration, unless we are making a master trace.
            process_directory(runtime_context, raw_path, [frame_type_to_stack])

        if frame_type in settings.CALIBRATION_IMAGE_TYPES:
            process_master_maker(runtime_context, instrument, frame_type_to_stack.upper(),
                                 min_date=min_date, max_date=max_date, master_frame_type=master_frame_type,
                                 use_masters=use_masters)
