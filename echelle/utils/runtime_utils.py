import datetime

from banzai_nres.utils.db_utils import get_raw_path

import logging




def get_telescope_filename(image):
    return image.header.get('TELESCOP', '').replace('nres', 'nrs')


def validate_raw_path(runtime_context, raw_path):
    if raw_path is not None and 'raw' not in raw_path.lower():
        raw_path = get_raw_path(base_raw_path=raw_path, runtime_context=runtime_context)
    return raw_path


def get_frame_types(runtime_context, default_frames_to_reduce):
    if getattr(runtime_context, 'frame_type', None) is None:
        frame_types = default_frames_to_reduce
    else:
        frame_types = [runtime_context.frame_type]
    return frame_types


def get_reduction_date_window(runtime_context):
    max_date = getattr(runtime_context, 'max_date', None)
    min_date = getattr(runtime_context, 'min_date', None)
    if max_date is None:
        max_date = datetime.datetime.utcnow()

    if min_date is None:
        min_date = max_date - datetime.timedelta(hours=24)

    if min_date > max_date:
        logger.error('The start cannot be after the end. Aborting reduction!')
        raise ValueError('min_date > max_date.')
    return min_date, max_date
