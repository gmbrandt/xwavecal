import abc
import os

from echelle.database import query_db_for_nearest

import logging
logger = logging.getLogger(__name__)


class Stage(object):
    def __init__(self, runtime_context):
        self.runtime_context = runtime_context

    @abc.abstractmethod
    def do_stage(self, image):
        return image  # pragma: no cover


class ApplyCalibration(Stage):
    def __init__(self, runtime_context):
        super(ApplyCalibration, self).__init__(runtime_context)

    @property
    def calibration_type(self):
        return 'None'  # pragma: no cover

    def do_stage(self, image):
        master_calibration_path = self.get_calibration_filename(image)
        if master_calibration_path is None:
            self.on_missing_master_calibration(image)
            return image
        elif not os.path.exists(master_calibration_path):
            self.on_missing_master_calibration(image)
            return image
        return self.apply_master_calibration(image, master_calibration_path)

    def on_missing_master_calibration(self, image):
        logger.error('Master calibration file not found.')

    def get_calibration_filename(self, image):
        return query_db_for_nearest(self.runtime_context.database_path,
                                    image, self.calibration_type.lower(),
                                    time_format=self.runtime_context.time_format)

    @abc.abstractmethod
    def apply_master_calibration(self, image, path):
        return image  # pragma: no cover
