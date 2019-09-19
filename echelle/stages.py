import abc
import logging as logger


class Stage(object):
    def __init__(self, runtime_context):
        self.runtime_context = runtime_context

    @abc.abstractmethod
    def do_stage(self, image):
        return image


class ApplyCalibration(Stage):
    def __init__(self, runtime_context):
        super(ApplyCalibration, self).__init__(runtime_context)

    def do_stage(self, image):
        master_calibration_path = self.get_calibration_filename(image)
        if master_calibration_path is None:
            self.on_missing_master_calibration(image)
            return image
        return self.apply_master_calibration(image, master_calibration_path)

    def on_missing_master_calibration(self, image):
        logger.error('Master calibration file not found.')

    def get_calibration_filename(self, image):
        return ''

    @abc.abstractmethod
    def apply_master_calibration(self, image, path):
        return image