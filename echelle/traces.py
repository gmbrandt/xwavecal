"""
traces.py: Driver scripts for finding echelle orders across a CCD.

Authors
    G. Mirek Brandt (gmbrandt@ucsb.edu)
"""
import sep
import os
import logging as logger

from echelle.utils.trace_utils import Trace, AllTraceFitter
from echelle.stages import Stage, ApplyCalibration

import echelle.settings as nres_settings


class TraceMaker(Stage):
    def __init__(self, runtime_context):
        super(TraceMaker, self).__init__(runtime_context)
        self.runtime_context = runtime_context
        self.order_of_poly_fit = nres_settings.TRACE_FIT_POLYNOMIAL_ORDER
        self.second_order_coefficient_guess = nres_settings.TRACE_FIT_INITIAL_DEGREE_TWO_GUESS
        self.trace_table_name = nres_settings.TRACE_TABLE_NAME
        self.xmin = nres_settings.WINDOW_FOR_TRACE_IDENTIFICATION['min']
        self.xmax = nres_settings.WINDOW_FOR_TRACE_IDENTIFICATION['max']
        self.min_peak_to_peak_spacing = nres_settings.MIN_FIBER_TO_FIBER_SPACING
        self.min_snr = nres_settings.MIN_SNR_FOR_TRACE_IDENTIFICATION

    @property
    def calibration_type(self):
        return 'TRACE'

    def do_stage(self, image):
        # need to set obstype to TRACE and set master filename.
        logger.info('fitting traces order by order', )
        bkg_subtracted_image_data = image.data - sep.Background(image.data).back()
        fitter = AllTraceFitter(xmin=self.xmin, xmax=self.xmax,
                                min_peak_to_peak_spacing=self.min_peak_to_peak_spacing,
                                min_snr=self.min_snr)
        trace = Trace(data=None, filepath=image.filepath, header=image.header,
                      num_centers_per_trace=image.data.shape[1], table_name=self.trace_table_name)
        trace = fitter.fit_traces(trace=trace, image_data=bkg_subtracted_image_data,
                                  poly_fit_order=self.order_of_poly_fit,
                                  second_order_coefficient_guess=self.second_order_coefficient_guess,
                                  image_noise_estimate=image.get_header_val('read_noise'))
        trace.set_header_val('type', self.calibration_type.lower())
        logger.info('Created master trace') # need to show information about the image.
        return trace


class LoadTrace(ApplyCalibration):
    """
    Loads trace coefficients from file and appends them onto the image object.
    """
    def __init__(self, runtime_context):
        super(LoadTrace, self).__init__(runtime_context)

    @property
    def calibration_type(self):
        return 'TRACE'

    def apply_master_calibration(self, image, master_calibration_path):
        image.trace = Trace.load(master_calibration_path, extension_name=nres_settings.TRACE_TABLE_NAME)
        master_trace_filename = os.path.basename(master_calibration_path)
        image.set_header_val('L1IDTRAC', (master_trace_filename, 'ID of trace centers file'))
        logger.info('Loading trace centers',   extra={'L1IDTRAC': image.get_header_val('L1IDTRAC')})
        return image
