"""
traces.py: Module for finding xwavecal orders across a CCD.
"""
import os
from copy import deepcopy, copy

from xwavecal.utils.trace_utils import Trace, AllTraceFitter
from xwavecal.stages import Stage, ApplyCalibration


class TraceMaker(Stage):
    def __init__(self, runtime_context):
        super(TraceMaker, self).__init__(runtime_context)
        self.order_of_poly_fit = runtime_context.trace_fit_polynomial_order
        self.second_order_coefficient_guess = runtime_context.trace_fit_initial_degree_two_guess
        self.trace_table_name = runtime_context.trace_table_name
        self.xmin = runtime_context.window_for_trace_identification['min']
        self.xmax = runtime_context.window_for_trace_identification['max']
        self.min_peak_to_peak_spacing = runtime_context.min_fiber_to_fiber_spacing
        self.min_snr = runtime_context.min_snr_for_trace_identification

    @property
    def calibration_type(self):
        return 'TRACE'

    def do_stage(self, image):
        self.logger.info('fitting traces order by order', )
        fitter = AllTraceFitter(xmin=self.xmin, xmax=self.xmax,
                                min_peak_to_peak_spacing=self.min_peak_to_peak_spacing,
                                min_snr=self.min_snr)
        trace = Trace(data=None, filepath=copy(image.filepath), header=deepcopy(image.header), translator=image.translator,
                      num_centers_per_trace=image.data.shape[1], table_name=self.trace_table_name)
        trace = fitter.fit_traces(trace=trace, image_data=image.data,
                                  poly_fit_order=self.order_of_poly_fit,
                                  second_order_coefficient_guess=self.second_order_coefficient_guess,
                                  image_noise_estimate=image.get_header_val('read_noise'))
        trace.set_header_val('type', self.calibration_type.lower())
        trace.fiber0_lit, trace.fiber1_lit, trace.fiber2_lit = image.fiber0_lit, image.fiber1_lit, image.fiber2_lit
        self.logger.info('Created master trace')
        image.trace = trace
        return image, trace


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
        image.trace = Trace.load(master_calibration_path, extension_name=self.runtime_context.trace_table_name)
        master_trace_filename = os.path.basename(master_calibration_path)
        image.set_header_val('IDTRACE', (master_trace_filename, 'ID of trace centers file'))
        self.logger.info('Applied loaded trace centers to image')
        return image
