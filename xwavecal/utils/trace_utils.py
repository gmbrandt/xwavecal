"""
trace_utils.py: Routines for finding xwavecal orders across a CCD.
Authors
    G. Mirek Brandt (gmbrandt@ucsb.edu)

    Tim Brandt (tbrandt@physics.ucsb.edu)

    Curtis McCully (cmccully@lcogt.net)
"""


import numpy as np
from scipy import ndimage, optimize, signal
from astropy.table import Table, Column
from xwavecal.images import DataProduct

import logging
logger = logging.getLogger(__name__)


class Trace(DataProduct):
    """
    :param data = {'id': ndarray, 'centers': ndarray}. 'centers' gives a 2d array, where
    the jth row are the y centers across the detector for the trace with identification trace_centers['id'][j]
    """
    def __init__(self, data=None, table_name=None, num_centers_per_trace=0, filepath=None,
                 header=None, translator=None):
        super(Trace, self).__init__(data=data, data_name=table_name, filepath=filepath,
                                    header=header, translator=translator)
        if data is None and num_centers_per_trace <= 0:
            raise ValueError('Trace object instantiated but no trace data given and num_centers_per_trace is not > 0')
        if data is None:
            data = Table([Column(name='id'), Column(name='centers', shape=(num_centers_per_trace,))])
            data['id'].description = 'Identification tag for trace'
            data['centers'].description = 'Vertical position of the center of the' \
                                          ' trace as a function of horizontal pixel'
            data['centers'].unit = 'pixel'
        self.data = Table(data)

    def get_centers(self, row):
        return self.data['centers'][row]

    def get_id(self, row):
        return self.data['id'][row]

    def add_centers(self, trace_centers, trace_id):
        self.data.add_row([trace_id, trace_centers])

    def del_centers(self, rows):
        self.data.remove_rows(rows)

    def num_traces_found(self):
        return len(self.data['id'])

    def sort(self):
        center = int(self.data['centers'].shape[1] / 2)
        self.data['centers'] = np.array(self.data['centers'])
        self.data['centers'] = self.data['centers'][self.data['centers'][:, center].argsort()]
        self.data['id'] = np.arange(self.data['centers'].shape[0])

    def remove_duplicates(self, thresh=3):
        c = int(self.data['centers'].shape[1] / 2)
        duplicates = np.arange(len(self.data['centers']) - 1)[np.isclose(self.data['centers'][1:, c],
                                                                         self.data['centers'][:-1, c], atol=thresh)]
        self.del_centers(duplicates)
        self.sort()


class AllTraceFitter(object):
    def __init__(self, xmin=2000, xmax=2100, min_peak_to_peak_spacing=10, min_snr=6):
        self.xmin = xmin
        self.xmax = xmax
        self.min_peak_to_peak_spacing = min_peak_to_peak_spacing
        self.min_snr = min_snr

    def fit_traces(self, trace, image_data, poly_fit_order, second_order_coefficient_guess, image_noise_estimate=10):
        """
        :param trace: Trace object with get_centers, num_traces_found, add_centers, del_centers, and sort method
        :param image_data: The image data array, e.g. image.data or hdu[1].data for a fits HDUlist
        :param poly_fit_order: degree of the polynomial which will fit the xwavecal orders across the detector
        :param second_order_coefficient_guess: guess for the 2nd order coefficient of the aforementioned polynomial
        :param image_noise_estimate: estimate of the images read noise.
        :return: trace object with the y coordinates of the trace centers as a function of x pixel.
        """
        trace_fitter = SingleTraceFitter(image_data=image_data,
                                         second_order_coefficient_guess=second_order_coefficient_guess,
                                         poly_fit_order=poly_fit_order)
        peak_xy_coordinates = self._identify_traces(image_data, image_noise_estimate)
        trace = self._step_through_detector(trace, trace_fitter, peak_xy_coordinates)
        trace.sort()
        trace.remove_duplicates(thresh=self.min_peak_to_peak_spacing/2)
        return trace

    @staticmethod
    def _step_through_detector(trace, trace_fitter, peak_xy_coordinates):
        for trace_id, peak_xy_coordinate in enumerate(peak_xy_coordinates):
            trace_fitter.update_initial_guess_to_run_through_pt(peak_xy_coordinate)
            trace_centers = trace_fitter.fit_trace()
            trace.add_centers(trace_centers, trace_id)
            trace_fitter.use_fit_as_initial_guess(-1)
        return trace

    def _identify_traces(self, bkg_subtracted_image_data, image_noise_estimate):
        """
        :param bkg_subtracted_image_data: Image with the background subtracted
        :param image_noise_estimate: Read noise on the image.
        :return:
        """
        flux_down_detector = self._filtered_flux_down_detector(bkg_subtracted_image_data, image_noise_estimate)
        peak_y_coordinates = signal.find_peaks(flux_down_detector, distance=self.min_peak_to_peak_spacing,
                                               height=self.min_snr)[0]
        peak_x_coordinates = np.ones_like(peak_y_coordinates) * np.mean((self.xmin, self.xmax), dtype=int)
        peak_xy_coordinates = list(zip(peak_x_coordinates, peak_y_coordinates))
        return peak_xy_coordinates

    def _filtered_flux_down_detector(self, image_data, image_noise_estimate):
        wedge_of_image = image_data[:, self.xmin:self.xmax]
        noise_normalized_wedge = wedge_of_image / np.sqrt(image_noise_estimate**2 + np.abs(wedge_of_image))
        median_slice = np.median(noise_normalized_wedge, axis=1)
        return median_slice


class SingleTraceFitter(object):
    def __init__(self, image_data=None, poly_fit_order=2, second_order_coefficient_guess=None,
                 extraargs=None):
        if extraargs is None:
            extraargs = {}
        if extraargs.get('coefficients') is None:
            extraargs['coefficients'] = []
        self.second_order_coefficient_guess = second_order_coefficient_guess
        self.image_data = image_data
        self.poly_fit_order = poly_fit_order
        self.filtered_image_data = extraargs.get('filtered_image_data')
        self.initial_guess_next_fit = extraargs.get('initial_guess_next_fit')
        self.coefficients = extraargs.get('coefficients')
        self.x = extraargs.get('x')
        self.x_norm = extraargs.get('xnorm')
        self.design_matrix = extraargs.get('design_matrix')

        if extraargs.get('initialize_fit_objects', True):
            if self.second_order_coefficient_guess is None:
                logger.error('The second order guess have been specified '
                             ', aborting trace fitting.')
                raise ValueError('Starting y position up the detector nor the second order guess have been specified')
            self._normalize_domain_coordinates()
            self.filtered_image_data = self._prefilter_image_data(self.image_data)
            self.design_matrix = self._generate_design_matrix(self.x_norm, self.poly_fit_order)
            self.initial_guess_next_fit = np.zeros(self.poly_fit_order + 1).astype(np.float64)
            self.initial_guess_next_fit[2] = self.second_order_coefficient_guess

    def fit_trace(self):
        refined_coefficients = optimize.minimize(SingleTraceFitter._trace_merit_function, self.initial_guess_next_fit,
                                                 args=(self,), method='Powell').x
        self.coefficients.append(refined_coefficients)
        return self._centers_from_coefficients(refined_coefficients)

    def use_fit_as_initial_guess(self, fit_id):
        self.initial_guess_next_fit = np.copy(self.coefficients[fit_id])

    @staticmethod
    def _generate_design_matrix(normalized_domain, poly_fit_order):
        design_matrix = np.ones((poly_fit_order + 1, normalized_domain.shape[0]))
        for i in range(poly_fit_order + 1):
            design_matrix[i] = legendre(order=i, values=normalized_domain)
        return design_matrix

    def _centers_from_coefficients(self, coefficients):
        """
        :param coefficients: a 1d array for the coefficients of a single trace in order of: [0th, 1st,.. mth degree coefficient]
        or a 2d array of shape (N, m+1) for N traces.
        :return: If coefficients is 1d of size m+1, this will return an array of shape image_data.shape[1] with the trace
        centers at each of those points. If coefficients is for many traces (shape (N, m+1)) then this will return a 2d
        array of shape (N, image_data.shape[1]) e.g. the trace centers for each of the traces.
        """
        trace_centers = np.dot(coefficients, self.design_matrix)
        return trace_centers

    def update_initial_guess_to_run_through_pt(self, xy):
        """
        :param xy: (x,y) point on the detector you want the initial_guess for the next fit to pass through. (0,0)
                    is the smallest coordinate as usual.
                    This method does not change the shape of the trace generated by initial_guess_next_fit, it only
                    shifts the trace so that at the pixel x, the trace center has a y value of y.
        """
        # NOTE: This method calls all the trace centers, which means it evalutes 4096 centers, when it only needs one.
        current_centers = self._centers_from_coefficients(self.initial_guess_next_fit)
        self.initial_guess_next_fit[0] += xy[1] - current_centers[xy[0]]

    def _flux_across_trace(self, trace_centers):
        total_flux = np.sum(ndimage.map_coordinates(self.filtered_image_data, [trace_centers, self.x], prefilter=False))
        return total_flux

    def _normalize_domain_coordinates(self):
        self.x = np.arange(self.image_data.shape[1])
        # note this normalization scheme only works if np.min(x) = 0, which it is in this case.
        self.x_norm = self.x * 2. / self.x[-1] - 1  # x normalized to run from -1 to 1

    @staticmethod
    def _trace_merit_function(single_trace_coefficients, cls):
        trace_centers = cls._centers_from_coefficients(single_trace_coefficients)
        flux = cls._flux_across_trace(trace_centers)
        return (-1)*flux

    @staticmethod
    def _prefilter_image_data(image_data):
        return ndimage.spline_filter(image_data)


def legendre(order, values):
    #TODO replace with scipy.legendre(order)(values) if it works the same.
    if order == 0:
        return np.ones_like(values)
    else:
        return np.polynomial.legendre.legval(values, [0 for j in range(order)] + [1])
