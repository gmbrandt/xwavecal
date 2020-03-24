from datetime import datetime
import numpy as np
from ast import literal_eval
from configparser import ConfigParser
from scipy.interpolate import UnivariateSpline

from xwavecal.utils.fiber_utils import lit_fibers, lit_wavecal_fibers
from xwavecal.utils.trace_utils import SingleTraceFitter


class FakeContext(object):
    # parses the test_config.ini (which is just the nres_config.ini, but a protected copy for tests only).
    def __init__(self):
        config = ConfigParser()
        config.read('xwavecal/tests/data/test_config.ini')
        dictionary = {key: literal_eval(item) for key, item in config.items('reduction')}
        for attribute, value in dictionary.items():
            setattr(self, attribute, value)


class FakeImage(object):
    def __init__(self, nx=102, ny=100, overscan_size=2, ccdsum='2 2', epoch='20180807'):
        self.nx = nx
        self.ny = ny
        self.telescope_id = -1
        self.site = 'lsc'
        self.instrument = 'nres01'
        self.camera = None
        self.ivar = None
        self.ccdsum = ccdsum
        self.epoch = epoch
        self.overscan_size = overscan_size
        self.data = np.ones((ny, nx))
        self.filename = 'test.fits'
        self.filter = 'U'
        self.dateobs = datetime(2018, 8, 7)
        self.header = {'read_noise': 11, 'gain': 1.0, 'OBSTYPE': 'LAMPFLAT',
                       'rdnoise': 11, 'type': 'lampflat',
                       'observation_date': '2019-04-10T12:56:44.466',
                       'instrument': 'nres03', 'site_name': 'test', 'unique_id': 77,
                       'fiber_state': 'none&tung&tung', 'instrument2': 'fa13'}
        self.filepath = 'None'
        self.caltype = ''
        self.bpm = np.zeros((ny, nx-overscan_size), dtype=np.uint8)
        self.request_number = '0000331403'
        self.block_id = '254478983'
        self.molecule_id = '544562351'
        self.exptime = 30.0
        self.obstype = 'TEST'
        self.data_tables = {}
        self.trace = None
        self.fiber0_lit, self.fiber1_lit, self.fiber2_lit = 0, 1, 1
        self.fiber0_wavecal, self.fiber1_wavecal, self.fiber2_wavecal = 0, 1, 1
        self.wavelength_solution = {}
        self.translator = None

    def get_header_val(self, key):
        return self.header[key]

    def set_header_val(self, key, value):
        self.header[key] = value

    def num_lit_fibers(self):
        return len(lit_fibers(self))

    def num_wavecal_fibers(self):
        return len(lit_wavecal_fibers(self))

    @classmethod
    def load(cls, *args):
        return FakeImage()

    def write(self, *args, **kwargs):
        pass


def gaussian(x, A, b, sigma):
    return A * np.exp(-(x - b) ** 2 / (2 * sigma ** 2))


def noisify_image(image):
    """
    :param image: Banzai_nres FakeImage object.
    This adds poisson and read noise to an image with traces already on it, in that order.
    """
    image.data = np.random.poisson(image.data) + np.random.normal(0, image.get_header_val('RDNOISE'), image.data.shape)


def noisify_spectrum(spectrum, rdnoise=10):
    for i in range(len(spectrum)):
        spectrum['flux'][i] = np.random.poisson(spectrum['flux'][i]) + \
                                np.random.normal(0, rdnoise, spectrum['flux'][i].shape)
    return spectrum


def array_with_peaks(x, centroids, amplitudes, stds):
    vectorized_gaussian = np.vectorize(gaussian)
    y = np.zeros_like(x)
    for centroid, amplitude, std in zip(centroids, amplitudes, stds):
        y += vectorized_gaussian(x, amplitude, centroid, std)
    return y


def fill_image_with_traces(image, poly_order_of_traces=4, order_width=1.5, fiber_intensity=1E4, max_num_traces=1000):
    """
    :param image: Banzai image object, where image.data is a 2d array of the image data
    :param poly_order_of_traces: the max order of the polynomial describing the traces. Maximum of 4.
    :param order_width: width of the traces in pixels
    :param fiber_intensity: peak intensity of the orders
    :param max_num_traces: max number of traces to try and fit onto the image
    :return: An image populated with semi-realistic traces which bend like parabolas, but are of degree
    min(4, poly_order_of_traces).
    """
    trace_fitter = SingleTraceFitter(image_data=image.data,
                                     second_order_coefficient_guess=0,
                                     poly_fit_order=4)
    num_fake_traces = min(int((image.data.shape[1] - 60)/20), max_num_traces)
    coefficients = np.zeros((num_fake_traces, 4+1))
    coefficients[:, 0] = np.linspace(30, image.data.shape[0] - 30, num=num_fake_traces)
    if poly_order_of_traces >= 2:
        coefficients[:, 2] = np.linspace(30, 40, num=num_fake_traces)
    if poly_order_of_traces >= 3:
        coefficients[:, 3] = np.linspace(1, 3, num=num_fake_traces)
    if poly_order_of_traces >= 4:
        coefficients[:, 4] = np.linspace(5, 10, num=num_fake_traces)

    trace_centers = trace_fitter._centers_from_coefficients(coefficients)
    trace_overlay = np.zeros_like(image.data).astype(np.float64)
    vectorized_gaussian = np.vectorize(gaussian)
    for x_pixel in range(trace_centers.shape[1]):
        for i in range(num_fake_traces):
            centroid = trace_centers[i, x_pixel]
            low, high = max(0, int(centroid - 5 * order_width)), min(trace_centers.shape[1] - 1,
                                                                     int(centroid + 5 * order_width)) + 1
            evalwindow = np.arange(low, high, 1)
            if len(evalwindow) > 0:
                trace_overlay[low: high, x_pixel] += vectorized_gaussian(evalwindow, 1, centroid, order_width)
    image.data += trace_overlay * fiber_intensity
    second_order_coefficient_guess = np.mean(coefficients[:, 2])
    return image, trace_centers, second_order_coefficient_guess


class SpectrumUtils:
    @staticmethod
    def _restrict(line_lambdas, lambdas):
        return np.logical_and(line_lambdas >= np.min(lambdas), line_lambdas <= np.max(lambdas))

    @staticmethod
    def uniform_min_space(low, high, size, minspace=0.4):
        out = np.sort(np.random.uniform(low, high, 5*size))
        for i in range(5):
            out = out[1:][~np.isclose(out[1:], out[:-1], atol=minspace)]
        np.random.shuffle(out)
        return out[:size]

    def populate_orders(self, line_lambdas, line_flux, wcs):
        measured_lines = {'pixel': [], 'order': [], 'corrected_flux': [], 'true_wavelength': []}
        pixel = np.arange(wcs.min_pixel, wcs.max_pixel + 1)
        for order in np.arange(wcs.min_order, wcs.max_order + 1):
            order_coords = np.ones_like(pixel) * order
            lam_to_x = UnivariateSpline(wcs(pixel, order_coords), pixel, s=0, k=3)
            in_order = self._restrict(line_lambdas, wcs(pixel, order_coords))
            lines_x_pos = lam_to_x(line_lambdas[in_order])
            measured_lines['pixel'].extend(list(lines_x_pos))
            measured_lines['order'].extend(list(np.ones_like(lines_x_pos) * order))
            measured_lines['corrected_flux'].extend(list(line_flux[in_order]))
            measured_lines['true_wavelength'].extend(list(line_lambdas[in_order]))
        return measured_lines

    def generate_measured_lines(self, n_list, n_true, n_garb, wcs, lam_range=(3800, 9000), min_amp=500, max_amp=1E5):
        """
        :param n_list: number of lines in the line list
        :param n_true: number of lines which are in the spectrum but not in the line list
        :param n_garb: number of lines which are not true emission lines. These will
        be injected randomly per order. E.g. cosmic rays.
        :param wcs: wavelength solution to follow to generate the line positions
        :return: Dictionary: measured_lines with n_list + n_true + n_garb lines.
                 measured_lines['pixel'] is the pixel position of the line centroid.
                 measured_lines['order'] is the order index of the line.
                 measured_lines['corrected_flux'] is the flux of the line, where duplicated lines
                 have the same flux.
                 measured_lines['true_wavelength'] is the true wavelength of the line. For garbage
                 lines, this is np.nan.
        """
        # line lambda locations
        full_list = self.uniform_min_space(min(lam_range), max(lam_range), size=n_list + n_true)
        line_list, line_other = full_list[:n_list], full_list[:n_true]
        # fluxes
        line_list_flux = np.random.randint(int(min_amp), int(max_amp), size=n_list)
        line_other_flux = np.random.randint(int(min_amp), int(max_amp), size=n_true)
        measured_lines = self.populate_orders(np.hstack([line_list, line_other]),
                                              np.hstack([line_list_flux, line_other_flux]), wcs)
        # adding garbage lines
        line_garbage_x = np.random.randint(wcs.min_pixel, wcs.max_pixel, size=n_garb).astype(np.float32)
        line_garbage_o = np.random.randint(wcs.min_order, wcs.max_order+1, size=n_garb)
        line_garbage_flux = np.random.randint(int(min_amp), int(max_amp), size=n_garb)
        measured_lines['pixel'].extend(list(line_garbage_x))
        measured_lines['order'].extend(list(line_garbage_o))
        measured_lines['corrected_flux'].extend(list(line_garbage_flux))
        measured_lines['true_wavelength'].extend(list(np.zeros(n_garb, dtype=float) * np.nan))
        # cast the elements of measured lines into arrays
        for key, item in measured_lines.items():
            measured_lines[key] = np.array(item)
        return measured_lines, line_list

