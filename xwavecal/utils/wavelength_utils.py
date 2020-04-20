import numpy as np

from xwavecal.utils.misc_utils import find_peaks, find_nearest
from astropy.stats import sigma_clip


class Model(dict):
    def add_polynomial_term(self, x_degree, i_degree):
        """
        :param self: dict
        :param x_degree: int, the degree of the x term to add to the polynomial model
        :param i_degree: int, the degree of the order term to add to the polynomial model
        :return: model such that model[x_degree] is a list with i_degree in the list.
        """
        if x_degree not in self:
            self[x_degree] = []
        if i_degree not in self[x_degree]:
            self[x_degree].append(i_degree)

    def is_missing_polynomial_term(self, x_degree, i_degree):
        if x_degree not in self or i_degree not in self[x_degree]:
            return True
        return False


def normalize_coordinates(coordinates, max_value, min_value):
    """
    :param coordinates: list or ndarray
    :param max_value: max_value coordinates can achieve. e.g. if normalizing pixels, we would have 4096 as max_value
    :return: coordinates normalized to run from -1 to 1.
    """
    coordinates = np.array(coordinates)
    if min_value > max_value:
        raise ValueError('min_value > max_value')
    coordinates = 2. * (coordinates - min_value)/(max_value - min_value) - 1.
    return coordinates


def estimate_global_scale(detector_range, n, m0):
    """
    :param detector_range: max wavelength minus min wavelength of the detector in angstroms, e.g. 5000
    :param n: the number of distinct diffraction orders on the detector, e.g. 67
    :param m0: the principle order number, e.g. 52
    :return: The global scale constant for the wavelength solution, in units of Angstrom.
    """
    return detector_range * (1/m0 - 1/(m0 + n))**(-1)


def identify_lines(spectrum, stderr, min_snr=5, min_ptp=5, order_key='ref_id'):
    """
    :param spectrum: Table or dictionary with 'flux', 'pixel' and some order column
    :param stderr: float, array with shape (len(spectrum['flux']), 1),
                   or array with shape spectrum['flux'].shape : specifying the standard error globally,
                   per order or per pixel respectively.
    :param min_snr: min signal to noise ratio that a peak must possess.
    :param min_ptp: min peak-to-peak spacing in pixels
    :param order_key: the dict key such that spectrum[order_key] gives a column of order id's.
    :return: dictionary
             with keys 'pixel', 'order' and 'flux' giving 1d arrays of the peak (emission lines)
             locations pixel, order and their peak flux.
    """
    lines = {'order': [], 'pixel': [], 'pixel_err': [], 'flux': []}
    std_err = np.ones((len(spectrum), 1)) * stderr
    for row in range(len(spectrum)):
        peak_coordinates, peak_err, peak_indices = find_peaks(spectrum['flux'][row], spectrum['pixel'][row],
                                                              yerr=std_err[row], height=std_err[row] * min_snr,
                                                              distance=min_ptp,
                                                              prominence=0.5 * np.abs(spectrum['flux'][row]),
                                                              peak_width=2, window=6)
        if len(peak_indices) > 0:
            lines['flux'].extend(list(spectrum['flux'][row][peak_indices]))
            lines['pixel'].extend(list(peak_coordinates))
            lines['pixel_err'].extend(list(peak_err))
            lines['order'].extend([spectrum[order_key][row]]*len(peak_coordinates))
    # casting as arrays for convenient use.
    lines['pixel'] = np.array(lines['pixel'], dtype=np.float32)
    lines['pixel_err'] = np.array(lines['pixel_err'], dtype=np.float32)
    lines['order'] = np.array(lines['order'], dtype=np.int)
    lines['flux'] = np.array(lines['flux'], dtype=np.float32)
    return lines


def calc_residuals(wavelengths, reference_wavelengths):
    closest_ref_wavelengths = find_nearest(wavelengths, array_b=reference_wavelengths)
    return wavelengths - closest_ref_wavelengths


def restrict(lines, low_order=-np.inf, high_order=np.inf):
    if low_order == -np.inf and high_order == np.inf:
        return lines
    order_mask = np.logical_and(lines['order'] <= high_order, lines['order'] >= low_order)
    restricted_lines = {key: coordinates[order_mask] for key, coordinates in lines.items()}
    return restricted_lines


def _sigma_clip(residuals, lines, sigma=3, stdfunc=np.std):
    if len(residuals) == 0:
        return lines
    clip_mask = ~ sigma_clip(residuals, sigma=sigma, stdfunc=stdfunc, masked=True).mask
    lines = {key: coordinates[clip_mask] for key, coordinates in lines.items()}
    return lines


def pixel_order_as_array(spectrum):
    # TODO check if you NEED to call .data here, I don't think you do.
    pixels = spectrum['pixel'].data
    orders = np.repeat(spectrum['ref_id'].data.reshape(-1, 1), pixels.shape[1], axis=1)
    return pixels, orders
