import numpy as np
import mock
from scipy import optimize

from xwavecal.tests.utils import array_with_peaks
from xwavecal.utils.misc_utils import find_peaks, brute_local_min, find_nearest,\
                                         overlap_region, normalize_by_brightest


def test_find_peaks():
    min_snr = 1
    peak_centers = np.array([15.2, 40.2, 60.2, 80])
    pixel = np.arange(100).astype(float)
    std = 1.5  # standard deviation of the gaussian peaks.
    flux = array_with_peaks(x=pixel, centroids=peak_centers,
                            amplitudes=[1000, 6, min_snr + 0.1, min_snr - 0.1], stds=[std, std, std, std])
    center_uncertainties = std/np.sqrt([np.sum(flux[np.isclose(pixel, cntr, atol=8*std)]) for cntr in peak_centers])
    for min_height in [min_snr, min_snr * np.ones_like(flux)]:
        found_peak_centers, errs = find_peaks(flux, pixel, np.sqrt(flux), height=min_height, distance=10, window=6)[:2]
        assert np.allclose(found_peak_centers, peak_centers[:-1], rtol=.01)
        assert np.allclose(errs, center_uncertainties[:-1], rtol=.1)


@mock.patch('xwavecal.utils.misc_utils.fit_peaks', return_value=np.array([[0.1, 0.01], [13, 2]]))
@mock.patch('scipy.signal.find_peaks', return_value=(np.array([0, 1]), {}))
def test_identify_lines(mock_peaks, mock_refine):
    peak_x, peak_err, peak_ind = find_peaks(None, np.array([0, 20]), None)
    assert np.allclose(peak_x, [0.1, 20])
    assert np.allclose(peak_err, [0.01, 1])


def test_brute_force_min():
    def dimple(x, dummy):
        return -x**2 - 50 * np.exp(-x ** 2/(2*0.3**2)) + 500
    root = optimize.minimize(dimple, x0=1.01, method=brute_local_min, args=None,
                             options={'step': 1/10, 'rrange': (-15, 15), 'filtw': 31, 'finish': None}).x
    assert np.isclose(root[0], 0, atol=0.1)
    root = optimize.minimize(dimple, x0=1.01, method=brute_local_min, args=None,
                             options={'step': 1/10, 'rrange': (-15, 15), 'filtw': 31, 'finish': 'Nelder-Mead'}).x
    assert np.isclose(root[0], 0)


class TestFindNearest:
    # the find_nearest function is at the heart of the wavelength solution and warrants exhaustive testing.
    @staticmethod
    def find_nearest_slow(array_a, array_b):
        nearest_value_idxs = [(np.abs(array_b - value)).argmin() for value in array_a]
        return array_b[nearest_value_idxs]

    def test_strenuous(self):
        b = np.sort(np.random.random(size=100))
        a = np.random.random(size=100) * 2 - 1
        assert np.allclose(self.find_nearest_slow(a, b), find_nearest(a, b))

    def test_simple(self):
        b = np.array([1, 2, 3])
        a = np.array([1.1, 1.9, 3.49])
        assert np.allclose(self.find_nearest_slow(a, b), b)
        assert np.allclose(find_nearest(a, b), b)


def test_overlap_region():
    assert np.allclose(overlap_region([1, 2], [1.5, 3]), (1.5, 2))
    assert overlap_region([1, 2], [2.001, 3]) is None


def test_normalize_by_brightest():
    spec = normalize_by_brightest(np.array([[1, 2, 3, 2], [1, 1, 1, 1]], dtype=float), n=3)
    assert np.allclose(spec, [[.5, 1, 3/2, 1], [1, 1, 1, 1]])
