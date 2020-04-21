import mock
import numpy as np
from astropy.table import Table

import xwavecal.utils.overlap_utils as overlapu
from xwavecal.utils.overlap_utils import OverlapFitter


class Utils:
    @staticmethod
    def fake_overlap(b_lines, r_lines, b_flux, r_flux, linearscale, *args, **kwargs):
        return {'pixel': r_lines, 'matched_pixel': b_lines, 'peaks': 1}


class TestOverlapFitter:
    def test_fit_overlap(self):
        true_coeffs = [100, 1.4, 2E-5]
        num_peaks = 10
        r_lines = np.random.uniform(0, 200, num_peaks)
        r_flux = np.random.randint(100, int(1E5), num_peaks)
        b_lines = overlapu.coordinate_transform(r_lines, true_coeffs)
        b_flux = 1. * r_flux
        coeffs = OverlapFitter()._fit_overlap(b_lines, r_lines, b_flux, r_flux, (0.5, 2), pixel_tol=2, flux_tol=0.2, deg=2)
        assert np.allclose(coeffs, true_coeffs, rtol=0.05)

    @mock.patch('xwavecal.utils.overlap_utils.OverlapFitter._fit_overlap', return_value=[100, 1.4, 2E-5])
    def test_peaks_save(self, mock_fit):
        r_lines = np.random.randint(0, 100, 20)
        b_lines = overlapu.coordinate_transform(r_lines, [100, 1.4, 2E-5])
        overlap = OverlapFitter().fit(b_lines, r_lines, None, None, None)
        assert np.allclose(overlap['matched_pixel'], b_lines)
        assert np.allclose(overlap['pixel'], r_lines)
        assert np.isclose(overlap['peaks'], 20)

    @mock.patch('xwavecal.utils.overlap_utils.OverlapFitter.fit', side_effect=Utils.fake_overlap)
    def test_fit_overlaps(self, mock_overlap):
        spectrum = {'flux': np.zeros((3, 3)), 'pixel': np.array([np.arange(3), np.arange(1, 4), np.arange(3)]),
                    'ref_id': np.arange(1, 4)}
        lines = {'corrected_flux': np.zeros(spectrum['pixel'].size), 'pixel': spectrum['pixel'].flatten(),
                 'order': (np.arange(1, 4).reshape(-1, 1) * np.ones_like(spectrum['pixel'])).flatten()}
        overlaps = overlapu.fit_overlaps(spectrum, lines, max_overlap_red=2, max_overlap_blue=2,
                                         linear_scale_range=(0.5, 2), deg=-1)
        for i in range(2):
            assert np.allclose(overlaps['pixel'][i][0], spectrum['pixel'][i][0])
            assert np.allclose(overlaps['matched_pixel'][i][0], spectrum['pixel'][i+1][-1])
            assert np.all(np.isnan(overlaps['pixel'][i][-1]))
            assert np.all(np.isnan(overlaps['matched_pixel'][i][-1]))
            assert np.isclose(overlaps['peaks'][i], 1)
            assert np.isclose(overlaps['ref_id'][i], i + 1)
            assert np.isclose(overlaps['matched_ref_id'][i], i + 2)

    def test_fit_overlaps_skips_on_no_lines(self):
        spectrum = {'flux': np.zeros((3, 3)), 'pixel': np.array([np.arange(3), np.arange(1, 4), np.arange(3)]),
                    'ref_id': np.arange(1, 4)}
        lines = {'corrected_flux': np.zeros(spectrum['pixel'].size), 'pixel': spectrum['pixel'].flatten(),
                 'order': (np.arange(1, 4).reshape(-1, 1) * np.ones_like(spectrum['pixel'])).flatten()}
        overlaps = overlapu.fit_overlaps(spectrum, lines, max_overlap_red=2, max_overlap_blue=2,
                                         linear_scale_range=(0.5, 2), deg=10)
        assert len(overlaps) == 0

    def test_flag_bad_overlaps(self):
        overlaps = Table({'peaks': np.array([[1], [0], [0], [1]]), 'good': [True, True, True, True]})
        overlaps = overlapu.flag_bad_overlaps(overlaps, 1)
        assert np.allclose(overlaps['good'], [True, False, False, True])
