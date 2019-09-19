import mock
import numpy as np
from scipy import interpolate
from astropy.table import Table

import echelle.utils.overlap_utils as overlapu
from echelle.utils.overlap_utils import OverlapFitter


class Utils:
    @staticmethod
    def fake_overlap(b_lines, r_lines, b_flux, r_flux, linearscale, *args, **kwargs):
        return {'pixel': r_lines, 'matched_pixel': b_lines, 'peaks': 1}

    def uniform_sample_stretched_signal(self, signal, coords, coefficients):
        signal_f, transformed_coords = self.stretch_signal(signal, coords, coefficients)
        uniformly_sampled_coords = np.arange(int(np.min(transformed_coords)), int(np.max(transformed_coords)) + 1)
        return signal_f(uniformly_sampled_coords), uniformly_sampled_coords

    @staticmethod
    def stretch_signal(signal, coords, coefficients, fill_value=0):
        transformed_coords = overlapu.coordinate_transform(coords, coefficients=coefficients)
        signal_f = interpolate.interp1d(x=transformed_coords,
                                        y=signal,
                                        bounds_error=False, kind='linear', fill_value=fill_value)
        return signal_f, transformed_coords


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

    @mock.patch('echelle.utils.overlap_utils.OverlapFitter._fit_overlap', return_value=[100, 1.4, 2E-5])
    def test_peaks_save(self, mock_fit):
        r_lines = np.random.randint(0, 100, 20)
        b_lines = overlapu.coordinate_transform(r_lines, [100, 1.4, 2E-5])
        overlap = OverlapFitter().fit(b_lines, r_lines, None, None, None, peaks_only=True)
        assert np.allclose(overlap['matched_pixel'], b_lines)
        assert np.allclose(overlap['pixel'], r_lines)
        assert np.isclose(overlap['peaks'], 20)

    @mock.patch('echelle.utils.overlap_utils.OverlapFitter.fit', side_effect=Utils.fake_overlap)
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

    @mock.patch('echelle.utils.overlap_utils._is_bad', side_effect=min)
    def test_flag_bad_overlaps(self, fake_overlap_check):
        overlaps = Table({'pixel': np.array([0, 1, 1, 0]), 'good': np.array([True, True, True, True])})
        overlaps = overlapu.flag_bad_overlaps(overlaps)
        assert np.allclose(overlaps['good'], [True, False, False, True])
