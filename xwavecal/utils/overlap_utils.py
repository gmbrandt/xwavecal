import numpy as np
from scipy import interpolate
from astropy.table import Table, Column
from xwavecal.utils.misc_utils import find_nearest
import numpy.polynomial.polynomial as poly
import itertools


class OverlapFitter:
    def fit(self, b_lines, r_lines, b_lines_flux, r_lines_flux, linear_scale_range, flux_tol=0.2, deg=2):
        """
        :param deg. int
               deg >= 2.
        :return: (x, xp), a tuple of 1d ndarrays of equal length. For a given index i, x[i] maps to xp[i] in wavelength.
                 I.e. x and xp are the set of overlapping pixels in wavelength space for this pair of diffraction orders.
        """
        coeffs = self._fit_overlap(b_lines, r_lines, b_lines_flux, r_lines_flux,
                                   linear_scale_range, pixel_tol=2, flux_tol=flux_tol, deg=deg)
        coeffs = self._refine_fit(coeffs, b_lines, r_lines, pixel_tol=10, deg=deg)
        red_peaks = match_peaks(b_lines, r_lines, coeffs, pixel_tol=1)
        return {'pixel': red_peaks, 'matched_pixel': coordinate_transform(red_peaks, coeffs), 'peaks': len(red_peaks)}

    def _fit_overlap(self, b_lines, r_lines_f, b_lines_flux, r_lines_flux_f, linear_scale_range, pixel_tol=2,
                     flux_tol=0.2, deg=2):
        best_coeffs = np.array([0, 1, 0])

        r_lines, r_lines_flux = np.copy(r_lines_f), np.copy(r_lines_flux_f)

        #r_lines = r_lines[np.argsort(r_lines_flux)[::-1]][:15]
        #r_lines_flux = np.sort(r_lines_flux)[::-1][:15]
        r_lines_flux = r_lines_flux[np.argsort(r_lines)][:15]
        r_lines = np.sort(r_lines)[:15]
        b_lines_flux = b_lines_flux[np.argsort(b_lines)][-20:]
        b_lines = np.sort(b_lines)[-20:]
        # build nearest function for finding matching blue peaks.
        nearest_blue = interpolate.interp1d(b_lines, b_lines, kind='nearest', bounds_error=False,
                                            fill_value=(np.min(b_lines), np.max(b_lines)))
        matches, clipped_arr1 = self.match_matrix(r_lines, b_lines, r_lines_flux, b_lines_flux,
                                                  tol=flux_tol)
        vander = poly.polyvander(clipped_arr1, deg=deg)
        peak_combinations = itertools.combinations(np.arange(len(clipped_arr1)), deg+2)

        match_count = []
        all_coeffs = []
        for indices in peak_combinations:
            vander_inv = np.linalg.pinv(vander[list(indices)])
            lines2 = self._blue_matches(b_lines, matches, indices)
            potential_matched_blue_lines = self.line_combinations(lines2)
            if len(potential_matched_blue_lines) > 0:
                coeffs = np.matmul(vander_inv, potential_matched_blue_lines)
                coeffs[:, 2:, 0] = np.abs(coeffs[:, 2:, 0])  # ensures higher order coeffs are positive and thus
                # that the mappings are one-to-one.
                mapped_lines = np.matmul(vander, coeffs)
                num_matched_peaks = np.count_nonzero(np.isclose(mapped_lines, nearest_blue(mapped_lines),
                                                                atol=pixel_tol), axis=1).flatten()
                in_range = np.logical_and(coeffs[:, 1, 0] > min(linear_scale_range),
                                          coeffs[:, 1, 0] < max(linear_scale_range))
                coeffs, num_matched_peaks = coeffs[in_range], num_matched_peaks[in_range]
                if len(num_matched_peaks) > 0:
                    match_count.append(np.max(num_matched_peaks))
                    all_coeffs.append(coeffs[np.argmax(num_matched_peaks)])
        if len(match_count) > 0:
            best_coeffs = all_coeffs[np.argmax(match_count)]
        return np.array(best_coeffs).flatten()

    def _refine_fit(self, coeffs, b_lines, r_lines, pixel_tol=10, deg=2):
        """
        :param coeffs: array-like.
                       array of coefficients [a, b, c, ...] for g(x) = a + b*x + c*x^2 etc..
        :param b_lines: array-like
                        x-locations (pixel-locations) of emission lines from the blue order.
        :param r_lines: array-like
                        x-locations (pixel-locations) of emission lines from the red order.
        :param pixel_tol: tolerance for g(xr) - xb for the peaks at xr and xb to be considered
                          a matched-peak (i.e. duplicates of one-another). Where xr is the pixel
                          coordinate of the peak from the red order, and xb is the pixel coordinate
                          of the peak from the blue order.
        :param deg: polynomial degree of the mapping g(x)
        :return: ndarray
                 coefficients for g(x), [a, b, c,...] which minimize the chi^2 of (g(xr) - xb) for all
                 matched peaks (xr, xb) within the overlap.
        """
        peaks = poly.polyval(r_lines, coeffs)
        blue_peaks = interpolate.interp1d(b_lines, b_lines, kind='nearest', bounds_error=False,
                                          fill_value=(np.min(b_lines), np.max(b_lines)))(peaks)
        if np.count_nonzero(np.isclose(blue_peaks, peaks, atol=pixel_tol)) >= deg + 1:
            coeffs = poly.polyfit(r_lines[np.isclose(blue_peaks, peaks, atol=pixel_tol)],
                                  blue_peaks[np.isclose(blue_peaks, peaks, atol=pixel_tol)], deg)
        return coeffs

    @staticmethod
    def _blue_matches(blue_lines, match_matrix, red_index):
        return [blue_lines[match_matrix[list(red_index)][i]] for i in range(len(red_index))]

    @staticmethod
    def line_combinations(iterable):
        """
        :param iterable: list.
               list of lists/ndarrays of x positions. E.g. [np.array([...]), np.array([...])...]
        :return: ndarray.
        all possible N-length combinations of the N arrays inside iterable, picking one element
        from each, subject to the constraint that each N-length combination is strictly increasing.
        Each combination is returned as a row vector. E.g. if iterable = [[1,5], [2, 6]],
        the return would be:
        [np.array([[1], [2]]), np.array([[1], [6]]), np.array([[5], [6]])]
        """
        paths = np.array(list(itertools.product(*iterable)))
        paths = paths[np.all(paths[:, :-1] < paths[:, 1:], axis=1)]
        return paths.reshape(-1, len(iterable), 1)

    @staticmethod
    def match_matrix(lines1, lines2, lines1_flux, lines2_flux, tol=0.2):
        # make boolean matrix which indicates where peak fluxes are close.
        flux1_matrix = lines1_flux.reshape(-1, 1) * np.ones((len(lines1), len(lines2)))
        flux2_matrix = np.ones((len(lines1), len(lines2))) * lines2_flux
        # check that the difference of each peak is less than tol times their mean value (representative of the `true` peak value)
        matches = np.abs(flux2_matrix - flux1_matrix) < tol * np.mean([flux2_matrix, flux1_matrix], axis=0)
        # do a cut on arr1 where arr2 fluxes do not match close enough
        # consider cutting elements of arr2 as well which have no flux matches.
        valid_lines1 = lines1[np.sum(matches, axis=1) > 0]
        matches = matches[np.sum(matches, axis=1) > 0]
        return matches, valid_lines1


def blank_overlap_table(max_overlap_red):
    t = Table([Column(name='ref_id', dtype=int),
               Column(name='fiber', dtype=int),
               Column(name='matched_ref_id', dtype=int),
               Column(name='pixel', shape=(max_overlap_red,)),
               Column(name='matched_pixel', shape=(max_overlap_red,)),
               Column(name='peaks'),
               Column(name='good', dtype=bool)])
    return t


def _select_lines(lines, spec, ridx, bidx, max_overlap_red, max_overlap_blue):
    r_select = np.logical_and(lines['order'] == spec['ref_id'][ridx],
                              lines['pixel'] < np.max(spec['pixel'][ridx][:max_overlap_red]))
    b_select = np.logical_and(lines['order'] == spec['ref_id'][bidx],
                              lines['pixel'] > np.min(spec['pixel'][bidx][-max_overlap_blue:]))
    fluxkey = 'corrected_flux' if lines.get('corrected_flux', None) is not None else 'flux'
    return lines['pixel'][r_select], lines['pixel'][b_select], lines[fluxkey][r_select], lines[fluxkey][b_select]


def fit_overlaps(spec, lines, max_overlap_red=1000, max_overlap_blue=2000, linear_scale_range=(0.5, 2),
                 fiber=np.nan, deg=2, flux_tol=0.2):
    """
    :param spec: astropy.table.Table : sf_spec['flux'][i] should be a 1d ndarray with flux for the
                    diffraction order labelled sf_spec['ref_id'][i].
    :param lines: astropy.table.Table : lines['flux'][i], lines['order'][i], lines['pixel'][i]
                  must be the flux, order id, and pixel of an emission line.
    :param max_overlap_red: int. the pixel size of the blue edge of the redder order to consider for overlap
    :param max_overlap_blue: int. the pixel size of the red edge of the bluer order to consider for overlap
    :param linear_scale_range: ordered tuple. The range of linear coefficients to consider for the mapping between
    the red order and the blue order.
    :param flux_tol: float.
                     float between (0, 1]. For two lines to be considered a matched peak, their fluxes must
                     agree to within flux_tol.
    :return: A table of overlaps.

    WARNING: sf_spec[i] must be redder than sf_spec[i+1] if __ is True.
    WARNING: The spectrum of each order must be contiguous and evenly
    sampled over pixel space.
    WARNING: Each order must go from blue to red proceeding from left to right (pixel 0 to pixel 4096)
    """
    max_overlap_red = min(max_overlap_red, spec['pixel'].data.shape[1])
    max_overlap_blue = min(max_overlap_blue, spec['pixel'].data.shape[1])
    overlaps = blank_overlap_table(max_overlap_red)
    fitter = OverlapFitter()
    for row in range(len(spec['flux']) - 1):
        ridx, bidx = row, row+1
        r_lines, b_lines, r_line_flux, b_line_flux = _select_lines(lines, spec, ridx,
                                                                   bidx, max_overlap_red,
                                                                   max_overlap_blue)
        if len(r_lines) < deg + 1 or len(b_lines) < deg + 1:
            # skip to next if there are not enough lines to constrain the fit.
            continue
        overlap = fitter.fit(b_lines, r_lines, b_line_flux, r_line_flux, linear_scale_range,
                             deg=deg, flux_tol=flux_tol)
        pad = max_overlap_red - len(overlap['pixel'])
        overlaps.add_row({'ref_id': spec['ref_id'][ridx], 'matched_ref_id': spec['ref_id'][bidx],
                          'pixel': np.pad(overlap['pixel'].astype(float), ((0, pad)),
                                          mode='constant', constant_values=np.nan),
                          'matched_pixel': np.pad(overlap['matched_pixel'].astype(float), ((0, pad)),
                                                  mode='constant', constant_values=np.nan),
                          'peaks': overlap['peaks'],
                          'good': np.ones_like(overlap['peaks'], dtype=bool),
                          'fiber': fiber * np.ones_like(overlap['peaks'], dtype=int)})
    return overlaps


def match_peaks(b_lines, r_lines, coeffs, pixel_tol=1):
    mapped_red_lines = coordinate_transform(r_lines, coeffs)
    return r_lines[np.isclose(mapped_red_lines, find_nearest(mapped_red_lines, b_lines), atol=pixel_tol)]


def coordinate_transform(coordinates, coefficients):
    return poly.polyval(x=coordinates, c=coefficients)


def flag_bad_overlaps(overlaps, min_num_matches):
    """
    :param overlaps: Table of overlaps with ref_id, matched_ref_id, pixel and
    matched_pixel columns
    :param min_num_matches: int. Minimum number of peak matches for an overlap to count as well fit.
    :return overlaps: Same is input but each bad overlap has had its entry in the 'good' column
                      set to False.
    """
    bad_overlaps = [idx for idx in range(len(overlaps)) if overlaps[idx]['peaks'] < min_num_matches]
    overlaps['good'][bad_overlaps] = False
    return overlaps


def flag_outlier_overlaps(overlaps):
    med = np.median([overlap['ref_id'] for overlap in overlaps if overlap['good']])
    outlier_overlaps = [idx for idx in range(len(overlaps)) if
                        overlaps[idx]['good'] and not np.isclose(overlaps[idx]['ref_id'], med, atol=10)]
    overlaps['good'][outlier_overlaps] = False
    return overlaps

