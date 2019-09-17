import numpy as np
import copy
from scipy import optimize
from astropy.stats import median_absolute_deviation
from astropy.table import Column, Table, vstack
from numpy.polynomial import legendre

from banzai.stages import Stage
from banzai.calibrations import ApplyCalibration
from banzai.images import DataTable
from scipy import interpolate

from banzai_nres.utils.wavelength_utils import identify_lines, calc_residuals, restrict, Model, _sigma_clip
from banzai_nres.utils.wavelength_utils import estimate_global_scale, normalize_coordinates, pixel_order_as_array
from banzai_nres.utils.overlap_utils import flag_bad_overlaps, fit_overlaps, blank_overlap_table, flag_outlier_overlaps
from banzai_nres.utils.misc_utils import brute_local_min, find_nearest, minmax
from banzai_nres.utils.fiber_utils import lit_wavecal_fibers
import banzai_nres.settings as nres

import logging as logger


class WavelengthSolution(object):
    def __init__(self, model=None, model_coefficients=None, measured_lines=None, reference_lines=None,
                 m0=30, max_pixel=100, max_order=2, min_order=0, min_pixel=0, overlap_range=(),
                 grating_eq=True):
        if model is None:
            model = {1: [0, 1],
                     2: [0, 1]}
        self.max_order = max_order
        self.max_pixel = max_pixel
        self.min_pixel = min_pixel
        self.min_order = min_order
        self.model = Model(model)
        self.m0 = m0
        self.overlap_range = overlap_range
        self.model_coefficients = model_coefficients
        self.measured_lines = measured_lines
        self.reference_lines = reference_lines
        self.grating_eq = grating_eq  # True to use a 1/(m0+i) prefactor.

    def wavelength(self, pixel, order):
        normed_pixel = normalize_coordinates(pixel, max_value=self.max_pixel, min_value=self.min_pixel)
        normed_order = normalize_coordinates(order, max_value=self.max_order, min_value=self.min_order)
        return self.wavelength_normed_input(normed_pixel, normed_order, order)

    def wavelength_normed_input(self, normed_pixel, normed_order, order, **kwargs):
        A, c = self._construct_wavelength_map_matrices(normed_pixel.flatten(), normed_order.flatten(), order.flatten())
        wavelengths = self._map_to_wavelength(A, c, self.model_coefficients)
        return wavelengths.reshape(normed_pixel.shape)

    def update_model(self, new_model, grating_eq=True):
        pixel, order = np.meshgrid(np.linspace(self.min_pixel, self.max_pixel, 30),
                                   np.arange(self.min_order, self.max_order + 1))
        coords = {'pixel': pixel.flatten(), 'order': order.flatten()}
        coords['normed_pixel'] = normalize_coordinates(coords['pixel'], self.max_pixel, self.min_pixel)
        coords['normed_order'] = normalize_coordinates(coords['order'], self.max_order, self.min_order)
        old_wavelengths = self.wavelength_normed_input(**coords)
        self.model = Model(new_model)
        self.grating_eq = grating_eq
        self.model_coefficients = self.solve(coords, old_wavelengths)

    def solve(self, line_coordinates, wavelengths_to_fit, weights=np.array([1])):
        """
        :param weights: array.
                        A 1d array of weights. E.g. the square root of the inverse variance.
        :return: array: best fit coefficients
        """
        A, c = self._construct_wavelength_map_matrices(**line_coordinates)
        c += np.array(wavelengths_to_fit).reshape(-1, 1)
        model_coefficients, residuals = np.linalg.lstsq(A * weights.reshape(-1, 1),
                                                        c * weights.reshape(-1, 1), rcond=None)[:2]
        return model_coefficients.flatten()

    def solve_from_overlaps(self, overlaps):
        if 0 in self.model:
            logger.error('Model contains parameters independent of pixel coordinate x. '
                         'Overlap fit will fail. Do not include 0: [...] terms in '
                         'the initial wavelength model.', extra_tags={'model': self.model})
            raise ValueError('0 in self.model.keys(). Not allowed for overlap fit.')
        coordinates = self._format_overlaps(overlaps, pixel_key='pixel', order_key='ref_id')
        matched_coordinates = self._format_overlaps(overlaps, pixel_key='matched_pixel', order_key='matched_ref_id')
        A1, c1 = self._construct_wavelength_map_matrices(**coordinates)
        A2, c2 = self._construct_wavelength_map_matrices(**matched_coordinates)
        A, c = A1 - A2, c1 - c2
        model_coefficients = np.linalg.lstsq(A, c, rcond=None)[0]
        self.model_coefficients = model_coefficients.flatten()

    def apply_scale(self, scale):
        self.model_coefficients *= scale

    def _format_overlaps(self, overlaps, pixel_key='pixel', order_key='ref_id'):
        """
        :param overlaps: An astropy table of overlaps with the 'pixel' and 'matched_pixel' columns
        giving 1d arrays of pixel points which must map to the same wavelength.
        :return: coordinates of overlaps.
        """
        orders = (overlaps[order_key].data.reshape(-1, 1) * np.ones_like(overlaps[pixel_key].data)).flatten()
        pixels = overlaps[pixel_key].data.flatten()
        pixels, orders = pixels[np.isfinite(pixels)], orders[np.isfinite(pixels)]
        coordinates = {'normed_pixel': normalize_coordinates(pixels, self.max_pixel, self.min_pixel), 'pixel': pixels,
                       'normed_order': normalize_coordinates(orders, self.max_order, self.min_order), 'order': orders}
        return coordinates

    def _construct_wavelength_map_matrices(self, normed_pixel, normed_order, order, **kwargs):
        if self.grating_eq:
            m_divisor = 1 / np.add(order, self.m0).astype(np.float64)
        else:
            m_divisor = np.ones_like(normed_pixel, dtype=np.float64)
        columns = ()
        for xpower in self.model:
            for orderpower in self.model[xpower]:
                coefficient_array = np.zeros((xpower+1, orderpower+1))
                coefficient_array[xpower, orderpower] = 1
                columns += (m_divisor * legendre.legval2d(normed_pixel, normed_order, coefficient_array),)
        A = np.column_stack(columns)
        if (0 not in self.model) or (0 not in self.model[0]):
            # we force 0,0 coefficient to 1 for the overlap fit.
            c = -1 * m_divisor.reshape(-1, 1)
        else:
            c = np.zeros_like(normed_pixel, dtype=np.float64).reshape(-1, 1)
        return A, c

    @staticmethod
    def _map_to_wavelength(A, c, coefficients):
        return np.dot(A, coefficients).flatten() + (-1) * c.flatten()


class WavelengthStage(Stage):
    """
    Base class for all wavelength stages
    """
    def __init__(self, runtime_context=None):
        super(WavelengthStage, self).__init__(runtime_context=runtime_context)

    def do_stage(self, image):
        valid_fibers = self._valid_fibers(image)
        for fiber in valid_fibers:
            image = self.do_stage_fiber(image, fiber)
        if len(valid_fibers) < 1:
            self.on_no_valid_fibers(image)
        return image

    def do_stage_fiber(self, image, fiber):
        return image

    @staticmethod
    def on_no_valid_fibers(image):
        logger.warning('No fibers found with non None WavelengthSolution objects. Aborting this '
                       'stage.', image=image)

    @staticmethod
    def _valid_fibers(image):
        return [fiber for fiber in lit_wavecal_fibers(image) if
                image.wavelength_solution[fiber] is not None]


class Initialize(WavelengthStage):
    """
    Stage 1/8 for the wavelength solution
    """
    def __init__(self, runtime_context=None):
        super(Initialize, self).__init__(runtime_context=runtime_context)

    def do_stage_fiber(self, image, fiber):
        logger.info('Appending blank WavelengthSolution object to image for this fiber.',
                    image=image, extra_tags={'fiber': str(fiber)})
        spectrum = image.data_tables[nres.BOX_SPECTRUM_NAME]
        single_fiber_spectrum = spectrum[spectrum['fiber'] == fiber]
        image.wavelength_solution[fiber] = WavelengthSolution(model=nres.initial_wavelength_model,
                                                              m0=nres.principle_order_number,
                                                              max_order=np.max(single_fiber_spectrum['ref_id']),
                                                              min_order=np.min(single_fiber_spectrum['ref_id']),
                                                              max_pixel=np.max(single_fiber_spectrum['pixel']),
                                                              min_pixel=np.min(single_fiber_spectrum['pixel']))
        return image

    @staticmethod
    def on_no_valid_fibers(image):
        for fiber in lit_wavecal_fibers(image):
            image.wavelength_solution[fiber] = None
        logger.error('Image has a length 0 spectrum. Aborting wavelength calibration', image=image)

    @staticmethod
    def _valid_fibers(image):
        return lit_wavecal_fibers(image) if len(image.data_tables[nres.BOX_SPECTRUM_NAME]['flux']) > 0 else []


class AddWavelengthColumn(WavelengthStage):
    """

    """
    def __init__(self, runtime_context=None):
        super(AddWavelengthColumn, self).__init__(runtime_context=runtime_context)
        self.spectrum_table_names = [nres.BOX_SPECTRUM_NAME, nres.BLAZE_CORRECTED_BOX_SPECTRUM_NAME]

    def do_stage(self, image):
        if len(self._valid_fibers(image)) > 0:
            for name in self.spectrum_table_names:
                logger.info('Appending a blank wavelengths column onto {0} data table'.format(name), image=image)
                image.data_tables[name].add_column(Column(np.zeros_like(image.data_tables[name]['flux'], dtype=np.float64),
                                                          unit='angstrom'), name='wavelength')
        return image


class LoadReferenceLineList(ApplyCalibration):
    """
    Stage 2/8 for the wavelength solution
    """
    def __init__(self, runtime_context=None):
        super(LoadReferenceLineList, self).__init__(runtime_context=runtime_context)

    def do_stage(self, image):
        line_list_path = self.get_calibration_filename(image)
        if line_list_path is None:
            self.on_missing_master_calibration(image)
            return image

        line_list = np.sort(np.genfromtxt(line_list_path, usecols=[1]).flatten())
        return self.apply_master_calibration(image, line_list)

    @property
    def calibration_type(self):
        return 'WAVELENGTH'

    def apply_master_calibration(self, image, reference_line_list):
        for fiber in [fiber for fiber in lit_wavecal_fibers(image) if
                      image.wavelength_solution[fiber] is not None]:
            image.wavelength_solution[fiber].reference_lines = reference_line_list
        return image

    def get_calibration_filename(self, image):
        return '/home/gmbrandt/Documents/banzai-nres/banzai_nres/data/ThAr_atlas_ESO.txt'


class FitOverlaps(WavelengthStage):
    """
    Stage 4/8 for the wavelength solution
    This should run on a blaze corrected calibration spectrum.
    """
    def __init__(self, runtime_context=None):
        super(FitOverlaps, self).__init__(runtime_context=runtime_context)

    def do_stage(self, image):
        if image.data_tables.get(nres.OVERLAP_TABLE_NAME, None) is None:
            image.data_tables[nres.OVERLAP_TABLE_NAME] = blank_overlap_table(nres.max_red_overlap)

        image = super(FitOverlaps, self).do_stage(image)

        name = nres.OVERLAP_TABLE_NAME
        image.data_tables[name] = DataTable(image.data_tables[name], name=name)
        return image

    def do_stage_fiber(self, image, fiber):
        logger.info('Fitting overlaps.', image=image, extra_tags={'fiber': str(fiber)})
        spectrum = image.data_tables[nres.BOX_SPECTRUM_NAME]
        single_fiber_spectrum = spectrum[spectrum['fiber'] == fiber]
        overlaps = fit_overlaps(spec=single_fiber_spectrum,
                                lines=image.wavelength_solution[fiber].measured_lines,
                                max_overlap_red=nres.max_red_overlap,
                                max_overlap_blue=nres.max_blue_overlap,
                                linear_scale_range=nres.overlap_linear_scale_range,
                                fiber=fiber)
        overlaps = flag_bad_overlaps(overlaps)
        logger.info('{0} overlaps verified.'.format(np.count_nonzero(overlaps['good'])),
                    image=image, extra_tags={'fiber': str(fiber)})
        overlaps = flag_outlier_overlaps(overlaps)
        logger.info('{0} overlaps will be used.'.format(np.count_nonzero(overlaps['good'])),
                    image=image, extra_tags={'fiber': str(fiber)})

        image.data_tables[nres.OVERLAP_TABLE_NAME] = vstack([overlaps,
                                                             image.data_tables[nres.OVERLAP_TABLE_NAME]])
        if np.count_nonzero(overlaps['good']) < nres.min_num_overlaps:
            logger.error('Less than {0} overlaps verified as good,'
                         'setting wavelength solution to None.'.format(nres.min_num_overlaps),
                         image=image, extra_tags={'fiber': str(fiber)})
            image.wavelength_solution[fiber] = None
        return image


class SolveFromOverlaps(WavelengthStage):
    def __init__(self, runtime_context=None):
        super(SolveFromOverlaps, self).__init__(runtime_context=runtime_context)

    def do_stage_fiber(self, image, fiber):
        image.wavelength_solution[fiber].model = Model(nres.initial_wavelength_model)
        overlaps = image.data_tables.get(nres.OVERLAP_TABLE_NAME, blank_overlap_table(1))
        overlaps = self._prune_overlaps(overlaps, fiber)
        logger.info('Initializing wavelength solution from overlaps.',
                    image=image, extra_tags={'fiber': str(fiber)})
        image.wavelength_solution[fiber].overlap_range = minmax([overlaps['ref_id'], overlaps['matched_ref_id']])
        image.wavelength_solution[fiber].solve_from_overlaps(overlaps)
        return image

    @staticmethod
    def _prune_overlaps(overlaps, fiber):
        overlaps = overlaps[overlaps['good']]
        overlaps = overlaps[overlaps['fiber'] == fiber]
        return overlaps


class IdentifyArcEmissionLines(WavelengthStage):
    """
    Stage 3/8 for the wavelength solution
    """
    def __init__(self, runtime_context=None):
        super(IdentifyArcEmissionLines, self).__init__(runtime_context=runtime_context)
        self.min_peak_snr = nres.min_peak_snr

    def do_stage_fiber(self, image, fiber):
        spectrum = image.data_tables[nres.BOX_SPECTRUM_NAME]
        single_fiber_spectrum = spectrum[spectrum['fiber'] == fiber]
        measured_lines = identify_lines(spectrum=single_fiber_spectrum,
                                        stderr=100,
                                        min_snr=self.min_peak_snr,
                                        order_key='ref_id')

        measured_lines['normed_order'] = normalize_coordinates(measured_lines['order'],
                                                               image.wavelength_solution[fiber].max_order,
                                                               image.wavelength_solution[fiber].min_order)
        measured_lines['normed_pixel'] = normalize_coordinates(measured_lines['pixel'],
                                                               image.wavelength_solution[fiber].max_pixel,
                                                               image.wavelength_solution[fiber].min_pixel)
        image.wavelength_solution[fiber].measured_lines = measured_lines
        logger.info('{0} emission lines identified from {1} unique '
                    'diffraction orders'.format(len(measured_lines['pixel']), len(set(measured_lines['order']))),
                    image=image, extra_tags={'fiber': str(fiber)})
        return image
    # TODO modify valid_fibers to include a check if nres.BOX_SPECTRUM_NAME is not None.


class BlazeCorrectArcEmissionLines(WavelengthStage):
    def __init__(self, runtime_context=None):
        super(BlazeCorrectArcEmissionLines, self).__init__(runtime_context=runtime_context)

    def do_stage_fiber(self, image, fiber):
        spectrum = image.data_tables.get(nres.BLAZE_CORRECTED_BOX_SPECTRUM_NAME)
        lines = image.wavelength_solution[fiber].measured_lines
        if spectrum is None:
            image.wavelength_solution[fiber].measured_lines['corrected_flux'] = lines['flux']
        else:
            spectrum = spectrum[spectrum['fiber'] == fiber]
            lines['corrected_flux'] = np.zeros_like(lines['flux'], dtype=float)
            for spec in spectrum:
                flux = interpolate.interp1d(spec['pixel'], spec['flux'], kind='nearest')
                in_order = np.where(lines['order'] == spec['ref_id'])
                lines['corrected_flux'][in_order] = flux(lines['pixel'][in_order])
            image.wavelength_solution[fiber].measured_lines = lines
        return image


class IdentifyArcEmissionLinesLowSN(IdentifyArcEmissionLines):
    def __init__(self, runtime_context=None):
        super(IdentifyArcEmissionLinesLowSN, self).__init__(runtime_context=runtime_context)
        self.min_peak_snr = 3


class FindGlobalScale(WavelengthStage):
    """
    Stage 5/8 for the wavelength solution
    """
    def __init__(self, runtime_context=None):
        super(FindGlobalScale, self).__init__(runtime_context=runtime_context)

    def do_stage_fiber(self, image, fiber):
        scale_guess = estimate_global_scale(detector_range=nres.approx_detector_range_angstroms,
                                            n=nres.approx_num_orders,
                                            m0=image.wavelength_solution[fiber].m0)
        scale = self._find_scale(image.wavelength_solution[fiber], scale_guess, nres.global_scale_range)
        image.wavelength_solution[fiber].update_model(nres.intermediate_wavelength_model)
        image.wavelength_solution[fiber].apply_scale(scale)
        logger.info('The scale guess was {:.6e} and the search yielded {:.6e}'.format(scale_guess, scale),
                    image=image, extra_tags={'fiber': str(fiber)})
        if not np.isclose(scale_guess, scale, rtol=2):
            logger.error('Global scale is more than a factor of two away from initial guess, '
                         'an error in the wavelength solution for this fiber is likely.',
                         image=image, extra_tags={'fiber': str(fiber)})
        return image

    @staticmethod
    def _find_scale(wcs, scale_guess, rrange=(0.1, 5), step=10):
        unscaled_wavelengths = wcs.wavelength_normed_input(**restrict(wcs.measured_lines, *wcs.overlap_range))
        args = (unscaled_wavelengths, np.sort(wcs.reference_lines))
        scale = optimize.minimize(fun=FindGlobalScale._chi_squared_safe, x0=scale_guess, args=args,
                                  method=brute_local_min, options={'step': step, 'rrange': rrange, 'filtw': 501,
                                                                   'finish': 'Nelder-Mead'}).x[0]
        return scale

    @staticmethod
    def _chi_squared_safe(scales, unscaled_wavelengths, reference_wavelengths, mem_limit=1E8):
        chi2 = np.zeros_like(scales, dtype=np.float32)
        step = max(int(len(scales) * len(unscaled_wavelengths)) // int(mem_limit * 8 / 32), 1)  # bits / bit_limit
        for i in range(step):
            chi2[i::step] = FindGlobalScale._chi_squared(scales[i::step], unscaled_wavelengths, reference_wavelengths)
        return chi2

    @staticmethod
    def _chi_squared(scales, unscaled_wavelengths, reference_wavelengths):
        wavelengths = np.array(scales).reshape(-1, 1) * unscaled_wavelengths
        residuals = calc_residuals(wavelengths.flatten(), reference_wavelengths).reshape(wavelengths.shape)
        return np.sum(residuals ** 2, axis=1).flatten()
        #return np.sum(np.sort(residuals**2, axis=1)[:, :nmax], axis=1).flatten()

class SolutionRefineInitial(WavelengthStage):
    """
    Stage 6/8 for the wavelength solution
    """
    def __init__(self, runtime_context=None):
        super(SolutionRefineInitial, self).__init__(runtime_context=runtime_context)

    def do_stage_fiber(self, image, fiber):
        image.wavelength_solution[fiber].update_model(new_model=nres.intermediate_wavelength_model)
        image.wavelength_solution[fiber], rsd = self.constrain_solution_over_detector(image.wavelength_solution[fiber])

        mad, std = median_absolute_deviation(rsd), np.std(rsd)
        logger.info('median absolute deviation is {0} and the standard deviation is {1}'.format(mad, std),
                    image=image, extra_tags={'fiber': str(fiber)})
        logger.info('{0} lines within 4.5 median absolute deviations and {1} lines within 4.5 standard deviations'
                    ''.format(np.count_nonzero(np.isclose(rsd, 0, atol=4.5*mad)),
                              np.count_nonzero(np.isclose(rsd, 0, atol=4.5*std))))
        return image

    @staticmethod
    def constrain_solution_over_detector(wcs):
        """
        :param wcs: WavelengthSolution
        :return:     (WavelengthSolution, array)
        wavelength solution, residuals. Each entry of residuals is the difference in
        wavelength between a measured line and its corresponding closest reference line.
        """
        reference_lines = np.sort(wcs.reference_lines)
        residuals = []
        for order_range in SolutionRefineInitial._ranges_to_evaluate(*wcs.overlap_range, wcs.min_order, wcs.max_order):
            wcs, residuals = refine_wcs(wcs, wcs.measured_lines, reference_lines,
                                        SolutionRefineInitial._converged_when_max_iter,
                                        SolutionRefineInitial._clip,
                                        kwargs={'sigma': 6,
                                                'stdfunc': median_absolute_deviation,
                                                'order_range': order_range},
                                        max_iter=5)
        return wcs, residuals

    @staticmethod
    def _converged_when_max_iter(*args, **kwargs):
        return False

    @staticmethod
    def _clip(wcs, lines, reference_lines, sigma=3, stdfunc=np.std, order_range=(), **kwargs):
        lines = restrict(lines, *order_range)
        residuals = calc_residuals(wcs.wavelength_normed_input(**lines), reference_lines)
        lines = _sigma_clip(residuals, lines, sigma=sigma, stdfunc=stdfunc)
        return lines

    @staticmethod
    def _ranges_to_evaluate(low, high, min_order, max_order):
        return [(low - i, high + i) for i in range(max(low - min_order, max_order - high) + 1)]


def refine_wcs(wcs, measured_lines, reference_lines, converged, clip_fun, kwargs=None, max_iter=20):
    """
    :param wcs: WavelengthSolution
    :param reference_lines: array. Wavelengths of reference calibration lines.
    :param measured_lines: dict: {'normed_pixel': array, 'normed_order': array}. Coordinates of measured lines.
    :param kwargs: keyword arguments to feed to clip_fun and converged.
    :param max_iter: the max number of iterations allowed to reach convergence.
    :param converged: function which accepts wcs, residuals, last_residuals, measured_lines
                      and **kwargs and returns True or False. If True, the solve loop will stop.
    :param clip_fun: function which accepts residuals, measured_lines, reference_lines and **kwargs and returns
    a dict in the same format as measured_lines.
    :return:     (WavelengthSolution, array)
    wavelength solution, residuals. Each entry of residuals is the difference in
    wavelength between a measured line and its corresponding closest reference line.
    """
    if kwargs is None:
        kwargs = {}
    residuals, last_residuals = [], np.inf
    # note we could evaluate the find_nearest function out here instead of in the loop.
    for iteration in range(max_iter):
        clipped_lines = clip_fun(wcs, measured_lines, reference_lines, **kwargs)
        closest_ref_lines = find_nearest(wcs.wavelength_normed_input(**clipped_lines), array_b=reference_lines)
        wcs.model_coefficients = wcs.solve(clipped_lines, closest_ref_lines)
        residuals = calc_residuals(wcs.wavelength_normed_input(**measured_lines), reference_lines)
        if converged(residuals, last_residuals, **kwargs):
            break
        last_residuals = np.copy(residuals)
    return wcs, residuals


class SolutionRefineFinal(WavelengthStage):
    """
    Stage 7/8 for the wavelength solution
    """
    def __init__(self, runtime_context=None):
        super(SolutionRefineFinal, self).__init__(runtime_context=runtime_context)

    def do_stage_fiber(self, image, fiber):
        image.wavelength_solution[fiber], rsd = self._refine(image.wavelength_solution[fiber],
                                                             nres.final_wavelength_model)

        mad, std = median_absolute_deviation(rsd), np.std(rsd)
        logger.info('median absolute deviation is {0} and the standard deviation is {1}'.format(mad, std),
                    image=image, extra_tags={'fiber': str(fiber)})
        logger.info('{0} lines within 4.5 median absolute deviations and {1} lines within 4.5 standard deviations'
                    ''.format(np.count_nonzero(np.isclose(rsd, 0, atol=4.5*mad)),
                              np.count_nonzero(np.isclose(rsd, 0, atol=4.5*std))))
        return image

    @staticmethod
    def _refine(wavelength_solution, final_model):
        """
        :param wavelength_solution: WavelengthSolution
        :param final_model: dict. The model which describes the form of the mapping from pixel and order
                            to wavelength. See __init__ of WavelengthSolution.
        :return:     (WavelengthSolution, array)
        wavelength solution, residuals. Each entry of residuals is the difference in
        wavelength between a measured line and its corresponding closest reference line.
        """
        reference_lines = np.sort(wavelength_solution.reference_lines)
        measured_lines = wavelength_solution.measured_lines
        residuals = []
        for x_degree in final_model:
            for i_degree in final_model[x_degree]:
                if wavelength_solution.model.is_missing_polynomial_term(x_degree, i_degree):
                    # This will only refine if the model is missing the term,
                    #  which means this stage does nothing if you want to rerun a model.
                    new_model = copy.deepcopy(wavelength_solution.model)
                    new_model.add_polynomial_term(x_degree, i_degree)
                    wavelength_solution.update_model(new_model=new_model)
                    wavelength_solution, residuals = refine_wcs(wavelength_solution,
                                                                measured_lines, reference_lines,
                                                                SolutionRefineFinal._converged,
                                                                SolutionRefineFinal._clip,
                                                                max_iter=20,
                                                                kwargs={'sigma': 4,
                                                                        'stdfunc': median_absolute_deviation})
        return wavelength_solution, residuals

    @staticmethod
    def _converged(residuals, last_residuals, **kwargs):
        if np.allclose(residuals, last_residuals, atol=1E-3):
            return True
        return False

    @staticmethod
    def _clip(wcs, lines, reference_lines, sigma=3, stdfunc=np.std, **kwargs):
        residuals = calc_residuals(wcs.wavelength_normed_input(**lines), reference_lines)
        lines = _sigma_clip(residuals, lines, sigma=sigma, stdfunc=stdfunc)
        return lines


class ApplyToSpectrum(WavelengthStage):
    """
    Stage 8/8 for the wavelength solution
    """
    def __init__(self, runtime_context=None):
        super(ApplyToSpectrum, self).__init__(runtime_context=runtime_context)

    def do_stage_fiber(self, image, fiber):
        spectrum = image.data_tables[nres.BOX_SPECTRUM_NAME]
        fiber_mask = np.where(spectrum['fiber'] == fiber)
        wcs = image.wavelength_solution[fiber]
        pixel_coordinates, order_coordinates = pixel_order_as_array(spectrum[fiber_mask])
        spectrum['wavelength'][fiber_mask] = wcs.wavelength(pixel=pixel_coordinates,
                                                            order=order_coordinates)
        image.data_tables[nres.BOX_SPECTRUM_NAME] = spectrum
        return image


class TabulateArcEmissionLines(WavelengthStage):
    """
    Optional stage which will create a table of all identified emission lines and their matches
    in the reference list. This is saved when the image reduction is completed, under the extension LINES.
    """
    def __init__(self, runtime_context=None):
        super(TabulateArcEmissionLines, self).__init__(runtime_context=runtime_context)

    def do_stage(self, image):
        valid_fibers = self._valid_fibers(image)
        if len(valid_fibers) > 0:
            lines = self._format_lines(image, valid_fibers)
            for fiber in valid_fibers:
                wcs = image.wavelength_solution[fiber]
                fib = np.where(lines['fiber'] == fiber)
                lines['wavelength'][fib] = wcs.wavelength(lines['pixel'][fib].data, lines['order'][fib].data)
                lines['reference_wavelength'][fib] = find_nearest(lines['wavelength'][fib], wcs.reference_lines)

            image.data_tables[nres.EMISSION_LINES_TABLE_NAME] = DataTable(lines, nres.EMISSION_LINES_TABLE_NAME)
        return image

    @staticmethod
    def _format_lines(image, fibers):
        all_lines = []
        for fiber in fibers:
            lines = Table(image.wavelength_solution[fiber].measured_lines)
            lines.add_column(Column(fiber * np.ones_like(lines['pixel']), dtype=int), name='fiber')
            lines.add_column(Column(np.zeros_like(lines['pixel'], dtype=np.float64), unit='angstrom'), name='wavelength')
            lines.add_column(Column(np.zeros_like(lines['pixel'], dtype=np.float64), unit='angstrom'), name='reference_wavelength')
            lines['reference_wavelength'].description = 'wavelength of the nearest match in the reference line list'
            lines['wavelength'].description = 'wavelength of measured line from the calibration frame'
            all_lines.append(lines)
        return vstack(all_lines)


class IdentifyPrincipleOrderNumber(WavelengthStage):
    """
    Stage for identifying the principle order number m0. This should be run between IdentifyArcEmissionLines and
    SolveFromOverlaps. This stage is meant to be run once on a batch of stages to find the principle order
    number. The principle order number should then be fixed in settings.py.
    """
    def __init__(self, runtime_context=None):
        super(IdentifyPrincipleOrderNumber, self).__init__(runtime_context=runtime_context)
        self.STAGES_TODO = [SolveFromOverlaps, FindGlobalScale, SolutionRefineInitial, SolutionRefineFinal]

    def do_stage_fiber(self, image, fiber):
        logger.info('Looking for the principle order number between {0} and {1}'.format(*nres.m0_range),
                     image=image, extra_tags={'fiber': str(fiber)})
        logger.disabled = True
        merits, m0_values = self.merit_per_m0(image, fiber, nres.m0_range)
        logger.disabled = False
        best_m0, merit = m0_values[np.argmin(merits)], np.min(merits)

        if not merit < 1/10 * np.median(merits):
            logger.warning('A definitive principle order number was not found. Aborting wavelength solution',
                            image=image, extra_tags={'fiber': str(fiber)})
            image.wavelength_solution[fiber] = None
        else:
            image.wavelength_solution[fiber].m0 = best_m0
            logger.info('The best principle order number is {0}'.format(best_m0),
                        image=image, extra_tags={'fiber': str(fiber)})
        return image

    def merit_per_m0(self, image, fiber, m0_range):
        m0_values = np.arange(*m0_range)
        merit = []
        for m0 in m0_values:
            image.wavelength_solution[fiber].m0 = m0
            for STAGE in self.STAGES_TODO:
                image = STAGE(self.runtime_context).do_stage_fiber(image, fiber)

            wcs = image.wavelength_solution[fiber]
            if wcs is not None:
                merit.append(median_absolute_deviation(calc_residuals(wcs.wavelength_normed_input(**wcs.measured_lines),
                                                             wcs.reference_lines)))
        return np.array(merit), m0_values
