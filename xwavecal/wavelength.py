import numpy as np
import copy
from scipy import optimize
from astropy.stats import median_absolute_deviation
from astropy.table import Column, Table, vstack
from numpy.polynomial import legendre

from scipy import interpolate

from xwavecal.images import Image
from xwavecal.stages import Stage, ApplyCalibration
from xwavecal.utils.wavelength_utils import identify_lines, calc_residuals, restrict, Model, _sigma_clip
from xwavecal.utils.wavelength_utils import estimate_global_scale, normalize_coordinates, pixel_order_as_array
from xwavecal.utils.overlap_utils import flag_bad_overlaps, fit_overlaps, blank_overlap_table, flag_outlier_overlaps
from xwavecal.utils.misc_utils import brute_local_min, find_nearest, minmax
from xwavecal.utils.fiber_utils import lit_wavecal_fibers

import logging
logger = logging.getLogger(__name__)


class WavelengthSolution(object):
    """
    Class to evaluate, contain and update the wavelength solution.

    Notable attributes:
        WavelengthSolution.model gives the current model of the wavelength solution.
        WavelengthSolution.model_coefficients gives the coefficients for the current model, in the order
        of WavelengthSolution.model.
    """
    def __init__(self, model=None, model_coefficients=None, measured_lines: dict = None, reference_lines=None,
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
        self.reference_lines = reference_lines
        self.grating_eq = grating_eq  # True to use a 1/(m0+i) prefactor.
        self._measured_lines = None
        self.measured_lines = measured_lines

    @property
    def measured_lines(self):
        return self._measured_lines

    @measured_lines.setter
    def measured_lines(self, lines):
        if lines is not None and lines.get('order', None) is not None and lines.get('pixel', None) is not None:
            lines['normed_pixel'] = normalize_coordinates(lines['pixel'], max_value=self.max_pixel, min_value=self.min_pixel)
            lines['normed_order'] = normalize_coordinates(lines['order'], max_value=self.max_order, min_value=self.min_order)
        self._measured_lines = lines

    def __call__(self, pixel, order):
        """
        :param pixel: ndarray. array (any shape) of pixel coordinates of spectral features.
        :param order: ndarray. Same length as pixel. array of order indices of spectral features.
        :return: ndarray. array of wavelengths. Same shape as pixel and order.
        """
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

    def solve(self, line_coordinates, wavelengths_to_fit, weights=None):
        """
        :param weights: array.
                        A 1d array of weights. E.g. the square root of the inverse variance.
        :return: array: best fit coefficients
        """
        weights = np.array([1]) if weights is None else weights
        A, c = self._construct_wavelength_map_matrices(**line_coordinates)
        c += np.array(wavelengths_to_fit).reshape(-1, 1)
        model_coefficients, residuals = np.linalg.lstsq(A * weights.reshape(-1, 1),
                                                        c * weights.reshape(-1, 1), rcond=None)[:2]
        return model_coefficients.flatten()

    def solve_from_overlaps(self, overlaps):
        if 0 in self.model:
            logger.error('Model contains parameters independent of pixel coordinate x. '
                         'Overlap fit will fail. Do not include 0: [...] terms in '
                         'the initial wavelength model. The model is {0}'.format(self.model))
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
        if coefficients is None:
            return np.ones_like(c.flatten(), dtype=float) * np.nan
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

    def on_no_valid_fibers(self, image):
        self.logger.warning('No fibers found with non None WavelengthSolution objects. Aborting this '
                            'stage.')

    @staticmethod
    def _valid_fibers(image):
        return [fiber for fiber in lit_wavecal_fibers(image) if
                image.wavelength_solution[fiber] is not None]


class Initialize(WavelengthStage):
    """
    Stage 1/9 for the wavelength solution
    """
    def __init__(self, runtime_context=None):
        super(Initialize, self).__init__(runtime_context=runtime_context)

    def do_stage_fiber(self, image, fiber):
        self.logger.info('Appending blank WavelengthSolution object to image for this fiber. '
                         'fiber={0}'.format(str(fiber)))
        spectrum = image.data_tables[self.runtime_context.main_spectrum_name]
        single_fiber_spectrum = spectrum[spectrum['fiber'] == fiber]
        image.wavelength_solution[fiber] = WavelengthSolution(model=self.runtime_context.initial_wavelength_model,
                                                              m0=self.runtime_context.principle_order_number,
                                                              max_order=np.max(single_fiber_spectrum['ref_id']),
                                                              min_order=np.min(single_fiber_spectrum['ref_id']),
                                                              max_pixel=np.max(single_fiber_spectrum['pixel']),
                                                              min_pixel=np.min(single_fiber_spectrum['pixel']))
        return image

    def on_no_valid_fibers(self, image):
        for fiber in lit_wavecal_fibers(image):
            image.wavelength_solution[fiber] = None
        self.logger.error('Aborting wavelength calibration and setting wavelength solution to None.')

    def _valid_fibers(self, image):
        if image.data_tables.get(self.runtime_context.main_spectrum_name) is None:
            self.logger.error('No main spectrum on image.')
            return []

        colnames_present = all([key in image.data_tables[self.runtime_context.main_spectrum_name].colnames
                                for key in ['ref_id', 'fiber']])
        if not colnames_present:
            self.logger.error('Main spectrum on image is missing ref_id and fiber columns')
            return []
        return lit_wavecal_fibers(image)


class LoadReferenceLineList(ApplyCalibration):
    """
    Stage 2/8 for the wavelength solution
    """
    def __init__(self, runtime_context=None):
        super(LoadReferenceLineList, self).__init__(runtime_context=runtime_context)

    def apply_master_calibration(self, image, reference_list_path):
        line_list = np.sort(np.genfromtxt(reference_list_path, usecols=[1]).flatten())
        for fiber in [fiber for fiber in lit_wavecal_fibers(image) if
                      image.wavelength_solution[fiber] is not None]:
            image.wavelength_solution[fiber].reference_lines = line_list
        image.set_header_val('IDLIST', (reference_list_path, 'ID of the reference line list'))
        return image

    def get_calibration_filename(self, image):
        return self.runtime_context.line_list_path


class FitOverlaps(WavelengthStage):
    """
    Stage 3/8 for the wavelength solution
    This should run on a blaze corrected calibration spectrum.
    """
    def __init__(self, runtime_context=None):
        super(FitOverlaps, self).__init__(runtime_context=runtime_context)

    def do_stage(self, image):
        if image.data_tables.get(self.runtime_context.overlap_table_name, None) is None:
            image.data_tables[self.runtime_context.overlap_table_name] = blank_overlap_table(self.runtime_context.max_red_overlap)

        image = super(FitOverlaps, self).do_stage(image)

        name = self.runtime_context.overlap_table_name
        image.data_tables[name] = Table(image.data_tables[name])
        return image

    def do_stage_fiber(self, image, fiber):
        self.logger.info('Fitting overlaps. fiber={0}'.format(str(fiber)))
        spectrum = image.data_tables[self.runtime_context.main_spectrum_name]
        single_fiber_spectrum = spectrum[spectrum['fiber'] == fiber]
        overlaps = fit_overlaps(spec=single_fiber_spectrum,
                                lines=image.wavelength_solution[fiber].measured_lines,
                                max_overlap_red=self.runtime_context.max_red_overlap,
                                max_overlap_blue=self.runtime_context.max_blue_overlap,
                                linear_scale_range=self.runtime_context.overlap_linear_scale_range,
                                fiber=fiber,
                                flux_tol=getattr(self.runtime_context, 'flux_tol', 0.2))
        overlaps = flag_bad_overlaps(overlaps, getattr(self.runtime_context, 'min_num_matches', 6))
        self.logger.info('{0} overlaps verified. fiber={1}'.format(np.count_nonzero(overlaps['good']), str(fiber)))
        overlaps = flag_outlier_overlaps(overlaps)
        self.logger.info('{0} overlaps will be used. fiber={1}'.format(np.count_nonzero(overlaps['good']), str(fiber)))

        image.data_tables[self.runtime_context.overlap_table_name] = vstack([overlaps,
                                                             image.data_tables[self.runtime_context.overlap_table_name]])
        if np.count_nonzero(overlaps['good']) < self.runtime_context.min_num_overlaps:
            self.logger.error('Less than {0} overlaps verified as good,'
                         'setting wavelength solution to None.'
                         ' fiber={1}'.format(self.runtime_context.min_num_overlaps, str(fiber)))
            image.wavelength_solution[fiber] = None
        return image


class SolveFromOverlaps(WavelengthStage):
    """
    Stage 4/8. Solves for the coefficients of the wavelength solution from the overlaps.
    """
    def __init__(self, runtime_context=None):
        super(SolveFromOverlaps, self).__init__(runtime_context=runtime_context)

    def do_stage_fiber(self, image, fiber):
        image.wavelength_solution[fiber].model = Model(self.runtime_context.initial_wavelength_model)
        overlaps = image.data_tables.get(self.runtime_context.overlap_table_name, blank_overlap_table(1))
        overlaps = self._prune_overlaps(overlaps, fiber)
        self.logger.info('Initializing wavelength solution from overlaps. fiber={0}'.format(str(fiber)))
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
    Stage 5/8 for the wavelength solution
    """
    def __init__(self, runtime_context=None):
        super(IdentifyArcEmissionLines, self).__init__(runtime_context=runtime_context)
        self.min_peak_snr = self.runtime_context.min_peak_snr

    def do_stage_fiber(self, image, fiber):
        spectrum = image.data_tables[self.runtime_context.main_spectrum_name]
        single_fiber_spectrum = spectrum[spectrum['fiber'] == fiber]
        measured_lines = identify_lines(spectrum=single_fiber_spectrum,
                                        stderr=single_fiber_spectrum['stderr'],
                                        min_snr=self.min_peak_snr,
                                        order_key='ref_id')
        image.wavelength_solution[fiber].measured_lines = measured_lines
        self.logger.info('{0} emission lines identified from {1} unique '
                    'diffraction orders. fiber={2}'
                    ''.format(len(measured_lines['pixel']), len(set(measured_lines['order'])), str(fiber)))
        return image


class BlazeCorrectArcEmissionLines(WavelengthStage):
    def __init__(self, runtime_context=None):
        super(BlazeCorrectArcEmissionLines, self).__init__(runtime_context=runtime_context)

    def do_stage_fiber(self, image, fiber):
        spectrum = image.data_tables.get(self.runtime_context.blaze_corrected_spectrum_name)
        lines = image.wavelength_solution[fiber].measured_lines
        if spectrum is None:
            image.wavelength_solution[fiber].measured_lines['corrected_flux'] = lines['flux']
        else:
            spectrum = spectrum[spectrum['fiber'] == fiber]
            lines['corrected_flux'] = np.zeros_like(lines['flux'], dtype=float)
            for spec in spectrum:
                flux = interpolate.interp1d(spec['pixel'], spec['flux'], kind='nearest', bounds_error=False, fill_value=0)
                in_order = np.where(lines['order'] == spec['ref_id'])
                lines['corrected_flux'][in_order] = flux(lines['pixel'][in_order])
            image.wavelength_solution[fiber].measured_lines = lines
        return image


class IdentifyArcEmissionLinesLowSN(IdentifyArcEmissionLines):
    def __init__(self, runtime_context=None):
        super(IdentifyArcEmissionLinesLowSN, self).__init__(runtime_context=runtime_context)
        self.min_peak_snr = self.runtime_context.overlap_min_peak_snr


class FindGlobalScale(WavelengthStage):
    """
    Stage 6/8 for the wavelength solution
    """
    def __init__(self, runtime_context=None):
        super(FindGlobalScale, self).__init__(runtime_context=runtime_context)
        self.step = getattr(self.runtime_context, 'global_scale_spacing', 10)

    def do_stage_fiber(self, image, fiber):
        scale_guess = estimate_global_scale(detector_range=self.runtime_context.approx_detector_range_angstroms,
                                            n=self.runtime_context.approx_num_orders,
                                            m0=image.wavelength_solution[fiber].m0)
        scale = self._find_scale(image.wavelength_solution[fiber], scale_guess, self.runtime_context.global_scale_range,
                                 step=self.step)
        image.wavelength_solution[fiber].update_model(self.runtime_context.intermediate_wavelength_model)
        image.wavelength_solution[fiber].apply_scale(scale)
        self.logger.info('The scale guess was {0:.6e} and the search yielded {1:.6e}. fiber={2}'
                    ''.format(scale_guess, scale, str(fiber)))
        if not np.isclose(scale, scale_guess, rtol=2):
            self.logger.error('Global scale is more than a factor of two away from initial guess, '
                         'an error in the wavelength solution for this fiber is likely. fiber={0}'.format(str(fiber)))
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
    Stage 7/8 for the wavelength solution
    """
    def __init__(self, runtime_context=None):
        super(SolutionRefineInitial, self).__init__(runtime_context=runtime_context)
        self.kappa = getattr(runtime_context, 'initial_mad_clip', 6)

    def do_stage_fiber(self, image, fiber):
        image.wavelength_solution[fiber].update_model(new_model=self.runtime_context.intermediate_wavelength_model)
        image.wavelength_solution[fiber], rsd = self.constrain_solution_over_detector(image.wavelength_solution[fiber],
                                                                                      self.kappa)

        mad, std = median_absolute_deviation(rsd), np.std(rsd)
        self.logger.info('median absolute deviation is {0} and the standard deviation is {1}.'
                    ' fiber={2}'.format(mad, std, str(fiber)))
        self.logger.info('{0} lines within 4.5 median absolute deviations and {1} lines within 4.5 standard deviations'
                    ''.format(np.count_nonzero(np.isclose(rsd, 0, atol=4.5*mad)),
                              np.count_nonzero(np.isclose(rsd, 0, atol=4.5*std))))
        return image

    @staticmethod
    def constrain_solution_over_detector(wcs, kappa=6):
        """
        :param wcs: WavelengthSolution
        :param kappa: outliers with residuals exceeding kappa m.a.d. from 0 will not be used to constrain the solution.
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
                                        kwargs={'sigma': kappa,
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
        wcs.model_coefficients = wcs.solve(clipped_lines, closest_ref_lines, clipped_lines.get('weight', None))
        residuals = calc_residuals(wcs.wavelength_normed_input(**measured_lines), reference_lines)
        if converged(residuals, last_residuals, **kwargs):
            break
        last_residuals = np.copy(residuals)
    return wcs, residuals


class SolutionRefineFinal(WavelengthStage):
    """
    Stage 8/8 for the wavelength solution
    """
    def __init__(self, runtime_context=None):
        super(SolutionRefineFinal, self).__init__(runtime_context=runtime_context)
        self.kappa = getattr(runtime_context, 'final_mad_clip', 4)

    def do_stage_fiber(self, image, fiber):
        image.wavelength_solution[fiber], rsd = self._refine(image.wavelength_solution[fiber],
                                                             self.runtime_context.final_wavelength_model, self.kappa)

        mad, std = median_absolute_deviation(rsd), np.std(rsd)
        self.logger.info('median absolute deviation is {0} and the standard deviation is {1}.'
                    ' fiber={2}'.format(mad, std, str(fiber)))
        self.logger.info('{0} lines within 4.5 median absolute deviations and {1} lines within 4.5 standard deviations'
                    ''.format(np.count_nonzero(np.isclose(rsd, 0, atol=4.5*mad)),
                              np.count_nonzero(np.isclose(rsd, 0, atol=4.5*std))))
        return image

    @staticmethod
    def _refine(wavelength_solution, final_model, kappa=4):
        """
        :param wavelength_solution: WavelengthSolution
        :param final_model: dict. The model which describes the form of the mapping from pixel and order
                            to wavelength. See __init__ of WavelengthSolution.
        :param kappa: outliers with residuals exceeding kappa m.a.d. from 0 will not be used to constrain the solution.
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
                                                                kwargs={'sigma': kappa,
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


class SolutionRefineOnce(SolutionRefineFinal):
    """
    Single iteration of refining the wavelength solution. Useful if more lines have been added
    to the line list after final refine.
    """
    def __init__(self, runtime_context=None):
        super(SolutionRefineOnce, self).__init__(runtime_context=runtime_context)

    def do_stage_fiber(self, image, fiber):
        image.wavelength_solution[fiber], rsd = refine_wcs(image.wavelength_solution[fiber],
                                                           image.wavelength_solution[fiber].measured_lines,
                                                           image.wavelength_solution[fiber].reference_lines,
                                                           self._converged, self._clip, max_iter=20,
                                                           kwargs={'sigma': 4,
                                                                   'stdfunc': median_absolute_deviation})

        mad, std = median_absolute_deviation(rsd), np.std(rsd)
        self.logger.info('median absolute deviation is {0} and the standard deviation is {1}.'
                    ' fiber={2}'.format(mad, std, str(fiber)))
        self.logger.info('{0} lines within 4.5 median absolute deviations and {1} lines within 4.5 standard deviations'
                    ''.format(np.count_nonzero(np.isclose(rsd, 0, atol=4.5*mad)),
                              np.count_nonzero(np.isclose(rsd, 0, atol=4.5*std))))
        return image


class ApplyToSpectrum(WavelengthStage):
    """
    Stage 8/9 for the wavelength solution
    """
    def __init__(self, runtime_context=None):
        super(ApplyToSpectrum, self).__init__(runtime_context=runtime_context)

    def do_stage_fiber(self, image, fiber):
        spectrum = image.data_tables[self.runtime_context.main_spectrum_name]
        fiber_mask = np.where(spectrum['fiber'] == fiber)
        wcs = image.wavelength_solution[fiber]
        pixel_coordinates, order_coordinates = pixel_order_as_array(spectrum[fiber_mask])
        spectrum['wavelength'][fiber_mask] = wcs(pixel=pixel_coordinates, order=order_coordinates)
        spectrum.meta['header'] = {'MODEL': str(wcs.model), 'MCOEFFS': str(list(wcs.model_coefficients))}
        image.data_tables[self.runtime_context.main_spectrum_name] = spectrum
        return image


class TabulateArcEmissionLines(WavelengthStage):
    """
    Optional stage which will create a table of all identified emission lines and their matches
    in the reference list. This is saved when the image reduction is completed, under the extension LINES.
    """
    def __init__(self, runtime_context=None):
        super(TabulateArcEmissionLines, self).__init__(runtime_context=runtime_context)

    def do_stage(self, image):
        # TODO if wavelength solution fails, then the emission lines are not saved. Make it so that they are
        #  always saved
        valid_fibers = self._valid_fibers(image)
        if len(valid_fibers) > 0:
            lines = self._format_lines(image, valid_fibers)
            for fiber in valid_fibers:
                wcs = image.wavelength_solution[fiber]
                fib = np.where(lines['fiber'] == fiber)
                lines['wavelength'][fib] = wcs(lines['pixel'][fib].data, lines['order'][fib].data)
                lines['reference_wavelength'][fib] = find_nearest(lines['wavelength'][fib], wcs.reference_lines)

            image.data_tables[self.runtime_context.emission_lines_table_name] = Table(lines)
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
    number. The principle order number should then be fixed in the config.ini file for your instrument.
    """
    def __init__(self, runtime_context=None):
        super(IdentifyPrincipleOrderNumber, self).__init__(runtime_context=runtime_context)
        self.STAGES_TODO = [SolveFromOverlaps, FindGlobalScale, SolutionRefineInitial, SolutionRefineFinal]

    def do_stage_fiber(self, image, fiber):
        self.logger.info('Looking for the principle order number between {0} and {1}.'
                         ' fiber={2}'.format(*self.runtime_context.m0_range, str(fiber)))
        self.logger.disabled = True
        merits, m0_values = self.merit_per_m0(image, fiber, self.runtime_context.m0_range)
        self.logger.disabled = False
        best_m0, merit = m0_values[np.argmin(merits)], np.min(merits)

        if not merit < 1/10 * np.median(merits):
            self.logger.warning('A definitive principle order number was not found. Aborting wavelength solution.'
                           ' fiber={0}'.format(str(fiber)))
            image.wavelength_solution[fiber] = None
        else:
            image.wavelength_solution[fiber].m0 = best_m0
            self.logger.info('The best principle order number is {0}. fiber={1}'.format(best_m0, str(fiber)))
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


def find_feature_wavelengths(measured_lines, reference_lines, m0_range: tuple, max_pixel, min_pixel=0,
                             wavelength_models: dict = None, overlap_settings: dict = None, scale_settings: dict = None,
                             stages_todo: list = None):
    """
    A convenience function for any users who want to be able to wavelength calibrate from a list of
    spectral feature (measured_lines) pixel and order positions and a reference line list (reference_lines).
    This is not used anywhere in xwavecal, nor the included pipeline, except in
    xwavecal.tests.test_wavelength.

    The measured_lines should be for one fiber only. This function should be looped over fibers to calibrate
    many fibers.

    This function was designed to be used in Banzai-NRES. It returns the wavelengths of the spectral features
    input in measured_lines.

    :param measured_lines: dict.
           dictionary of ndarrays of pixel and order positions of measured spectral features (e.g. emission lines).
           Example:
               measured_lines = {'pixel': np.array([1, 2.5, 6.1]), 'order': np.array([1, 1, 2])}
               If the principle order number is 52, then these measured_lines represents 3 spectral features,
               with (pixel, diffraction order) coordinates of (1, 53), (2.5, 53), (6.1, 54).
               respectively. The wavelength solution will calibrate each fiber separately.
    :param reference_lines: list or ndarray. List of reference wavelengths in Angstroms
           for what features you expect are in measured_lines.
           Example:
               reference_lines = [5000, 5001, 5002, 5003.1, 6001.2]
    :param m0_range: tuple.
           The range of integers possible for the principle order number: (low, high). low is an inclusive bound
           and high is exclusive.
           The principle order number gives the real diffraction order index of the i=0 labelled diffraction order.
           See the xwavecal README or Brandt et al. (2019) for more.
    :param max_pixel: int. Maximum pixel coordinate on the detector. Used for finding the overlaps. E.g. 4096
    :param min_pixel: int. Minimum pixel coordinate on the detector. Used for finding the overlaps. E.g. 0
    :param wavelength_models: dict.
           Dictionary of models. See the xwavecal README, or wavelength_models below.
    :param overlap_settings: dict.
           Dictionary of settings for the overlap fitting algorithms.
           See the xwavecal README, or overlap_settings below.
    :param scale_settings: dict.
           Dictionary of settings for the global scale fitting algorithms.
           See the xwavecal README, or scale_settings below.
    :param stages_todo: list.
           list of xwavecal.wavelength.WavelengthSolution stages to run on the input list of measured lines.
           If building a new wavelength solution from scratch, then leave the default stages_todo.
    :return: measured_line_wavelengths: ndarray. if return_all==False.
             measured_line_wavelengths are the wavelengths of the measured_lines under the calibrated model.
             measured_line_wavelengths is an array with the same length as measured_lines['pixel']
              and measured_lines['order'].
    Notes
    -----
    See xwavecal.tests.test_wavelength.TestOnSyntheticData.test_performance for a use example of this function.
     the test can be run from an ipython instance and can provide insight as to how to use this wrapper function.
    """
    if stages_todo is None:
        stages_todo = [FitOverlaps, IdentifyPrincipleOrderNumber, SolveFromOverlaps,
                       FindGlobalScale, SolutionRefineInitial, SolutionRefineFinal]
    if wavelength_models is None:
        wavelength_models = {'initial_wavelength_model': {1: [0, 1, 2], 2: [0, 1, 2]},
                             'intermediate_wavelength_model': {0: [0, 1, 2], 1: [0, 1, 2], 2: [0, 1, 2]},
                             'final_wavelength_model': {0: [0, 1, 2, 3, 4, 5],
                                                        1: [0, 1, 2, 3, 4, 5],
                                                        2: [0, 1, 2, 3, 4, 5],
                                                        3: [0, 1, 2, 3, 4, 5],
                                                        4: [0]}}
    if overlap_settings is None:
        # settings relevant to fitting the overlaps.
        overlap_settings = {'min_num_overlaps': 5, 'overlap_linear_scale_range': (0.5, 2),
                            'flux_tol': 0.2, 'max_red_overlap': 1000, 'max_blue_overlap': 2000}
    if scale_settings is None:
        # settings relevant to finding the global scale.
        scale_settings = {'global_scale_range': (0.8, 1.2), 'approx_detector_range_angstroms': 5000,
                          'approx_num_orders': 67}
    # instantiate the context which every WavelengthStage will use.
    context = type('context', (), {**wavelength_models, **overlap_settings, **scale_settings,
                                   'm0_range': m0_range, 'main_spectrum_name':
                                   'spectrum', 'overlap_table_name': 'overlap'})

    # make dummy spectrum so that FitOverlaps will run.
    # TODO: the spectrum is only used in fit_overlaps to get the reference id's etc. It in principle is not
    #  needed at all. Fix this deeper in the code so that we don't have to make a fake spectrum here.
    pixel = np.arange(min_pixel, max_pixel + 1)
    orders = np.arange(np.min(measured_lines['order']), np.max(measured_lines['order'])).astype(int)
    spectrum_shape = (len(orders), len(pixel))
    spectrum = Table({'ref_id': orders, 'fiber': np.ones_like(orders),
                      'flux': np.zeros(spectrum_shape), 'pixel': pixel * np.ones(spectrum_shape)})
    # Initialize the WavelengthSolution
    wavelength_solution = WavelengthSolution(model=wavelength_models.get('initial_wavelength_model'),
                                             min_order=np.min(orders), max_order=np.max(orders),
                                             min_pixel=np.min(pixel), max_pixel=np.max(pixel),
                                             measured_lines=measured_lines,
                                             reference_lines=np.sort(reference_lines))
    # make a container (e.g. the Image object) for the spectrum and wavelength solution
    image = Image(header={'fiber_state': 'none&thar&none'},
                  wavelength_solution={1: wavelength_solution},
                  data_tables={context.main_spectrum_name: spectrum})
    # run the wavelength calibration stages
    for stage in stages_todo:
        image = stage(context).do_stage(image)
    measured_lines['wavelength'] = np.zeros_like(measured_lines['pixel'], dtype=float)
    wavelengths = None
    if image.wavelength_solution[1] is not None:
        wavelengths = image.wavelength_solution[1](measured_lines['pixel'], measured_lines['order'])
    measured_lines['wavelength'] = wavelengths
    return measured_lines['wavelength']

