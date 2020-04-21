import pytest
import mock
import numpy as np
from scipy.optimize import OptimizeResult
import copy
from astropy.table import Table

import xwavecal.utils.wavelength_utils as wcsu
from xwavecal.utils.trace_utils import legendre as legendre_val
from xwavecal.utils.overlap_utils import blank_overlap_table
from xwavecal.wavelength import WavelengthSolution, FindGlobalScale, SolutionRefineInitial, SolutionRefineFinal
from xwavecal.wavelength import refine_wcs, FitOverlaps, WavelengthStage, SolveFromOverlaps, IdentifyArcEmissionLines
from xwavecal.wavelength import ApplyToSpectrum, TabulateArcEmissionLines, BlazeCorrectArcEmissionLines, Initialize
from xwavecal.wavelength import IdentifyPrincipleOrderNumber, SolutionRefineOnce, find_feature_wavelengths
from xwavecal.tests.utils import SpectrumUtils, FakeImage, FakeContext


class TestWavelengthSolution:
    def test_format_overlaps(self):
        overlaps = Table({'pixel': [np.arange(3, 6), np.arange(3, 6)], 'id': [1, 2]})
        wcs = WavelengthSolution(min_order=1, max_order=2, min_pixel=3, max_pixel=5)
        coordinates = wcs._format_overlaps(overlaps, pixel_key='pixel', order_key='id')
        assert np.allclose(coordinates['normed_order'], [-1, -1, -1, 1, 1, 1])
        assert np.allclose(coordinates['normed_pixel'], [-1, 0, 1, -1, 0, 1])
        assert np.allclose(coordinates['order'], [1, 1, 1, 2, 2, 2])

    def test_format_overlaps_rejects_nan(self):
        overlaps = Table({'pixel': [np.array([3, 4, np.nan]), np.array([np.nan, 4, 5])], 'id': [1, 2]})
        wcs = WavelengthSolution(min_order=1, max_order=2, min_pixel=3, max_pixel=5)
        coordinates = wcs._format_overlaps(overlaps, pixel_key='pixel', order_key='id')
        assert np.allclose(coordinates['normed_order'], [-1, -1, 1, 1])
        assert np.allclose(coordinates['normed_pixel'], [-1, 0, 0, 1])
        assert np.allclose(coordinates['order'], [1, 1, 2, 2])

    @pytest.mark.integration
    @mock.patch('xwavecal.wavelength.WavelengthSolution._format_overlaps')
    def test_solve_from_overlaps(self, mock_format):
        m0 = 52
        m_coordinates = {'normed_pixel': np.array([1]), 'normed_order': np.array([1]), 'order': np.array([21])}
        coordinates = {'normed_pixel': np.array([-1]), 'normed_order': np.array([-1]), 'order': np.array([20])}

        def fake_format(overlaps, pixel_key=None, order_key=None):
            return m_coordinates if pixel_key == 'matched_pixel' else coordinates
        mock_format.side_effect = fake_format
        wcs = WavelengthSolution(m0=m0, model={1: [0]})  # (1 + b * x)/(m0 + order_num) model

        wcs.solve_from_overlaps(None)
        coefficients = wcs.model_coefficients
        expected_coefficient = 1/(legendre_val(1, 1)/(m0 + 21) - legendre_val(1, -1)/(m0 + 20)) *\
                                                                 (1/(m0 + 20) - 1/(m0 + 21))
        assert np.isclose(coefficients[0], expected_coefficient)

    def test_invalid_model_raises_error(self):
        with pytest.raises(ValueError):
            WavelengthSolution(model={0: [0]}).solve_from_overlaps(None)

    def test_construct_wavelength_map_matrices_grating_eq(self):
        m0 = 52
        coordinates = {'normed_pixel': np.array([-1]), 'normed_order': np.array([-1]), 'order': np.array([20])}
        wcs = WavelengthSolution(m0=m0, model={1: [0]})  # (1 + b * x)/(m0 + order_num) model
        A, c = wcs._construct_wavelength_map_matrices(coordinates['normed_pixel'], coordinates['normed_order'],
                                                      coordinates['order'])
        assert np.allclose(c, [[-1/(m0 + 20)]])
        assert np.allclose(A, [[1/(m0 + 20) * coordinates['normed_pixel'][0]]])

    def test_construct_wavelength_map_matrices_no_grating(self):
        m0 = 52
        coordinates = {'normed_pixel': np.array([-1]), 'normed_order': np.array([-1]), 'order': np.array([20])}
        wcs = WavelengthSolution(m0=m0, model={1: [0]}, grating_eq=False)  # (1 + b * x) model
        A, c = wcs._construct_wavelength_map_matrices(coordinates['normed_pixel'], coordinates['normed_order'],
                                                      coordinates['order'])
        assert np.allclose(c, [[-1]])
        assert np.allclose(A, [[1 * coordinates['normed_pixel'][0]]])

    def test_construct_wavelength_map_matrices_no_force_unity(self):
        m0 = 52
        coordinates = {'normed_pixel': np.array([-1]), 'normed_order': np.array([-1]), 'order': np.array([20])}
        wcs = WavelengthSolution(m0=m0, model={0: [0]})  # (1 + b * x)/(m0 + order_num) model
        A, c = wcs._construct_wavelength_map_matrices(coordinates['normed_pixel'], coordinates['normed_order'],
                                                      coordinates['order'])
        assert np.allclose(c, [[0]])

    def test_solve(self):
        coordinates = {'normed_pixel': np.array([-1, 1]), 'normed_order': np.array([1, 1]),
                       'order': np.array([20, 20])}
        wcs = WavelengthSolution(m0=52, model={0: [0], 1: [1]})  # (a + b*x*i)/(m0+i) model
        a, b = 1.31, 5.42
        wavelengths = 1/(52 + 20) * np.array([a + b * (-1)*1, a + b * 1*1])
        coeffs = wcs.solve(coordinates, wavelengths)
        assert np.allclose(coeffs, [a, b])

    def test_weighted_solve(self):
        coordinates = {'normed_pixel': np.array([-1, 1, .5]), 'normed_order': np.array([1, 1, .5]),
                       'order': np.array([20, 20, 22])}
        wcs = WavelengthSolution(m0=52, model={0: [0], 1: [1]})  # (a + b*x*i)/(m0+i) model
        a, b = 1.31, 5.42
        wavelengths = 1/(52 + 20) * np.array([a + b * (-1)*1, a + b * 1*1, 10])
        coeffs = wcs.solve(coordinates, wavelengths, weights=np.array([1, 1, 0]))
        assert np.allclose(coeffs, [a, b])

    def test_update_model(self):
        wcs = WavelengthSolution(m0=52, model={0: [0], 1: [1, 2]}, model_coefficients=np.arange(1, 4),
                                 min_order=0, max_order=5, min_pixel=0, max_pixel=400)
        new_model = {0: [0, 1], 1: [1, 2, 3]}
        wcs.update_model(new_model=new_model)
        assert new_model == wcs.model
        assert np.allclose(np.array([1, 0, 2, 3, 0]), wcs.model_coefficients)

    def test_update_model_complex(self):
        model = {0: [0, 1, 2, 3, 4, 5], 1: [0, 1, 2, 3, 4, 5], 2: [0, 1, 2, 3, 4, 5], 3: [0, 1, 2, 3, 4]}
        wcs = WavelengthSolution(m0=52, model=model, model_coefficients=np.ones(23),
                                 min_order=0, max_order=67, min_pixel=0, max_pixel=4095)
        new_model = wcsu.Model(copy.deepcopy(model))
        new_model.add_polynomial_term(3, 5)
        wcs.update_model(new_model=new_model)
        assert new_model == wcs.model
        assert np.allclose([1]*23 + [0], wcs.model_coefficients)

    def test_apply_scale(self):
        wcs = WavelengthSolution(model_coefficients=np.arange(5))
        wcs.apply_scale(2)
        assert (wcs.model_coefficients, np.arange(5) * 2)

    def test_wavelength_is_nan_for_no_coefficients(self):
        wcs = Utils.simple_wcs()
        wcs.model_coefficients = None
        assert np.all(np.isnan(wcs(np.arange(5), np.arange(5))))


class TestModel:
    def test_add_term(self):
        model = wcsu.Model({0: [0]})
        model.add_polynomial_term(x_degree=3, i_degree=2)
        assert model == wcsu.Model({0: [0], 3: [2]})

    def test_check_term_exists(self):
        assert not wcsu.Model({0: [1]}).is_missing_polynomial_term(0, 1)
        assert wcsu.Model({0: [1]}).is_missing_polynomial_term(0, 0)


class TestStatistics:
    @mock.patch('xwavecal.utils.wavelength_utils.find_nearest', return_value=np.array([1, 2.]))
    def test_calc_residuals(self, fake_match):
        assert np.allclose(wcsu.calc_residuals(np.array([1., 3]), None), [0, 1.0])
        assert np.allclose(wcsu.calc_residuals(np.array([1., 2]), None), [0, 0])


class TestWavelengthStage:
    def test_do_stage(self):
        image = FakeImage()
        image.fiber0_wavecal, image.fiber1_wavecal, image.fiber2_wavecal = 0, 1, 0
        image.wavelength_solution = {0: None, 1: 'placeholder', 2: None}
        image = WavelengthStage(FakeContext()).do_stage(image)
        assert image.wavelength_solution[1] == 'placeholder'

    @mock.patch('xwavecal.wavelength.WavelengthStage.do_stage_fiber', return_value=None)
    def test_do_stage_skips_if_no_valid_fibers(self, fake_fiber_stage):
        image = FakeImage()
        image.fiber0_wavecal, image.fiber1_wavecal, image.fiber2_wavecal = 0, 0, 1
        image.wavelength_solution = {0: None, 1: 'placeholder', 2: None}
        assert WavelengthStage(FakeContext()).do_stage(image) is not None

    def test_valid_fibers(self):
        image = FakeImage()
        image.fiber0_wavecal, image.fiber1_wavecal, image.fiber2_wavecal = 0, 1, 0
        image.wavelength_solution = {0: None, 1: 'placeholder', 2: None}
        valid_fibers = WavelengthStage._valid_fibers(image)
        assert np.allclose(valid_fibers, [1])
        image.wavelength_solution = {0: None, 1: None, 2: 'placeholder'}
        assert len(WavelengthStage._valid_fibers(image)) == 0


class TestInitialize:
    CONTEXT = FakeContext()

    def test_on_no_valid_fibers(self):
        image = FakeImage()
        image.fiber0_wavecal, image.fiber1_wavecal, image.fiber2_wavecal = 0, 1, 1
        image.wavelength_solution = {1: 'not none', 2: 'not none'}
        Initialize(None).on_no_valid_fibers(image)
        assert image.wavelength_solution == {1: None, 2: None}

    def test_all_fibers_invalid_without_spectrum(self):
        image = FakeImage()
        image.data_tables = {}
        assert len(Initialize(self.CONTEXT)._valid_fibers(image)) == 0

    def test_all_fibers_invalid_without_ref_id_and_fiber_cols(self):
        image = FakeImage()
        image.data_tables[self.CONTEXT.main_spectrum_name] = Table({'col1': [1]})
        assert len(Initialize(self.CONTEXT)._valid_fibers(image)) == 0


class TestFitOverlaps:
    @mock.patch('xwavecal.wavelength.fit_overlaps')
    def test_do_stage(self, mock_overlaps):
        num = 2
        image = FakeImage()
        context = FakeContext()
        spectrum = Table({'ref_id': [0, 1], 'fiber': [0, 1], 'flux': [[0, 1, 2], [0, 1, 2]],
                          'pixel': [[-1, 0, 1], [-1, 0, 1]]})
        lines = Table({'order': [0, 1], 'fiber': [0, 1], 'flux': [[0, 1, 2], [0, 1, 2]],
                       'pixel': [[-1, 0, 1], [-1, 0, 1]]})
        mock_overlaps.return_value = Table({'fiber': [1, 0], 'ref_id': [0, 1], 'matched_ref_id': [1, 2],
                                            'good': [True, True], 'peaks': [10, 10]})
        image.data_tables = {context.main_spectrum_name: spectrum}
        wcs = WavelengthSolution(min_order=0, max_order=num, min_pixel=0, max_pixel=500, model={1: [1], 2: [1]})
        image.wavelength_solution = {0: wcs, 1: wcs}
        image.fiber0_wavecal, image.fiber1_wavecal, image.fiber2_wavecal = 1, 1, 0
        image = FitOverlaps(FakeContext()).do_stage(image)
        assert isinstance(image.data_tables[context.overlap_table_name], Table)

    @mock.patch('xwavecal.wavelength.fit_overlaps')
    def test_do_stage_single_fiber_stacks_overlap_tables(self, fake_fit):
        image = FakeImage()
        context = FakeContext()
        image.data_tables = {context.overlap_table_name: Table({'ref_id': [0], 'fiber': [0], 'matched_ref_id': [1],
                                                                      'good': [False], 'peaks': [10]})}
        image.data_tables[context.main_spectrum_name] = Table({'fiber': [0, 1]})
        wcs = WavelengthSolution(min_order=0, max_order=10, min_pixel=0, max_pixel=500, model={1: [1], 2: [1]})
        image.wavelength_solution = {0: wcs, 1: wcs}
        fake_fit.return_value = Table({'ref_id': [0], 'fiber': [1], 'matched_ref_id': [1], 'good': [False], 'peaks': [2]})
        image = FitOverlaps(FakeContext()).do_stage_fiber(image, fiber=1)
        assert len(image.data_tables[context.overlap_table_name]['fiber']) == 2

    @mock.patch('xwavecal.wavelength.fit_overlaps', return_value=blank_overlap_table(3))
    def test_do_stage_sets_wcs_to_none(self, fake_overlaps):
        image = FakeImage()
        context = FakeContext()
        image.data_tables = {context.main_spectrum_name: Table({'ref_id': [0, 1],
                                                                     'fiber': [0, 1],
                                                                     'flux': [[0, 1, 2], [0, 1, 2]],
                                                                     'pixel': [[-1, 0, 1], [-1, 0, 1]]}),
                             context.overlap_table_name: blank_overlap_table(3)}
        image.wavelength_solution = {0: WavelengthSolution()}
        image = FitOverlaps(FakeContext()).do_stage_fiber(image, fiber=0)
        assert image.wavelength_solution[0] is None


class TestSolveFromOverlaps:
    @mock.patch('xwavecal.wavelength.WavelengthSolution.solve_from_overlaps')
    def test_do_stage_does_not_crash(self, fake_solve):
        image = FakeImage()
        context = FakeContext()
        context.min_num_overlaps = 0
        image.data_tables = {context.overlap_table_name: Table({'ref_id': [0, 1],
                                                                     'matched_ref_id': [1, 2],
                                                                     'good': [True, True],
                                                                     'fiber': [0, 0]})}
        image.wavelength_solution = {0: WavelengthSolution()}
        image = SolveFromOverlaps(FakeContext()).do_stage_fiber(image, fiber=0)
        assert np.allclose(image.wavelength_solution[0].overlap_range, (0, 2))

    @mock.patch('xwavecal.wavelength.WavelengthSolution.solve_from_overlaps')
    def test_do_stage_does_not_overwrite_overlaps(self, fake_solve):
        image = FakeImage()
        context = FakeContext()
        context.min_num_overlaps = 0
        image.data_tables = {context.overlap_table_name: Table({'fiber': [1, 0],
                                                                     'ref_id': [0, 1],
                                                                     'matched_ref_id': [1, 2],
                                                                     'good': [True, True]})}
        image.fiber0_wavecal, image.fiber1_wavecal, image.fiber2_wavecal = 1, 1, 0
        image.wavelength_solution = {0: Utils.simple_wcs(), 1: Utils.simple_wcs()}
        image = SolveFromOverlaps(FakeContext()).do_stage(image)
        assert np.allclose(image.data_tables[context.overlap_table_name]['fiber'], [1, 0])


class TestIdentifyArcEmissionLines:
    @pytest.mark.integration
    @mock.patch('xwavecal.wavelength.identify_lines')
    def test_do_stage_fiber(self, fake_measured_lines):
        fake_measured_lines.return_value = {'pixel': np.array([0, 20]), 'order': np.array([0, 5]), 'pixel_err': np.array([1, 2])}
        image = FakeImage()
        image.wavelength_solution = {1: WavelengthSolution(min_pixel=0, max_pixel=20, min_order=-5,
                                                           max_order=5)}
        image.data_tables = {FakeContext().main_spectrum_name: Table({'fiber': np.array([0, 1]),
                                                                     'stderr': [0, 0]})}
        stage = IdentifyArcEmissionLines(FakeContext())
        image = stage.do_stage_fiber(image, fiber=1)
        measured_lines = image.wavelength_solution[1].measured_lines
        assert np.allclose([measured_lines['pixel'], measured_lines['order']], [[0, 20], [0, 5]])
        assert np.allclose([measured_lines['normed_pixel'], measured_lines['normed_order']],
                           [[-1, 1], [0, 1]])


class TestFindGlobalScale:
    @mock.patch('xwavecal.wavelength.WavelengthSolution.update_model')
    @mock.patch('xwavecal.wavelength.FindGlobalScale._find_scale', return_value=10)
    @mock.patch('xwavecal.wavelength.estimate_global_scale', return_value=1)
    def test_do_stage_fiber(self, fake_guess, fake_result, fake_change):
        image = FakeImage()
        image.wavelength_solution = {1: WavelengthSolution(model_coefficients=1)}
        image = FindGlobalScale(FakeContext()).do_stage_fiber(image, fiber=1)
        assert np.isclose(image.wavelength_solution[1].model_coefficients, 10)

    def test_estimate_global_scale(self):
        assert np.isclose(wcsu.estimate_global_scale(5000, n=67, m0=52), 461791, atol=2)
        assert np.isclose(wcsu.estimate_global_scale(6000, n=67, m0=50), 523880, atol=2)

    @mock.patch('xwavecal.wavelength.calc_residuals', return_value=np.arange(10))
    def test_chi_squared_many_scales(self, fake_chi_squared):
        chi_squared = FindGlobalScale._chi_squared(scales=np.array([1, 2]), unscaled_wavelengths=np.ones(5),
                                                   reference_wavelengths=None)
        assert chi_squared.shape == (2,)
        assert np.allclose(chi_squared[0], np.sum(np.arange(5) ** 2))
        assert np.allclose(chi_squared[1], np.sum(np.arange(5, 10) ** 2))

    def test_chi_squared_safe(self):
        ref = 5 * np.random.random(30)
        chi_squared_s = FindGlobalScale._chi_squared_safe(scales=np.arange(5), unscaled_wavelengths=np.ones(5),
                                                          reference_wavelengths=ref, mem_limit=32)
        chi_squared = FindGlobalScale._chi_squared(scales=np.arange(5), unscaled_wavelengths=np.ones(5),
                                                   reference_wavelengths=ref)
        assert np.allclose(chi_squared, chi_squared_s)

    @mock.patch('xwavecal.wavelength.calc_residuals', return_value=np.arange(5))
    def test_chi_squared_single_scale(self, fake_chi_squared):
        chi_squared = FindGlobalScale._chi_squared(scales=1, unscaled_wavelengths=np.ones(5),
                                                   reference_wavelengths=None)
        assert chi_squared.shape == (1,)
        assert np.allclose(chi_squared[0], np.sum(np.arange(5) ** 2))

    @mock.patch('xwavecal.wavelength.WavelengthSolution.wavelength_normed_input', return_value=1)
    @mock.patch('xwavecal.wavelength.restrict')
    @mock.patch('xwavecal.wavelength.optimize.minimize')
    def test_find_scale(self, fake_minimize, fake_restrict, fake_wavelength):
        wcs = WavelengthSolution(reference_lines=np.arange(10)[::-1])

        def check_minimize_args(fun=None, args=(), x0=None, method=None, options=None, **kwargs):
            args_passed = (x0 == 'bla' and np.allclose(args[1], np.arange(10)) and args[0] == 1)
            return OptimizeResult(x=[args_passed])
        fake_minimize.side_effect = check_minimize_args
        assert FindGlobalScale._find_scale(wcs, scale_guess='bla')

    @staticmethod
    def random_wcs_with_random_reference_lines(num_lines=200, cntd_error=0.0, true_scale=1E5, ml=None):
        # note centroid error cntd_error is the error in wavelength space, not pixel space.
        m0 = np.random.randint(10, 100)
        true_scale = np.random.normal(true_scale, true_scale/10)
        pix = np.random.randint(10, 1000, size=num_lines)
        orders = np.random.randint(0, 4, size=num_lines)
        normed_pix, normed_orders = wcsu.normalize_coordinates(pix, 1000, 0), wcsu.normalize_coordinates(orders, 4, 0)
        measured_lines = {'normed_pixel': normed_pix[ml:], 'normed_order': normed_orders[ml:], 'order': orders[ml:]}
        wcs = WavelengthSolution(min_order=0, max_order=4, min_pixel=0, max_pixel=1000,
                                 measured_lines=measured_lines, model={0: [0, 1, 2], 1: [0, 1]},
                                 model_coefficients=np.random.normal(2, .3, size=5), m0=m0)
        reference_lines = wcs(pix, orders) * true_scale
        wcs.reference_lines = reference_lines + np.random.normal(0, cntd_error)  # centroiding errors, sort of.
        return wcs, true_scale

    @pytest.mark.integration
    def test_accuracy_best_case(self):
        wcs, true_scale = self.random_wcs_with_random_reference_lines(num_lines=20, cntd_error=0.0, true_scale=1E5)
        scale = FindGlobalScale._find_scale(wcs, scale_guess=1.1*true_scale, rrange=(0.1, 1.5))
        assert np.isclose(scale, true_scale)

    @pytest.mark.integration
    def test_accuracy_with_centroiding_error(self):
        err = 10
        wcs, true_scale = self.random_wcs_with_random_reference_lines(num_lines=80, cntd_error=err, true_scale=1E5)
        scale = FindGlobalScale._find_scale(wcs, scale_guess=1.1*true_scale, rrange=(0.1, 1.5))
        assert np.isclose(scale, true_scale, atol=100 * err)

    @pytest.mark.integration
    def test_accuracy_with_missed_lines(self):
        err = 10
        missed_lines = 40
        wcs, true_scale = self.random_wcs_with_random_reference_lines(num_lines=80, cntd_error=err,
                                                                      true_scale=1E5, ml=missed_lines)
        scale = FindGlobalScale._find_scale(wcs, scale_guess=1.1*true_scale, rrange=(0.1, 1.5))
        assert np.isclose(scale, true_scale, rtol=0.01)


class TestSolutionRefineInitial:
    @mock.patch('xwavecal.wavelength.SolutionRefineInitial.constrain_solution_over_detector')
    def test_do_stage_fiber(self, fake_constrain):
        image = FakeImage()
        context = FakeContext()
        context.intermediate_wavelength_model = {1: [0], 2: [0, 1]}
        image.wavelength_solution = {1: WavelengthSolution(model_coefficients=np.arange(2), model={1: [0], 2: [0]},
                                                           min_order=0, max_order=5, min_pixel=0, max_pixel=400)}
        def constrain(wcs, *args, **kwargs):
            return wcs, 1
        fake_constrain.side_effect = constrain
        image = SolutionRefineInitial(context).do_stage_fiber(image, fiber=1)
        assert image.wavelength_solution[1].model == {1: [0], 2: [0, 1]}

    def test_selecting_orders(self):
        wcs = WavelengthSolution(min_order=0, max_order=7, overlap_range=(2, 4))
        order_ranges = SolutionRefineInitial._ranges_to_evaluate(*wcs.overlap_range, wcs.min_order, wcs.max_order)
        assert np.allclose(order_ranges, [(2, 4), (1, 5), (0, 6), (-1, 7)])

    @pytest.mark.integration
    def test_accuracy_constrain_solution_over_detector(self):
        m0 = np.random.randint(10, 100)
        wcs = Utils.simple_wcs(m0)
        wcs, residuals = SolutionRefineInitial.constrain_solution_over_detector(wcs)
        assert np.isclose(np.mean(residuals), 0, atol=1E-5)

    def test_convergence_criterion(self):
        assert not SolutionRefineInitial._converged_when_max_iter()


class TestSolutionRefineFinal:
    @mock.patch('xwavecal.wavelength.SolutionRefineFinal._refine')
    def test_do_stage_fiber(self, fake_refine):
        image = FakeImage()
        context = FakeContext()
        context.final_wavelength_model = {1: [0], 2: [0, 1]}
        image.wavelength_solution = {1: WavelengthSolution(model={1: [0], 2: [0]})}

        def refine(wcs, final_model, *args, **kwargs):
            wcs.model = final_model
            return wcs, 1
        fake_refine.side_effect = refine
        image = SolutionRefineFinal(context).do_stage_fiber(image, fiber=1)
        assert image.wavelength_solution[1].model == {1: [0], 2: [0, 1]}

    @mock.patch('xwavecal.wavelength.refine_wcs')
    def test_refine(self, fake_refine):
        def refine(wcs, *args, **kwargs):
            wcs.models_tried.append(wcs.model)
            return wcs, 'residuals_passed'
        fake_refine.side_effect = refine
        wcs = WavelengthSolution(min_order=1, max_order=3, min_pixel=0, max_pixel=30,
                                 model={0: [0], 1: [0]},
                                 model_coefficients=np.array([1.0, 1.0]),
                                 reference_lines=np.array([1, 2]), m0=52, grating_eq=True)
        wcs.models_tried = []
        final_model = {0: [0, 1], 1: [0, 1]}
        wcs, residuals = SolutionRefineFinal._refine(wcs, final_model)
        assert residuals is 'residuals_passed'
        assert wcs.models_tried[0] == {0: [0, 1], 1: [0]}
        assert wcs.models_tried[1] == final_model

    def test_convergence(self):
        converged = False
        for residuals in [[1, 2, 3], [1, 2, 2], [2, 2, 2], [10, 2, 2]]:
            if SolutionRefineFinal._converged(residuals, [2, 2, 2]):
                converged = True
                break
        assert converged


class TestSolutionRefineOnce:
    @mock.patch('xwavecal.wavelength.refine_wcs')
    def test_do_stage_fiber(self, fake_refine):
        image = FakeImage()
        context = FakeContext()
        image.wavelength_solution = {1: WavelengthSolution(model={1: [0], 2: [0]})}

        def refine(wcs, *args, **kwargs):
            wcs.model = {1: [0], 2: [0, 1]}
            return wcs, 1
        fake_refine.side_effect = refine
        image = SolutionRefineOnce(context).do_stage_fiber(image, fiber=1)
        assert image.wavelength_solution[1].model == {1: [0], 2: [0, 1]}


class TestApplyToSpectrum:
    @pytest.mark.integration
    def test_do_stage(self):
        image = FakeImage()
        context = FakeContext()
        image.data_tables = {context.main_spectrum_name: Table({'ref_id': [0, 1],
                                                                     'fiber': [1, 2],
                                                                     'flux': [[0, 1, 2], [0, 1, 2]],
                                                                     'pixel': [[-1, 0, 1], [-1, 0, 1]],
                                                                     'wavelength': np.zeros((2, 3))})}
        image.fiber0_wavecal, image.fiber1_wavecal, image.fiber2_wavecal = 0, 1, 1
        image.wavelength_solution = {0: None, 1: Utils.simple_wcs(), 2: Utils.simple_wcs()}
        image = ApplyToSpectrum(context).do_stage(image)
        spectrum = image.data_tables[context.main_spectrum_name]
        assert not np.allclose(spectrum['wavelength'], 0)
        assert len(spectrum.colnames) == 5


class TestBlazeCorrectArcEmissionLines:
    def test_fluxes_are_assigned(self):
        image = FakeImage()
        context = FakeContext()
        image.data_tables = {context.blaze_corrected_spectrum_name: Table({'ref_id': [0, 1, 0],
                                                                     'fiber': [1, 1, 0],
                                                                     'flux': [[0, 10, 2], [0, 20, 2], [0, 15, 2]],
                                                                     'pixel': [[0, 1, 2], [0, 1, 2], [0, 1, 2]]})}
        image.fiber0_wavecal, image.fiber1_wavecal, image.fiber2_wavecal = 0, 1, 0
        image.wavelength_solution = {1: WavelengthSolution()}
        image.wavelength_solution[1].measured_lines = {'pixel': np.array([1, 1]), 'order': np.array([0, 1]),
                                                       'flux': np.array([-9, -8]),
                                                       'corrected_flux': np.array([80, 60])}
        image = BlazeCorrectArcEmissionLines(context).do_stage(image)
        assert np.allclose(image.wavelength_solution[1].measured_lines['corrected_flux'], [10, 20])

    def test_fluxes_are_assigned_on_missing_blaze_correction(self):
        image = FakeImage()
        context = FakeContext()
        image.data_tables = {}
        image.fiber0_wavecal, image.fiber1_wavecal, image.fiber2_wavecal = 0, 1, 0
        image.wavelength_solution = {1: WavelengthSolution()}
        image.wavelength_solution[1].measured_lines = {'pixel': np.array([1, 1]), 'order': np.array([0, 1]),
                                                       'flux': np.array([-9, -8])}
        image = BlazeCorrectArcEmissionLines(context).do_stage(image)
        assert np.allclose(image.wavelength_solution[1].measured_lines['corrected_flux'],
                           image.wavelength_solution[1].measured_lines['flux'])


class TestTabulateArcEmissionLines:
    def test_format_lines(self):
        image = FakeImage()
        image.wavelength_solution = {0: WavelengthSolution(measured_lines={'pixel': [0, 1]}),
                                     1: WavelengthSolution(measured_lines={'pixel': [1, 2, 3]})}
        lines = TabulateArcEmissionLines._format_lines(image, fibers=[0, 1])
        assert all([name in lines.colnames for name in ['fiber', 'pixel', 'wavelength', 'reference_wavelength']])
        assert np.allclose(lines['pixel'], [0, 1, 1, 2, 3])
        assert np.allclose(lines['fiber'], [0, 0, 1, 1, 1])

    def test_do_stage(self):
        image = FakeImage()
        context = FakeContext()
        image.data_tables = {}
        image.fiber0_wavecal, image.fiber1_wavecal, image.fiber2_wavecal = 1, 1, 0
        image.wavelength_solution = {0: Utils.simple_wcs(50), 1: Utils.simple_wcs(52)}
        image = TabulateArcEmissionLines(context).do_stage(image)
        lines = image.data_tables[context.emission_lines_table_name]
        length = len(image.wavelength_solution[0].measured_lines['pixel']) * 2
        assert np.allclose([len(lines['fiber']), len(lines['pixel']), len(lines['wavelength'])], length)
        assert np.allclose([len(lines['reference_wavelength']), len(lines['order'])], length)


def test_normalize_coordinates():
    normed = wcsu.normalize_coordinates(np.arange(3, 6), max_value=5, min_value=3)
    assert np.allclose(normed, [-1, 0, 1])


def test_normalize_coordinates_raises_error():
    with pytest.raises(Exception):
        wcsu.normalize_coordinates(np.arange(3, 6), max_value=3, min_value=5)


class TestLineIdentification:
    @mock.patch('xwavecal.utils.wavelength_utils.find_peaks',
                return_value=(np.array([0.1, 1.1]), np.array([0.01, 0.02]), np.array([0, 1])))
    def test_identify_lines(self, fake_peaks):
        spectrum = Table({'flux': [[10, 20]], 'pixel': [[1, 2]], 'id': [32]})
        lines = wcsu.identify_lines(spectrum, stderr=10, order_key='id')
        assert np.allclose(lines['flux'], [10, 20])
        assert np.allclose(lines['order'], [32, 32])
        assert np.allclose(lines['pixel'], [0.1, 1.1])
        assert np.allclose(lines['pixel_err'], [0.01, 0.02])


class TestRefineUtils:
    def test_restrict(self):
        lines = {'order': np.arange(10), 'misc': np.arange(10, 20)}
        r_lines = wcsu.restrict(lines, 2, high_order=8)
        assert np.allclose(r_lines['order'], np.arange(2, 9))
        assert np.allclose(r_lines['misc'], np.arange(12, 19))

    def test_restrict_returns_original(self):
        lines = {'order': np.arange(10), 'misc': np.arange(10, 20)}
        r_lines = wcsu.restrict(lines)
        assert np.allclose(r_lines['order'], lines['order'])
        assert np.allclose(r_lines['misc'], lines['misc'])

    @staticmethod
    def no_clip(wcs, measured_lines, reference_lines, **kwargs):
        return measured_lines

    @staticmethod
    def converge(residuals, *args, **kwargs):
        return True

    def test_refine_wcs_accuracy(self):
        # this fails very rarely, but the occasional failure is annoying.
        # So we run the test twice and look for at least one success.
        success = []
        for i in range(2):
            m0 = np.random.randint(10, 100)
            wcs = Utils.simple_wcs(m0)
            wcs, residuals = refine_wcs(wcs, wcs.measured_lines, wcs.reference_lines, self.converge, self.no_clip)
            success.append(np.isclose(np.mean(residuals), 0, atol=1E-5))
        assert any(success)

    def test_clip_returns_lines_without_residuals(self):
        assert np.allclose([1, 2, 3], wcsu._sigma_clip([], [1, 2, 3]))

    def test_clip(self):
        clipped_lines = wcsu._sigma_clip([10, 10, *np.ones(30)], {'a': np.arange(32)}, sigma=3, stdfunc=np.std)
        assert np.allclose(clipped_lines['a'], np.arange(32)[2:])


class TestIdentifyPrincipleOrderNumber:
    CONTEXT = FakeContext()

    @mock.patch('xwavecal.wavelength.IdentifyPrincipleOrderNumber.merit_per_m0', return_value=([1, .09, 1], np.array([2, 3, 4])))
    def test_do_stage_fiber(self, mock_merit):
        image = FakeImage()
        image.wavelength_solution = {1: Utils.simple_wcs(m0=1000)}
        image = IdentifyPrincipleOrderNumber(self.CONTEXT).do_stage_fiber(image, 1)
        assert image.wavelength_solution[1].m0 == 3

    @mock.patch('xwavecal.wavelength.IdentifyPrincipleOrderNumber.merit_per_m0', return_value=([1, .5, 1], np.array([2, 3, 4])))
    def test_failure(self, mock_merit):
        image = FakeImage()
        image.wavelength_solution = {1: Utils.simple_wcs(m0=1000)}
        image = IdentifyPrincipleOrderNumber(self.CONTEXT).do_stage_fiber(image, 1)
        assert image.wavelength_solution[1] is None

    def test_merit_per_m0(self):
        image = FakeImage()
        image.wavelength_solution = {1: Utils.simple_wcs(m0=10)}
        identifier = IdentifyPrincipleOrderNumber(self.CONTEXT)
        identifier.STAGES_TODO = []
        merit, m0_vals = identifier.merit_per_m0(image, 1, (5, 11, 5)) # try m0=5 and 10.
        assert merit[1] < merit[0]
        assert np.allclose(m0_vals, [5, 10])


class TestOnSyntheticData:
    @pytest.mark.integration
    def test_find_feature_wavelengths(self):
        # Generate emission line positions from a real wavelength solution on NRES.
        num_orders = 67
        wcs = WavelengthSolution(model={0: [0, 1, 2, 3, 4, 5],
                                        1: [0, 1, 2, 3, 4, 5],
                                        2: [0, 1, 2, 3, 4, 5],
                                        3: [0, 1, 2, 3, 4, 5],
                                        4: [0]},
                                 min_order=0, max_order=num_orders - 1,
                                 min_pixel=0, max_pixel=4095, m0=52,
                                 model_coefficients=np.array([4.66146591e+05, -7.74683540e+02, -3.54722161e+02, -8.99350361e+01,
                                                              -1.70864902e+01, -2.49239167e+00,  4.67592938e+03,  1.23351661e+02,
                                                              6.16440684e+01,  2.03181960e+01,  3.20985475e+00,  1.26399172e-02,
                                                              -5.20952098e+02, -2.52553564e+01, -1.41207709e+01, -1.80775282e+00,
                                                              -8.57769848e-01,  1.75480247e+00,  1.95194420e+01,  3.81273242e+00,
                                                              7.61239901e-01, -8.78661806e-01, -1.73232452e+00, -1.96449284e+00,
                                                              -3.57462413e+00]))

        measured_lines, line_list = SpectrumUtils().generate_measured_lines(n_list=1000, n_true=0,
                                                                            n_garb=0, wcs=wcs)
        measured_lines['fiber'] = np.ones_like(measured_lines['order'])
        wavelength_models = {'initial_wavelength_model': {1: [0, 1, 2], 2: [0, 1, 2]},
                             'intermediate_wavelength_model': {0: [0, 1, 2], 1: [0, 1, 2], 2: [0, 1, 2]},
                             'final_wavelength_model': copy.copy(wcs.model)}
        measured_lines['wavelength'] = find_feature_wavelengths(measured_lines, line_list,
                                                                m0_range=(51, 54), max_pixel=4096, min_pixel=0,
                                                                wavelength_models=wavelength_models)
        assert np.allclose(measured_lines['wavelength'], measured_lines['true_wavelength'], rtol=1e-8)


class Utils:
    @staticmethod
    def simple_wcs(m0=10):
        line_orders, line_normed_pixels = np.array([1., 1, 2, 2, 3, 3]), np.array([-1., 1, -1, 1, -1, 1])
        measured_lines = {'normed_pixel': line_normed_pixels, 'normed_order': np.array([-1., -1, 0, 0, 1, 1]),
                          'order': line_orders, 'pixel': line_normed_pixels}
        reference_lines = 1./(m0 + line_orders) * (1. + line_normed_pixels)
        wcs = WavelengthSolution(min_order=1, max_order=3, overlap_range=(2, 3), measured_lines=measured_lines,
                                 model={0: [0], 1: [0]}, model_coefficients=np.array([1.01, 1.02]),
                                 reference_lines=reference_lines, m0=m0, grating_eq=True)
        return wcs
