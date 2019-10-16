import pytest
import numpy as np

from astropy.table import Table
from astropy.io import fits
from scipy.optimize import OptimizeResult
from unittest import mock

from xwavecal.traces import TraceMaker, LoadTrace
from xwavecal.tests.utils import array_with_peaks, FakeImage, noisify_image, FakeContext
from xwavecal.utils.trace_utils import Trace, SingleTraceFitter, AllTraceFitter
from xwavecal.tests.utils import fill_image_with_traces

import logging
logger = logging.getLogger(__name__)


class FakeTraceImage(FakeImage):
    """
    Image of 500x500 is recommended. Drastic changes to that dimension may break trace testing.
    """
    def __init__(self, nx=500, ny=502, *args, **kwargs):
        super(FakeTraceImage, self).__init__(*args, **kwargs)
        self.caltype = 'TRACE'
        self.header['OBJECTS'] = 'tung&tung&none'
        self.nx = nx
        self.ny = ny
        self.bpm = np.zeros((self.ny, self.nx), dtype=np.uint8)
        self.data = np.zeros((self.ny, self.nx))
        self.fiber0_lit, self.fiber1_lit, self.fiber2_lit = False, True, True


class TestTrace:
    """
    Unit tests for the Trace class.
    """
    def test_trace_instantiates_from_num_centers(self):
        trace = Trace(num_centers_per_trace=5)
        assert trace.data.colnames == ['id', 'centers']
        assert len(trace.data['id']) == 0
        assert trace.data['centers'].shape == (0, 5)
        assert trace.data['centers'].description is not None
        assert trace.data['id'].description is not None

    def test_trace_instantiates_from_data(self):
        data = Table({'id': [1], 'centers': [[1, 2, 3]]})
        data['id'].description = 'test'
        data['centers'].description = 'test_2'
        trace = Trace(data=data)
        assert trace.filepath is None
        for name in ['id', 'centers']:
            assert name in trace.data.colnames
        assert len(trace.data.colnames) == 2
        assert len(trace.data['id']) == 1
        assert trace.data['centers'].shape == (1, 3)
        assert trace.data['centers'].description == 'test_2'
        assert trace.data['id'].description == 'test'

    def test_trace_raises_exception(self):
        with pytest.raises(Exception):
            Trace(data=None, num_centers_per_trace=0)

    def test_getting_trace_centers(self):
        trace = Trace(data={'id': [0, 1], 'centers': [[0, 1], [1, 2]]})
        assert np.allclose(trace.get_centers(0), [0, 1])

    def test_getting_trace_id(self):
        trace = Trace(data={'id': [0, 1], 'centers': [[0, 1], [1, 2]]})
        assert trace.get_id(-1) == 1
        assert trace.get_id(0) == 0

    def test_getting_num_found_traces(self):
        trace = Trace(data={'id': [0, 1], 'centers': [[0, 1], [1, 2]]})
        assert trace.num_traces_found() == 2

    def test_add_centers_to_empty(self):
        centers = np.arange(3)
        trace = Trace(data=None, num_centers_per_trace=len(centers))
        trace.add_centers(trace_centers=centers, trace_id=1)
        assert np.allclose(trace.data['centers'], [centers])
        assert np.allclose(trace.data['id'], [1])

    def test_add_centers_to_existing(self):
        centers = np.arange(3)
        trace = Trace(data={'id': [1], 'centers': [centers]})
        trace.add_centers(trace_centers=centers, trace_id=2)
        assert np.allclose(trace.data['centers'], [centers, centers])
        assert np.allclose(trace.data['id'], [1, 2])

    def test_sorting_trace_centers(self):
        centers = np.array([1, 2, 3])
        data = {'id': [1, 2, 3, 4],
                'centers': [centers, centers+5, centers-10, centers+2]}
        trace = Trace(data=data)
        trace.sort()
        assert np.allclose(trace.data['id'], np.arange(4))
        assert np.allclose(trace.data['centers'],
                           np.array([centers-10, centers, centers+2, centers+5]))

    def test_delete_centers(self):
        centers = np.array([1, 2, 3, 4])
        data = {'id': [1, 2, 3, 4],
                'centers': [centers, centers+5, centers+10, centers+11]}
        trace = Trace(data=data)
        trace.del_centers([])
        assert np.allclose(trace.data['id'], data['id'])
        assert np.allclose(trace.data['centers'], data['centers'])
        trace.del_centers(-1)
        assert np.allclose(trace.data['id'], [1, 2, 3])
        assert np.allclose(trace.data['centers'], [centers, centers+5, centers+10])
        trace.del_centers([-1, -2])
        assert np.allclose(trace.data['id'], [1])
        assert np.allclose(trace.data['centers'], [centers])

    def test_remove_duplicates(self):
        centers = np.array([1, 2, 3, 4])
        data = {'id': [0, 1, 2, 3],
                'centers': [centers, centers+1, centers+10, centers+15]}
        trace = Trace(data=data)
        trace.remove_duplicates(thresh=3)
        assert np.allclose(trace.data['id'], [0, 1, 2])
        assert np.allclose(trace.data['centers'].data, [centers+1, centers+10, centers+15])


class TestAllTraceFitter:
    @mock.patch('xwavecal.utils.trace_utils.SingleTraceFitter.update_initial_guess_to_run_through_pt')
    @mock.patch('xwavecal.utils.trace_utils.SingleTraceFitter.fit_trace')
    @mock.patch('xwavecal.utils.trace_utils.SingleTraceFitter.use_fit_as_initial_guess')
    def test_step_through_detector(self, use_last, fit, update):
        fit.return_value = np.arange(10)
        num_traces = 3
        fitter = SingleTraceFitter(extraargs={'initialize_fit_objects': False})
        trace = Trace(num_centers_per_trace=10)
        peak_xy_coordinates = list(zip(np.arange(num_traces), np.arange(num_traces)))
        trace = AllTraceFitter()._step_through_detector(trace, fitter, peak_xy_coordinates)
        for i in range(num_traces):
            assert np.allclose(trace.get_centers(i), np.arange(10))

    @mock.patch('xwavecal.utils.trace_utils.AllTraceFitter._step_through_detector')
    @mock.patch('xwavecal.utils.trace_utils.AllTraceFitter._identify_traces')
    @mock.patch('xwavecal.utils.trace_utils.SingleTraceFitter.__init__', return_value=None)
    def test_fit_traces(self, fitter, identify, step):
        step.return_value = Trace(data={'id': [0, 2, 1], 'centers': np.array([np.arange(3),
                                                                              np.arange(2, 5),
                                                                              np.arange(1, 4)])})
        trace = AllTraceFitter(min_peak_to_peak_spacing=0).fit_traces(None, None, None, None, None)
        for i in range(trace.num_traces_found()):
            assert np.allclose(trace.get_centers(i), np.arange(3)+i)
            assert np.isclose(trace.get_id(i), i)

    def test_getting_flux_down_detector(self):
        fake_trace_image = np.ones((5, 5))
        fake_trace_image[1] = 10
        fake_trace_image[3] = 10
        image_noise_estimate = np.sqrt(2)
        flux_signal = AllTraceFitter(xmin=3, xmax=5)._filtered_flux_down_detector(fake_trace_image,
                                                                                  image_noise_estimate)
        expected_slice = np.array([1, 10, 1, 10, 1])/np.sqrt(image_noise_estimate**2 + np.array([1, 10, 1, 10, 1]))
        assert np.allclose(flux_signal, expected_slice)

    @mock.patch('xwavecal.utils.trace_utils.AllTraceFitter._filtered_flux_down_detector')
    def test_identify_traces(self, flux):
        min_snr = 2
        peak_centers = np.array([15, 40, 60, 80])
        flux.return_value = array_with_peaks(x=np.arange(100).astype(float), centroids=peak_centers,
                                             amplitudes=[10, 6, min_snr+0.1, min_snr-0.1], stds=[3, 3, 3, 3])
        fitter = AllTraceFitter(min_snr=min_snr, min_peak_to_peak_spacing=5, xmin=0, xmax=2)
        peak_xy_coordinates = fitter._identify_traces(None, None)
        expected_peak_xy_coordinates = list(zip(np.ones_like(peak_centers[:-1]), peak_centers[:-1]))
        assert np.allclose(peak_xy_coordinates, expected_peak_xy_coordinates, atol=3)


class TestSingleTraceFitter:
    def test_class_attributes(self):
        fitter = SingleTraceFitter(extraargs={'initialize_fit_objects': False})
        assert fitter.second_order_coefficient_guess is None
        assert fitter.image_data is None
        assert fitter.filtered_image_data is None
        assert fitter.initial_guess_next_fit is None
        assert fitter.x is None
        assert fitter.x_norm is None
        assert fitter.design_matrix is None

    def test_default_class_attributes(self):
        fitter = SingleTraceFitter(extraargs={'initialize_fit_objects': False})
        assert fitter.poly_fit_order == 2
        assert fitter.coefficients == []

    def test_fit_initilization(self):
        """
        tests that SingleTraceFitter calls _initialize_fit_objects correctly upon
        instantiation of the class.
        """
        poly_fit_order = 2
        fitter = SingleTraceFitter(image_data=np.zeros((2, 2)),
                                   poly_fit_order=poly_fit_order,
                                   second_order_coefficient_guess=90)
        assert np.allclose(fitter.x_norm, np.array([-1, 1]))
        assert np.allclose(fitter.x, np.arange(2))
        assert fitter.filtered_image_data is not None
        assert fitter.design_matrix.shape == (poly_fit_order+1, 2)

    def test_generating_initial_guess(self):
        fitter = SingleTraceFitter(image_data=np.zeros((2, 2)),
                                   poly_fit_order=2,
                                   second_order_coefficient_guess=90)
        assert np.allclose(fitter.initial_guess_next_fit, np.array([0, 0, 90]))

    def test_generating_initial_guess_raises_error(self):
        with pytest.raises(Exception):
            SingleTraceFitter(image_data=np.zeros((2, 2)),
                              poly_fit_order=2,
                              second_order_coefficient_guess=None)

    def test_changing_initial_guesses(self):
        coefficients = [np.array([0, 0])]
        fitter = SingleTraceFitter(extraargs={'initialize_fit_objects': False,
                                              'coefficients': coefficients})
        fitter.use_fit_as_initial_guess(-1)
        assert np.allclose(fitter.initial_guess_next_fit, coefficients[-1])
        fitter.initial_guess_next_fit += 1
        assert not np.allclose(fitter.initial_guess_next_fit, coefficients[-1])

    @staticmethod
    def _generate_fitter():
        design_matrix = np.ones((2, 5))
        design_matrix[1] = np.linspace(-1, 1, 5)
        offset, linear_coefficient = 1, 0.5
        x = np.arange(design_matrix.shape[1])
        fitter = SingleTraceFitter(extraargs={'initialize_fit_objects': False,
                                              'design_matrix': design_matrix,
                                              'x': x})
        single_trace_coefficients = np.array([offset, linear_coefficient])
        fitter.initial_guess_next_fit = single_trace_coefficients
        a_line = np.linspace(offset - linear_coefficient, offset + linear_coefficient, 5)

        return fitter, a_line

    def test_centers_from_coefficients(self):
        fitter, a_line = self._generate_fitter()
        assert np.allclose(fitter._centers_from_coefficients(fitter.initial_guess_next_fit), a_line)

    def test_update_initial_guess_to_run_through_pt(self):
        fitter, a_line = self._generate_fitter()
        fitter.update_initial_guess_to_run_through_pt(xy=(1, 400))
        shifted_line = a_line - a_line[1] + 400
        assert np.allclose(fitter._centers_from_coefficients(fitter.initial_guess_next_fit), shifted_line)

    def test_flux_across_trace(self):
        x = np.arange(5)
        fake_data = np.zeros((9, len(x)))
        fake_data[3] += 1
        trace_centers = np.ones(len(x))*3
        filtered_fake_data = SingleTraceFitter._prefilter_image_data(fake_data)
        fitter = SingleTraceFitter(extraargs={'initialize_fit_objects': False,
                                              'x': x,
                                              'filtered_image_data': filtered_fake_data})
        expected_flux = len(x)
        flux = fitter._flux_across_trace(trace_centers)
        assert np.isclose(flux, expected_flux)

    def test_trace_merit_function_returns_negative_flux(self):
        x = np.arange(5)
        fake_data = np.zeros((9, len(x)))
        fake_data[3] += 1
        coefficients_for_line = np.array([3, 0])
        design_matrix = np.array([np.ones(len(x)),
                                  np.linspace(-1, 1, len(x))])
        filtered_fake_data = SingleTraceFitter._prefilter_image_data(fake_data)
        fitter = SingleTraceFitter(extraargs={'initialize_fit_objects': False,
                                              'x': x,
                                              'design_matrix': design_matrix,
                                              'filtered_image_data': filtered_fake_data})
        negative_flux = (-1)*len(x)
        merit = fitter._trace_merit_function(single_trace_coefficients=coefficients_for_line,
                                             cls=fitter)
        assert np.isclose(merit, negative_flux)

    def test_normalizing_coordinates(self):
        x = np.arange(5)
        x_norm = np.linspace(-1, 1, len(x))
        fake_data = np.zeros((9, len(x)))
        fitter = SingleTraceFitter(image_data=fake_data, extraargs={'initialize_fit_objects': False})
        fitter._normalize_domain_coordinates()
        assert np.allclose(fitter.x, x)
        assert np.allclose(fitter.x_norm, x_norm)

    def test_generating_legendre_design_matrix(self):
        x = np.arange(5)
        xnorm = np.linspace(-1, 1, len(x))
        design_matrix = np.array([np.ones(len(x)),
                                  xnorm])
        assert np.allclose(design_matrix,
                           SingleTraceFitter._generate_design_matrix(xnorm, poly_fit_order=1))

    @mock.patch('xwavecal.utils.trace_utils.SingleTraceFitter._centers_from_coefficients')
    @mock.patch('xwavecal.utils.trace_utils.optimize.minimize')
    def test_fit_trace(self, minimize, centers):
        coefficients = np.arange(4)
        centers.return_value = None
        minimize.return_value = OptimizeResult(x=coefficients)
        fitter = SingleTraceFitter(extraargs={'initialize_fit_objects': False,
                                              'coefficients': [],
                                              'initial_guess_next_fit': coefficients})

        assert fitter.fit_trace() is None
        assert fitter.coefficients == [coefficients]


class TestTraceMaker:
    def test_properties(self):
        assert TraceMaker(FakeContext()).calibration_type is 'TRACE'

    @pytest.mark.integration
    def test_trace_fit_does_not_crash_on_blank_frame(self):
        order_of_poly_fit = 4
        image = FakeTraceImage(nx=100, ny=100)
        image.set_header_val('RDNOISE', 11)
        noisify_image(image)
        fake_context = FakeContext()
        fake_context.db_address = ''
        trace_fitter = TraceMaker(fake_context)
        trace_fitter.order_of_poly_fit = order_of_poly_fit
        trace_fitter.xmin, trace_fitter.xmax = 50, 60
        trace_fitter.do_stage(image)
        assert True

    @mock.patch('xwavecal.utils.trace_utils.AllTraceFitter.fit_traces')
    def test_trace_maker(self, fit_traces):
        trace_table_name = 'test'
        data = {'id': [1], 'centers': [np.arange(3)]}
        fit_traces.return_value = Trace(data=data)
        expected_trace = Trace(data=data)
        trace_maker = TraceMaker(FakeContext())
        trace_maker.xmin = 5
        trace_maker.xmax = 10
        trace_maker.trace_table_name = trace_table_name
        loaded_trace = trace_maker.do_stage(image=FakeImage())[1]
        assert np.allclose(loaded_trace.get_centers(0), expected_trace.get_centers(0))
        assert np.allclose(loaded_trace.get_id(0), expected_trace.get_id(0))

    @pytest.mark.integration
    def test_accuracy_of_trace_fitting(self):
        """
        test type: Mock Integration Test with metrics for how well trace fitting is doing.
        info: This tests trace making via a blind fit.
        WARNING: Because trace fitting is defined with polynomials which are normalized from -1 to 1, if one squeezes
        the x axis of the image, then the traces bend more drastically. Thus it is recommended you do not change the
        size of the FakeTraceImage.
        """
        read_noise = 11.0
        poly_fit_order = 4

        image = FakeTraceImage()
        image.fiber0_lit, image.fiber1_lit, image.fiber2_lit = False, True, True
        image.set_header_val('RDNOISE', read_noise)

        image, trace_centers, second_order_coefficient_guess = fill_image_with_traces(image,
                                                                                      poly_order_of_traces=poly_fit_order)
        noisify_image(image)
        fake_context = FakeContext()

        trace_maker = TraceMaker(fake_context)
        trace_maker.xmin = image.data.shape[1]//2 - 20
        trace_maker.xmax = image.data.shape[1]//2 + 20
        trace_maker.order_of_poly_fit = poly_fit_order
        trace_maker.second_order_coefficient_guess = second_order_coefficient_guess
        trace = trace_maker.do_stage(image)[1]
        assert trace.data['centers'].shape[0] == trace_centers.shape[0]
        difference = trace.data['centers'] - trace_centers
        logger.debug('median absolute deviation in unit-test trace fitting is {0} pixels'
                     .format(np.median(np.abs(difference - np.median(difference)))))
        logger.debug('standard deviation in unit-test trace fitting is {0} pixels'
                     .format(np.std(difference)))
        logger.debug('worst error (max of true minus found) in unit-test trace fitting is {0} pixels'
                     .format(np.max(np.abs(difference))))
        logger.debug('median error (median of true minus found) in unit-test trace fitting is {0} pixels'
                     .format(np.abs(np.median(difference))))
        assert np.median(np.abs(difference - np.median(difference))) < 1/100
        assert np.abs(np.median(difference)) < 1/100


class TestLoadTrace:
    def test_properties(self):
        assert LoadTrace(FakeContext()).calibration_type is 'TRACE'

    @mock.patch('xwavecal.traces.LoadTrace.get_calibration_filename', return_value=None)
    def test_load_trace_flags_images_without_calibration(self, mock_get_cal):
        fake_context = FakeContext()
        setattr(fake_context, 'db_address', None)
        trace_loader = LoadTrace(fake_context)
        image = trace_loader.do_stage(image=FakeImage())
        assert getattr(image, 'trace') is None

    @mock.patch('xwavecal.traces.LoadTrace.get_calibration_filename', return_value='/path/to/master_trace.fits')
    @mock.patch('os.path.exists', return_value=True)
    @mock.patch('astropy.io.fits.open', return_value=None)
    @mock.patch('xwavecal.utils.trace_utils.Trace.load')
    def test_load_trace(self, mock_load, mock_open, mock_os, mock_get_cal):
        data = {'id': [1], 'centers': [np.arange(3)]}
        expected_trace = Trace(data=data)
        mock_load.return_value = Trace(data=data)
        fake_context = FakeContext()
        setattr(fake_context, 'db_address', None)
        trace_loader = LoadTrace(fake_context)
        image = trace_loader.do_stage(image=FakeImage())
        assert np.allclose(image.trace.get_centers(0), expected_trace.get_centers(0))
        assert np.allclose(image.trace.get_id(0), expected_trace.get_id(0))
