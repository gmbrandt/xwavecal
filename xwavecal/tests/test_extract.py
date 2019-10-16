import numpy as np
import pytest

from xwavecal.tests.test_traces import FakeTraceImage
from xwavecal.tests.utils import fill_image_with_traces, FakeContext
from xwavecal.utils.trace_utils import Trace
from xwavecal.utils import extract_utils
from xwavecal.extract import BoxExtract, RectifyTwodSpectrum, IVarExtract, safe_pow


class TestRectify:
    CONTEXT = FakeContext()
    def test_rectify_orders(self):
        image = FakeTraceImage()
        hw = 10
        peak_intensity = 1E4
        image, trace_centers, second_order_coefficient_guess = fill_image_with_traces(image, poly_order_of_traces=4,
                                                                                      max_num_traces=3,
                                                                                      fiber_intensity=peak_intensity)
        trace = Trace(data={'id': np.arange(trace_centers.shape[0]), 'centers': trace_centers})
        rectified_orders, zeroed_image_data = extract_utils.rectify_orders(image.data, trace,
                                                                           half_window=hw,
                                                                           debug=True)
        assert np.allclose(zeroed_image_data, 0)
        assert not np.allclose(image.data, 0)
        for key, rect_spect in rectified_orders.items():
            assert np.isclose(np.median(rect_spect['val'][hw]), peak_intensity, rtol=0.02)

    def test_rectify_orders_assigns_coordinates_at_edge(self):
        trace = Trace(data={'id': [1], 'centers': [np.zeros(3)]})
        image_data = np.zeros((3, 3))
        hw = 1
        rectified_orders, zeroed_image_data = extract_utils.rectify_orders(image_data, trace,
                                                                           half_window=hw,
                                                                           debug=True)
        assert np.allclose(rectified_orders[1]['y'], np.array([[0, 0, 0],
                                                               [0, 0, 0],
                                                               [1, 1, 1]]))
        assert np.allclose(rectified_orders[1]['x'], np.array([[0, 1, 2],
                                                               [0, 1, 2],
                                                               [0, 1, 2]]))

    def test_rectify_orders_assigns_coordinates(self):
        trace = Trace(data={'id': [1], 'centers': [np.ones(3)]})
        image = FakeTraceImage()
        image.data = np.zeros((3, 3))
        hw = 1
        rectified_orders, zeroed_image_data = extract_utils.rectify_orders(image.data, trace,
                                                                           half_window=hw,
                                                                           debug=True)
        assert np.allclose(rectified_orders[1]['y'], np.array([[0, 0, 0],
                                                               [1, 1, 1],
                                                               [2, 2, 2]]))
        assert np.allclose(rectified_orders[1]['x'], np.array([[0, 1, 2],
                                                               [0, 1, 2],
                                                               [0, 1, 2]]))

    def test_rectify_curved_order_maps_all_values(self):
        image = FakeTraceImage()
        image, trace_centers, second_order_coefficient_guess = fill_image_with_traces(image, poly_order_of_traces=4,
                                                                                      max_num_traces=1)
        x_coordinates, y_coordinates = np.meshgrid(np.arange(image.data.shape[1]), np.arange(image.data.shape[0]))
        image_coordinates = {'x': x_coordinates, 'y': y_coordinates}
        single_order_centers = trace_centers[0]
        rectified_order, zeroed_image_data = extract_utils.rectify_order(image.data, image_coordinates,
                                                                         single_order_centers, half_window=10,
                                                                         nullify_mapped_values=True)
        assert np.allclose(zeroed_image_data, 0)

    def test_rectify_flat_order(self):
        hw = 10
        image = FakeTraceImage()
        image, trace_centers, second_order_coefficient_guess = fill_image_with_traces(image, poly_order_of_traces=1,
                                                                                      max_num_traces=1)
        x_coordinates, y_coordinates = np.meshgrid(np.arange(image.data.shape[1]), np.arange(image.data.shape[0]))
        image_coordinates = {'x': x_coordinates, 'y': y_coordinates}
        single_order_centers = trace_centers[0]
        rectified_order, image_data = extract_utils.rectify_order(image.data, image_coordinates,
                                                                  single_order_centers, half_window=hw,
                                                                  nullify_mapped_values=False)
        trace_y_value = int(trace_centers[0][0])
        assert np.allclose(rectified_order['val'], image_data[trace_y_value - hw: trace_y_value + hw + 1, :])

    def test_rectification_does_not_change_box_extract(self):
        image = FakeTraceImage()
        image, trace_centers, second_order_coefficient_guess = fill_image_with_traces(image, poly_order_of_traces=4,
                                                                                      max_num_traces=1)
        x_coordinates, y_coordinates = np.meshgrid(np.arange(image.data.shape[1]), np.arange(image.data.shape[0]))
        image_coordinates = {'x': x_coordinates, 'y': y_coordinates}
        single_order_centers = trace_centers[0]
        rectified_order, image_data = extract_utils.rectify_order(image.data, image_coordinates,
                                                                  single_order_centers, half_window=10,
                                                                  nullify_mapped_values=False)
        rectified_spectrum = BoxExtract(self.CONTEXT).extract_order(rectified_order['val'])
        spectrum = BoxExtract(self.CONTEXT).extract_order(image_data)
        assert np.allclose(spectrum / np.median(spectrum), 1)
        assert np.allclose(rectified_spectrum, spectrum)

    def test_empty_spectrum_on_missing_trace(self):
        image = FakeTraceImage()
        image.trace = None
        image = RectifyTwodSpectrum(self.CONTEXT).do_stage(image=image)
        assert image.rectified_2d_spectrum == {}


class TestBoxExtract:
    CONTEXT = FakeContext()
    def test_box_extract_accuracy(self):
        image = FakeTraceImage()
        image, trace_centers, second_order_coefficient_guess = fill_image_with_traces(image, poly_order_of_traces=4,
                                                                                      max_num_traces=1)
        spectrum = BoxExtract(self.CONTEXT).extract_order(image.data)
        assert np.allclose(spectrum / np.median(spectrum), 1)

    def test_box_extract_trims_rectified_2d_spectrum(self):
        max_extract_window = 10
        for half_window in [2, 6, 10, 15]:
            fake_spectrum = np.zeros((2 * max_extract_window + 1, 5))
            fake_spectrum[max_extract_window] = 1
            extractor = BoxExtract(self.CONTEXT)
            extractor.extraction_half_window = half_window
            extractor.max_extraction_half_window = max_extract_window
            trimmed_spectrum = extractor._trim_rectified_2d_spectrum(rectified_2d_spectrum={1: {'val': fake_spectrum}})
            hw = min(half_window, max_extract_window)
            assert np.isclose(trimmed_spectrum[1]['val'].shape[0], 2*hw + 1)
            assert np.allclose(trimmed_spectrum[1]['val'][hw], 1)

    @pytest.mark.integration
    def test_box_extract(self):
        image = FakeTraceImage()
        image, trace_centers, second_order_coefficient_guess = fill_image_with_traces(image, poly_order_of_traces=4,
                                                                                      max_num_traces=2,
                                                                                      fiber_intensity=1E4)
        image.trace = Trace(data={'id': np.arange(trace_centers.shape[0]), 'centers': trace_centers})
        image = RectifyTwodSpectrum(self.CONTEXT).do_stage(image)
        image = BoxExtract(self.CONTEXT).do_stage(image)
        for spectrum in image.data_tables[self.CONTEXT.box_spectrum_name]['flux']:
            assert np.median(spectrum) > 1E4


class TestIVarExtract:
    CONTEXT = FakeContext()

    def test_weights_are_normalized(self):
        w = IVarExtract(self.CONTEXT)._weights(None, np.random.random((10, 20)))
        assert np.allclose(w.shape, (10, 20))
        assert np.allclose(np.sum(w, axis=0), 1)

    def test_weights_handle_div_by_zero(self):
        w = IVarExtract(self.CONTEXT)._weights(None, np.zeros((10, 20)))
        assert np.allclose(w.shape, (10, 20))
        assert np.allclose(w, 0)


def test_safe_pow():
    power = np.random.randint(1, 4)
    assert safe_pow(None, power) is None
    assert np.isclose(safe_pow(5, power), np.power(5, power))
