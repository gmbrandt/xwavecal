import mock
import numpy as np
from astropy.table import Table

from echelle.blaze import ApplyBlaze, BackgroundSubtractSpectrum, BlazeMaker
from echelle.tests.utils import FakeContext, FakeImage
from echelle.utils.blaze_utils import normalize_orders


class TestApplyBlaze:
    CONTEXT = FakeContext()
    @mock.patch('echelle.blaze.ApplyBlaze.get_calibration_filename', return_value='fake')
    @mock.patch('echelle.images.Image.load')
    def test_apply_calibration(self, mock_load, mock_get_cal):
        blaze = FakeImage()
        blaze.data = np.random.randint(1, 9, size=(10, 10))
        mock_load.return_value = blaze
        image = ApplyBlaze(self.CONTEXT).do_stage(blaze)
        assert np.allclose(image.data, 1)

    @mock.patch('echelle.blaze.ApplyBlaze.get_calibration_filename', return_value='fake')
    @mock.patch('echelle.images.Image.load')
    def test_apply_calibration_aborts_on_incorrect_shape(self, mock_load, mock_get_cal):
        blaze = FakeImage()
        blaze.data = np.random.randint(1, 9, size=(10, 11))
        mock_load.return_value = blaze
        image = ApplyBlaze(self.CONTEXT).do_stage(np.ones((10, 10)).astype(float))
        assert np.allclose(image.data, 1)


class TestBlazeMaker:
    CONTEXT = FakeContext()
    @mock.patch('echelle.blaze.normalize_orders')
    def test_do_stage(self, mock_normalize):
        mock_normalize.return_value = np.zeros((3, 3))
        image = FakeImage()
        image.data = np.random.randint(1, 9, size=(3, 3))
        image, blaze = BlazeMaker(self.CONTEXT).do_stage(image)
        assert not np.allclose(image.data, blaze.data)
        assert image.get_header_val('type') == 'lampflat'
        assert blaze.get_header_val('type') == BlazeMaker(self.CONTEXT).calibration_type.lower()
        assert np.allclose(blaze.data, 0)


class TestBackgroundSubtractSpectrum:
    def test_do_stage(self):
        stage = BackgroundSubtractSpectrum(FakeContext())
        image = FakeImage()
        image.data_tables = {'SPECBOX': Table({'fiber': [1, 1], 'flux': np.ones((2, 10))})}
        image = stage.do_stage(image)
        assert np.allclose(image.data_tables['SPECBOX']['flux'].data, np.zeros((2, 10)))


def test_normalize_orders():
    data = np.ones((10, 100)).astype(float)
    data[4:7] *= 3
    data[1:3] *= 0.5
    trace = type('Trace', (), {'data': {'centers': [5 * np.ones((100))]}})
    data = normalize_orders(data, trace, minval=0.6, half_window=4, n=10)
    assert np.allclose(data[4:7], 1)
    assert np.allclose(data[1:3], 0.6/3)
    # check that we do not divide values which are not near the trace and which are not < minval
    assert np.allclose([data[-1], data[0]], 1)