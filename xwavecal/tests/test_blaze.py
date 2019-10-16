import mock
import numpy as np
from astropy.table import Table
from copy import deepcopy

from xwavecal.blaze import ApplyBlaze, BlazeMaker
from xwavecal.tests.utils import FakeContext, FakeImage
from xwavecal.utils.blaze_utils import normalize_orders


class TestApplyBlaze:
    CONTEXT = FakeContext()

    @mock.patch('os.path.exists', return_value=True)
    @mock.patch('xwavecal.blaze.ApplyBlaze.get_calibration_filename', return_value='fake')
    @mock.patch('xwavecal.images.Image.load')
    def test_apply_calibration(self, mock_load, mock_get_cal, mock_os):
        blaze = FakeImage()
        self.CONTEXT.min_blaze_sn = np.inf
        blaze.data = np.random.randint(1, 9, size=(10, 10))
        image = deepcopy(blaze)
        image.ivar = 4 * np.ones((10, 10))  # stderr = 1/2
        mock_load.return_value = blaze
        image = ApplyBlaze(self.CONTEXT).do_stage(image)
        assert np.allclose(image.data, 1)
        assert np.allclose(np.sqrt(image.ivar)**(-1), (1/2) / blaze.data)

    @mock.patch('xwavecal.blaze.ApplyBlaze.get_calibration_filename', return_value='fake')
    @mock.patch('xwavecal.images.Image.load')
    def test_apply_calibration_aborts_on_incorrect_shape(self, mock_load, mock_get_cal):
        blaze = FakeImage()
        blaze.data = np.random.randint(1, 9, size=(10, 11))
        image = FakeImage()
        image.data = np.ones((10, 10)).astype(float)
        mock_load.return_value = blaze
        image = ApplyBlaze(self.CONTEXT).do_stage(image)
        assert np.allclose(image.data, 1)


class TestBlazeMaker:
    CONTEXT = FakeContext()
    @mock.patch('xwavecal.blaze.normalize_orders')
    def test_do_stage(self, mock_normalize):
        mock_normalize.return_value = 2 * np.ones((3, 3), dtype=float)
        image = FakeImage()
        image.data = np.random.randint(1, 9, size=(3, 3))
        image, blaze = BlazeMaker(self.CONTEXT).do_stage(image)
        assert np.allclose(image.data / 2, blaze.data)
        assert image.get_header_val('type') == 'lampflat'
        assert blaze.get_header_val('type') == BlazeMaker(self.CONTEXT).calibration_type.lower()


def test_normalize_orders():
    data = np.ones((10, 100)).astype(float)
    data[4:7] *= 3
    data[1:3] *= 0.5
    trace = type('Trace', (), {'data': {'centers': [5 * np.ones((100))]}})
    normalization_factor = normalize_orders(data, trace, half_window=4, n=10)
    data /= normalization_factor
    assert np.allclose(data[4:7], 1)
    assert np.allclose(data[1:3], 0.5/3)
    # check that we set values which are not near the trace and which are not < minval to 0
    assert np.allclose([data[-1], data[0]], 0)
