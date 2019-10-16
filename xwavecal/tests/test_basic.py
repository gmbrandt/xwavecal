import numpy as np
from astropy.table import Table

from xwavecal.tests.utils import FakeContext, FakeImage
from xwavecal.utils.basic_utils import median_subtract_channels_y
from xwavecal import basic


class TestBasic:
    CONTEXT = FakeContext()

    def test_gain_normalizer(self):
        image = FakeImage()
        image.data = np.ones((10, 10))
        image.set_header_val('gain', 2)
        image = basic.GainNormalizer(self.CONTEXT).do_stage(image)
        assert image.get_header_val('gain') == 1
        assert np.allclose(image.data, 2)

    def test_overscan_subtractor(self):
        image = FakeImage()
        image.data = np.ones((10, 12)).astype(float)
        image.data[:, 10:] = .5
        image.set_header_val('data_section', '[1:10,1:10]')
        image.set_header_val('overscan_section', '[1:10, 11:12]')
        image = basic.OverscanSubtractor(self.CONTEXT).do_stage(image)
        assert np.allclose(image.data[:10, :10], 0.5)

    def test_overscan_trimmer(self):
        image = FakeImage()
        image.data = np.ones((10, 12)).astype(float)
        image.data[:, 10:] = .5
        image.set_header_val('data_section', '[1:10,1:10]')
        image = basic.Trimmer(self.CONTEXT).do_stage(image)
        assert np.allclose(image.data.shape, (10, 10))


class TestBackgroundSubtract1dSpectrum:
    def test_do_stage(self):
        stage = basic.BackgroundSubtractSpectrum(FakeContext())
        image = FakeImage()
        image.data_tables = {'SPECBOX': Table({'fiber': [1, 1], 'flux': np.ones((2, 10)),
                                               'stderr': np.ones((2, 10))})}
        image = stage.do_stage(image)
        assert np.allclose(image.data_tables['SPECBOX']['flux'].data, np.zeros((2, 10)))


class TestMedianSubtractReadoutsAlongY:
    def test_do_stage(self):
        image = FakeImage()
        image.data = np.ones((4, 4))
        image.header['num_rd_channels'] = 1
        image = basic.MedianSubtractReadoutsAlongY(None).do_stage(image)
        assert np.allclose(image.data, 0)


class TestUtils:
    def test_median_subtract_channels(self):
        a = (np.arange(3) * np.ones((3, 3))).T
        assert np.allclose(median_subtract_channels_y(a, 3), 0)
        a = (np.array([1, 1, 1, 2, 2, 2]) * np.ones((6, 6))).T
        assert np.allclose(median_subtract_channels_y(a, 2), 0)
