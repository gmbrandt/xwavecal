import numpy as np

from xwavecal.tests.utils import FakeImage, FakeContext
from xwavecal.variance import CalcInverseVariance


class TestCalcInverseVariance:
    def test_calc_zero_rdnoise(self):
        image = FakeImage()
        image.header = {'read_noise': 0}
        image.data = np.random.random((10, 10)) * 100
        image = CalcInverseVariance(FakeContext()).do_stage(image)
        assert np.allclose(image.ivar, image.data ** (-1))

    def test_calc_zero_signal(self):
        image = FakeImage()
        image.header = {'read_noise': 10}
        image.data = np.zeros((10, 10))
        image = CalcInverseVariance(FakeContext()).do_stage(image)
        assert np.allclose(image.ivar, image.get_header_val('read_noise') ** (-2))

    def test_calc(self):
        image = FakeImage()
        image.header = {'read_noise': 10}
        image.data = np.random.random((10, 10)) * 100 - 50
        image = CalcInverseVariance(FakeContext()).do_stage(image)
        assert np.min(image.ivar ** -1) >= image.get_header_val('read_noise') ** (2)

