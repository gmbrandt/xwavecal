import mock
import numpy as np
import copy
from astropy.table import Table

from echelle.blaze import ApplyBlaze, BackgroundSubtractSpectrum
import echelle.settings as settings
from echelle.images import DataProduct
from echelle.tests.utils import FakeContext, FakeImage


@mock.patch('echelle.blaze.DataProduct.load')
class TestApplyBlaze:
    def test_apply_calibration(self, fake_load):
        stage = ApplyBlaze(FakeContext())
        fake_blaze = Table({'id': [1, 2], 'flux': [np.arange(1, 11), np.arange(2, 12)]})
        fake_load.return_value = DataProduct(data=fake_blaze)
        image = FakeImage()
        image.fiber0_wavecal, image.fiber1_wavecal, image.fiber2_wavecal = 0, 1, 1
        image.fiber0_lit, image.fiber1_lit, image.fiber2_lit = 0, 1, 1
        image.data_tables = {settings.BOX_SPECTRUM_NAME: copy.deepcopy(fake_blaze)}
        image = stage.apply_master_calibration(image, 'fake/path')
        assert np.allclose(image.data_tables[settings.BLAZE_CORRECTED_BOX_SPECTRUM_NAME]['flux'].data, 1)
        assert not np.allclose(image.data_tables[settings.BLAZE_CORRECTED_BOX_SPECTRUM_NAME]['flux'].data,
                               image.data_tables[settings.BOX_SPECTRUM_NAME]['flux'].data)


class TestBackgroundSubtractSpectrum:
    def test_do_stage(self):
        stage = BackgroundSubtractSpectrum(FakeContext())
        image = FakeImage()
        image.data_tables = {'SPECBOX': Table({'fiber': [1, 1], 'flux': np.ones((2, 10))})}
        image = stage.do_stage(image)
        assert np.allclose(image.data_tables['SPECBOX']['flux'].data, np.zeros((2, 10)))
