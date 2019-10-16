import mock
import numpy as np

import tempfile
import os
from astropy.table import Table
from astropy.io.fits import Header


from xwavecal.images import Image, DataProduct
from xwavecal.utils.fits_utils import Translator
from xwavecal.tests.utils import FakeContext, FakeImage


def test_image_class_loads():
    image = Image(header={'fiber_state': 'tung&tung&none'})
    assert image.trace is None
    assert image.rectified_2d_spectrum is None
    assert image.wavelength_solution == {}


def test_get_num_lit_fibers():
    image = Image(header={'fiber_state': 'tung&tung&none'})
    assert image.num_lit_fibers() == 2
    image = Image(header={'fiber_state': 'none&tung&none'})
    assert image.num_lit_fibers() == 1


def test_get_num_lit_wavecal_fibers():
    image = Image(header={'fiber_state': 'thar&tung&none'})
    assert image.num_wavecal_fibers() == 1
    assert image.fiber0_wavecal == 1
    assert image.fiber1_wavecal == 0
    image = Image(header={'fiber_state': 'none&thar&thar'})
    assert image.num_wavecal_fibers() == 2
    assert image.fiber1_wavecal == 1
    assert image.fiber2_wavecal == 1


class TestDataProduct:
    def test_load_and_write(self):
        name = 'trace'
        image = DataProduct(data={'id': [1], 'centers': [np.arange(3)]}, data_name='trace')
        with tempfile.TemporaryDirectory() as tmp_directory_name:
            path = os.path.join(tmp_directory_name, 'test_trace_table.fits')
            image.filepath = path
            image.header = {'bla': 1}
            image.write(False)
            loaded_image = DataProduct.load(path=path, extension_name=name)
            assert np.allclose(loaded_image.data['centers'][0], image.data['centers'][0])
            assert np.allclose(loaded_image.data['id'][0], image.data['id'][0])
            assert np.isclose(loaded_image.get_header_val('bla'), 1)

    def test_write_gets_correct_filename(self):
        name = 'trace'
        image = DataProduct(data=Table({'id': [1], 'centers': [np.arange(3)]}), data_name=name)
        with tempfile.TemporaryDirectory() as tmp_directory_name:
            for fpack, extension in zip([True, False], ['.fz', 'its']):
                path = os.path.join(tmp_directory_name, 'test_trace_table.fits')
                image.filepath = path
                image.header = {'bla': 1}
                image._update_filepath(fpack)
                assert image.filepath[-3:] == extension

    def test_translator(self):
        header_keys = {'type': 'OBSTYPE',
                       'gain': 'GAIN',
                       'read_noise': 'RDNOISE',
                       'multi': ('obstype', 'gain')}
        type_keys = {'LAMPFLAT': 'lampflat',
                     'DOUBLE': 'wavecal'}
        translator = Translator(header_keys, type_keys)
        header = Header({'OBSTYPE': 'DOUBLE', 'GAIN': 1, 'RDNOISE': 10})
        image = DataProduct(header=header, translator=translator)
        assert image.get_header_val('type') == 'wavecal'
        assert image.get_header_val('gain') == 1
        assert image.get_header_val('read_noise') == 10
        assert image.get_header_val('multi') == ('DOUBLE', 1)
        image.set_header_val('gain', 5)
        assert image.get_header_val('gain') == 5


class TestImage:
    def test_load_and_write(self):
        image = Image(data=np.ones((10, 10)), header={'bla': 1, 'fiber_state': 'none&none&none'})
        with tempfile.TemporaryDirectory() as tmp_directory_name:
            path = os.path.join(tmp_directory_name, 'test_trace_table.fits')
            image.filepath = path
            image.write(fpack=False)
            image = Image.load(path=path, extension_name=0)
            assert np.allclose(image.data, 1)
            assert image.header['bla'] == 1

    def test_load_and_write_with_name(self):
        name = 'a_name'
        image = Image(data=np.ones((10, 10)), data_name='a_name', header={'bla': 1, 'fiber_state': 'none&none&none'})
        with tempfile.TemporaryDirectory() as tmp_directory_name:
            path = os.path.join(tmp_directory_name, 'test_trace_table.fits')
            image.filepath = path
            image.write(fpack=False)
            image = Image.load(path=path, extension_name=name)
            assert np.allclose(image.data, 1)
            assert image.header['bla'] == 1
