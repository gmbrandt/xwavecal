import mock
import numpy as np

import tempfile
import os
from astropy.table import Table
from astropy.io.fits import Header


from echelle.images import Image, DataProduct
from echelle.tests.utils import FakeContext, FakeImage


@mock.patch('banzai.images.Image._init_instrument_info')
def test_image_class_loads(mock_instrument):
    mock_instrument.return_value = None, None, None
    image = Image(runtime_context=FakeContext(),
                      header={'OBJECTS': 'tung&tung&none'})
    assert image.trace is None
    assert image.rectified_2d_spectrum is None
    assert image.wavelength_solution == {}


@mock.patch('banzai.images.Image._init_instrument_info')
def test_get_num_lit_fibers(mock_instrument):
    mock_instrument.return_value = None, None, None
    image = Image(runtime_context=FakeContext(),
                      header={'OBJECTS': 'tung&tung&none'})
    assert image.num_lit_fibers() == 2
    image = Image(runtime_context=FakeContext(),
                      header={'OBJECTS': 'none&tung&none'})
    assert image.num_lit_fibers() == 1


@mock.patch('banzai.images.Image._init_instrument_info')
def test_get_num_lit_wavecal_fibers(mock_instrument):
    mock_instrument.return_value = None, None, None
    image = Image(runtime_context=FakeContext(),
                      header={'OBJECTS': 'thar&tung&none'})
    assert image.num_wavecal_fibers() == 1
    assert image.fiber0_wavecal == 1
    assert image.fiber1_wavecal == 0
    image = Image(runtime_context=FakeContext(),
                      header={'OBJECTS': 'none&thar&thar'})
    assert image.num_wavecal_fibers() == 2
    assert image.fiber1_wavecal == 1
    assert image.fiber2_wavecal == 1


class TestDataProduct:
    def test_load_and_write(self):
        name = 'trace'
        image = DataProduct(data={'id': [1], 'centers': [np.arange(3)]}, data_name=name)
        with tempfile.TemporaryDirectory() as tmp_directory_name:
            path = os.path.join(tmp_directory_name, 'test_trace_table.fits')
            image.filepath = path
            image.header = {'bla': 1}
            image.write(False)
            loaded_image = DataProduct.load(path=path, extension_name=name)
            assert np.allclose(loaded_image.data['centers'][0], image.data['centers'][0])
            assert np.allclose(loaded_image.data['id'][0], image.data['id'][0])
            assert np.isclose(loaded_image.header['bla'], 1)

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

    def test_translate_header(self):
        header_keys = {'type': 'OBSTYPE',
                       'gain': 'GAIN',
                       'read_noise': 'RDNOISE'}
        type_translator = {'LAMPFLAT': 'lampflat',
                           'DOUBLE': 'wavecal'}
        header = Header({'OBSTYPE': 'DOUBLE', 'GAIN': 1, 'RDNOISE': 10})
        new = DataProduct.translate_header(header, header_keys, type_translator)
        assert new['type'] == 'wavecal'

