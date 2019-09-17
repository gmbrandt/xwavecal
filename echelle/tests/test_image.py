import mock
import numpy as np
from banzai_nres.images import NRESImage
from banzai.tests.utils import FakeContext, FakeImage

import tempfile
import os
from astropy.table import Table


from banzai_nres.images import NRESImage, ImageBase
from banzai.tests.utils import FakeContext, FakeImage
import banzai_nres.settings as nres_settings
from banzai import settings


@mock.patch('banzai.images.Image._init_instrument_info')
def test_image_class_loads(mock_instrument):
    mock_instrument.return_value = None, None, None
    image = NRESImage(runtime_context=FakeContext(),
                      header={'OBJECTS': 'tung&tung&none'})
    assert image.trace is None
    assert image.rectified_2d_spectrum is None
    assert image.wavelength_solution == {}


@mock.patch('banzai.images.Image._init_instrument_info')
def test_get_num_lit_fibers(mock_instrument):
    mock_instrument.return_value = None, None, None
    image = NRESImage(runtime_context=FakeContext(),
                      header={'OBJECTS': 'tung&tung&none'})
    assert image.num_lit_fibers() == 2
    image = NRESImage(runtime_context=FakeContext(),
                      header={'OBJECTS': 'none&tung&none'})
    assert image.num_lit_fibers() == 1


@mock.patch('banzai.images.Image._init_instrument_info')
def test_get_num_lit_wavecal_fibers(mock_instrument):
    mock_instrument.return_value = None, None, None
    image = NRESImage(runtime_context=FakeContext(),
                      header={'OBJECTS': 'thar&tung&none'})
    assert image.num_wavecal_fibers() == 1
    assert image.fiber0_wavecal == 1
    assert image.fiber1_wavecal == 0
    image = NRESImage(runtime_context=FakeContext(),
                      header={'OBJECTS': 'none&thar&thar'})
    assert image.num_wavecal_fibers() == 2
    assert image.fiber1_wavecal == 1
    assert image.fiber2_wavecal == 1


class TestImageBase:
    def test_load_and_write(self):
        name = 'trace'
        image = ImageBase(data={'id': [1], 'centers': [np.arange(3)]}, table_name=name)
        runtime_context = FakeContext()
        with tempfile.TemporaryDirectory() as tmp_directory_name:
            runtime_context.fpack = False
            path = os.path.join(tmp_directory_name, 'test_trace_table.fits')
            image.filepath = path
            image.header = {'bla': 1}
            image.write(runtime_context, update_db=False)
            loaded_image = ImageBase.load(path=path, extension_name=name)
            assert np.allclose(loaded_image.data['centers'][0], image.data['centers'][0])
            assert np.allclose(loaded_image.data['id'][0], image.data['id'][0])
            assert np.isclose(loaded_image.header['bla'], 1)

    def test_write_gets_correct_filename(self):
        name = 'trace'
        image = ImageBase(data=Table({'id': [1], 'centers': [np.arange(3)]}), table_name=name)
        runtime_context = FakeContext()
        with tempfile.TemporaryDirectory() as tmp_directory_name:
            for fpack, extension in zip([True, False], ['.fz', 'its']):
                runtime_context.fpack = fpack
                path = os.path.join(tmp_directory_name, 'test_trace_table.fits')
                image.filepath = path
                image.header = {'bla': 1}
                image._update_filepath(runtime_context)
                assert image.filepath[-3:] == extension

    def test_instantiates_db_context_from_image(self):
        image = FakeImage()
        possible_attributes = ['dateobs', 'datecreated', 'instrument', 'is_master', 'is_bad']
        counter = np.arange(len(possible_attributes))
        for attribute, i in zip(possible_attributes, counter):
            setattr(image, attribute, str(i))
        imagebase = ImageBase(data={'id': [1], 'centers': [np.arange(3)]}, image=image, obstype='TRACE')
        for attribute in possible_attributes:
            assert getattr(imagebase, attribute) == getattr(image, attribute)
        assert imagebase.attributes == settings.CALIBRATION_SET_CRITERIA.get('TRACE', {})
