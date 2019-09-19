import tempfile
import os
from configparser import ConfigParser
import mock

from echelle.main import select_data, run, reduce_data, organize_config
from echelle.tests.utils import FakeImage


def test_select_data():
    header_keys = {'type': 'OBSTYPE',
                   'gain': 'GAIN',
                   'read_noise': 'RDNOISE'}
    type_translator = {'LAMPFLAT': 'lampflat',
                       'DOUBLE': 'wavecal'}
    with tempfile.TemporaryDirectory() as temp_directory:
        files = []
        for i in range(5):
            tmpfile = tempfile.NamedTemporaryFile(suffix='.fits', delete=False, dir=temp_directory)
            files.append(tmpfile.name)
        paths = select_data(temp_directory, 'any', ['.fits'], FakeImage, 0, header_keys, type_translator)
        assert len(paths) == 5
        paths = select_data(temp_directory, 'lampflat', ['.fits'], FakeImage, 0, header_keys, type_translator)
        assert len(paths) == 5
        paths = select_data(temp_directory, 'wavecal', ['.fits'], FakeImage, 0, header_keys, type_translator)
        assert len(paths) == 0
    pass


def test_organize_config():
    config = ConfigParser()
    config.read('echelle/tests/data/test_config.ini')
    runtime_context, data_class, extension, header_keys, type_translator = organize_config(config)
    assert data_class == 'echelle.images.Image'
    assert extension == 1
    print(type_translator)
    assert type_translator['LAMPFLAT'] == 'lampflat'
    assert type_translator['DOUBLE'] == 'wavecal'
    assert header_keys['read_noise'] == 'RDNOISE'
    assert type(runtime_context.final_wavelength_model) is dict
    assert type(runtime_context.ref_id) is int
    assert type(runtime_context.m0_range) is tuple


@mock.patch('echelle.main.parse_args')
def test_run(mock_args):
    with tempfile.TemporaryDirectory() as temp_directory:
        for i in range(5):
            tempfile.NamedTemporaryFile(suffix='.wrong', delete=False, dir=temp_directory)
        mock_args.return_value = type('', (), {'input_dir': temp_directory, 'output_dir': temp_directory,
                                               'config_file': 'echelle/tests/data/test_config.ini',
                                               'frame_type': 'any'})
        run()
    assert True


@mock.patch('echelle.main.add_data_to_db', return_value=None)
@mock.patch('echelle.main.format_db_info', return_value=None)
@mock.patch('echelle.main.parse_args')
def test_reduce_data_does_not_error(mock_args, mock_format, mock_add):
    config = ConfigParser()
    config.read('echelle/tests/data/test_config.ini')
    config.set('stages', 'lampflat', '[]')
    config.set('data', 'data_class', 'echelle.tests.utils.FakeImage')
    with tempfile.TemporaryDirectory() as temp_directory:
        tmp = tempfile.NamedTemporaryFile(suffix='.wrong', delete=False, dir=temp_directory)
        mock_args.return_value = type('', (), {'data_paths': [tmp], 'output_dir': temp_directory,
                                               'frame_type': 'any', 'fpack': False})
        reduce_data(config=config)
    assert True