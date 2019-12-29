import tempfile
from configparser import ConfigParser
import mock
import pytest

from xwavecal.main import select_data, run, reduce_data, organize_config, make_output_path, import_obj
from xwavecal.main import RuntimeContext
from xwavecal.tests.utils import FakeImage


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
    config.read('xwavecal/tests/data/test_config.ini')
    runtime_context, data_class, extension, header_keys, type_keys = organize_config(config)
    assert data_class == 'xwavecal.images.Image'
    assert extension == 1
    assert type_keys['LAMPFLAT'] == 'lampflat'
    assert type_keys['DOUBLE'] == 'wavecal'
    assert header_keys['read_noise'] == 'RDNOISE'
    assert type(runtime_context.final_wavelength_model) is dict
    assert type(runtime_context.ref_id) is int
    assert type(runtime_context.m0_range) is tuple
    assert runtime_context.time_format == '%Y-%m-%dT%H:%M:%S.%f'


@mock.patch('xwavecal.main.parse_args')
def test_run(mock_args):
    with tempfile.TemporaryDirectory() as temp_directory:
        for i in range(5):
            tempfile.NamedTemporaryFile(suffix='.wrong', delete=False, dir=temp_directory)
        mock_args.return_value = type('', (), {'input_dir': temp_directory, 'output_dir': temp_directory,
                                               'config_file': 'xwavecal/tests/data/test_config.ini',
                                               'frame_type': 'any'})
        run()
    assert True


@mock.patch('xwavecal.main.order_data')
@mock.patch('xwavecal.main.add_data_to_db', return_value=None)
@mock.patch('xwavecal.main.format_db_info', return_value=None)
@mock.patch('xwavecal.main.parse_args')
def test_reduce_data_does_not_err(mock_args, mock_format, mock_add, mock_order):
    config = ConfigParser()
    config.read('xwavecal/tests/data/test_config.ini')
    config.set('stages', 'lampflat', '[]')
    config.set('data', 'data_class', 'xwavecal.tests.utils.FakeImage')
    with tempfile.TemporaryDirectory() as temp_directory:
        tmp = tempfile.NamedTemporaryFile(suffix='.wrong', delete=False, dir=temp_directory)
        mock_order.return_value = [tmp]
        mock_args.return_value = type('', (), {'data_paths': [tmp], 'output_dir': temp_directory,
                                               'frame_type': 'any', 'fpack': False})
        reduce_data(config=config)
    assert True


@mock.patch('xwavecal.main.order_data', return_value=[])
@mock.patch('xwavecal.main.organize_config', return_value=(None, 'xwavecal.images.Image', None, {}, {}))
@mock.patch('configparser.ConfigParser.read')
@mock.patch('xwavecal.main.parse_args')
def test_reduce_data_calls_config(mock_args, mock_config, mock_organize, mock_order):
    mock_args.return_value = type('', (), {'config_file': 'file', 'data_paths': 'path'})
    reduce_data()
    mock_config.assert_called_with('file')
    assert True


def test_make_output_path():
    data = FakeImage()
    path = make_output_path('outdir', data)
    data.set_header_val('site_name', '(test, )')
    path2 = make_output_path('outdir', data)
    assert path == 'outdir/test_nres03_20190410_0077_lampflat_011.fits'
    assert path2 == 'outdir/_test_nres03_20190410_0077_lampflat_011.fits'


def test_import_from_string():
    numpyones = import_obj('numpy.ones')
    assert numpyones((5, 5)).shape == (5, 5)


def test_runtime_context_raises_attr_error():
    with pytest.raises(AttributeError):
        RuntimeContext({}).missing_attribute
