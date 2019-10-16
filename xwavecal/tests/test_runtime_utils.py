import pytest

from xwavecal.utils.runtime_utils import parse_args


def test_parse_args():
    args = parse_args(['--output-dir', 'output', '--data-paths', 'data', 'data2', '--config-file', 'config',
                       '--fpack', '--input-dir', 'input', '--frame-type', 'lampflat'])
    assert args.output_dir == 'output'
    assert args.data_paths == ['data', 'data2']
    assert args.config_file == 'config'
    assert args.input_dir == 'input'
    assert args.frame_type == 'lampflat'
    assert args.fpack is True


def test_parse_args_raises_error():
    with pytest.raises(ValueError):
        parse_args(['--config-file', 'config', '--output-dir', 'out'])
