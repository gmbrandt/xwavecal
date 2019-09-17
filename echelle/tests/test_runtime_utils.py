import datetime
import pytest
import mock

from banzai.tests.utils import FakeContext
from banzai_nres.utils import runtime_utils
from banzai_nres.tests.utils import FakeImage


def test_get_telescope_filename():
    image = FakeImage()
    image.header['TELESCOP'] = 'nres01'
    assert runtime_utils.get_telescope_filename(image) == 'nrs01'


def test_min_date_auto_fills():
    expected_max_date = datetime.datetime.utcnow()
    expected_min_date = expected_max_date - datetime.timedelta(hours=24)
    min_date, max_date = runtime_utils.get_reduction_date_window(runtime_context=None)
    assert min_date > expected_min_date - datetime.timedelta(seconds=10)
    assert min_date < expected_min_date + datetime.timedelta(seconds=10)
    assert max_date > expected_max_date - datetime.timedelta(seconds=10)
    assert max_date < expected_max_date + datetime.timedelta(seconds=10)


def test_min_date_is_correct_for_given_max_date():
    fake_context = FakeContext()
    fake_context.max_date = datetime.datetime.utcnow()
    min_date, max_date = runtime_utils.get_reduction_date_window(runtime_context=fake_context)
    assert min_date == fake_context.max_date - datetime.timedelta(hours=24)


def test_max_date_is_correct_for_given_min_date():
    fake_context = FakeContext()
    expected_max_date = datetime.datetime.utcnow()
    fake_context.min_date = datetime.datetime.utcnow()
    min_date, max_date = runtime_utils.get_reduction_date_window(runtime_context=fake_context)
    assert max_date > expected_max_date - datetime.timedelta(seconds=10)
    assert max_date < expected_max_date + datetime.timedelta(seconds=10)


def test_raise_exception_if_min_date_greater_than_max_date():
    fake_context = FakeContext()
    fake_context.min_date = datetime.datetime.utcnow()
    fake_context.max_date = fake_context.min_date - datetime.timedelta(hours=1)
    with pytest.raises(Exception):
        runtime_utils.get_reduction_date_window(runtime_context=fake_context)


def test_frame_types_auto_fill():
    assert runtime_utils.get_frame_types(None, default_frames_to_reduce=['TYPE', 'TYPE2']) == ['TYPE', 'TYPE2']


def test_get_single_frame_types():
    fake_context = FakeContext()
    fake_context.frame_type = 'TYPE'
    assert runtime_utils.get_frame_types(fake_context, None) == [fake_context.frame_type]


def test_validate_good_path():
    fake_context = FakeContext()
    raw_path = runtime_utils.validate_raw_path(runtime_context=fake_context, raw_path='test/raw')
    assert raw_path == 'test/raw'


@mock.patch('banzai_nres.utils.runtime_utils.get_raw_path', return_value='test/raw')
def test_complete_raw_path(mock_get_raw_path):
    fake_context = FakeContext()
    raw_path = runtime_utils.validate_raw_path(raw_path='test/', runtime_context=fake_context)
    assert raw_path == 'test/raw'
