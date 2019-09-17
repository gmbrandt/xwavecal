from banzai_nres.images import NRESImage
from banzai.tests.utils import FakeContext
import mock

from banzai.utils.import_utils import import_attribute


def test_frame_class_is_set():
    import banzai_nres.settings
    import banzai.settings as settings
    assert settings.FRAME_CLASS == 'banzai_nres.images.NRESImage'


@mock.patch('banzai.images.Image._init_instrument_info')
def test_frame_class_is_maintained(mock_instrument):
    mock_instrument.return_value = None, None, None
    import banzai_nres.settings
    import banzai.settings as settings
    frame_class = import_attribute(settings.FRAME_CLASS)
    nres_image = NRESImage(FakeContext(), header={'OBJECTS': 'none&tung&none'})
    frame = frame_class(FakeContext(), header={'OBJECTS': 'none&tung&none'})
    assert type(frame) is type(nres_image)
