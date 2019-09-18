from banzai.tests.utils import FakeContext
from echelle.flats import FlatStacker


def test_flatstacker_caltype():
    assert FlatStacker(FakeContext()).calibration_type == 'LAMPFLAT'
