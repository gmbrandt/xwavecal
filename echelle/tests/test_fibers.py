import mock
import numpy as np
import pytest
from astropy.table import Table

from echelle.utils.fiber_utils import fiber_states_from_header, fibers_state_to_filename, \
                                          lit_fibers, lit_wavecal_fibers
from echelle.fibers import IdentifyFibers
from echelle.images import Image

from echelle.tests.utils import FakeImage, FakeContext


def test_creation_from_header():
    header = {'OBJECTS': 'targ&ThAr&none'}
    assert (1, 1, 0) == fiber_states_from_header(header)

    header = {'OBJECTS': 'none&ThAr&targ'}
    assert (0, 1, 1) == fiber_states_from_header(header)

    header = {'OBJECTS': 'none&ThAr&none'}
    assert (0, 1, 0) == fiber_states_from_header(header)


def test_fiber_state_to_filename():
    image = Image(header={'OBJECTS': 'tung&tung&none'})
    assert fibers_state_to_filename(image) == '110'


class TestIdentifyFibers:
    def test_build_fiber_column_double(self):
        # frames with two arc fibers lit.
        image = FakeImage()
        image.fiber0_wavecal, image.fiber1_wavecal, image.fiber2_wavecal = 0, 1, 1
        image.fiber0_lit, image.fiber1_lit, image.fiber2_lit = 0, 1, 1
        spectrum = {'id': np.arange(5)}
        matched_ids = [2, 3]
        fiber_ids = IdentifyFibers.build_fiber_column(matched_ids, image, spectrum)
        assert np.allclose(fiber_ids, [1, 2, 1, 2, 1])
        image.fiber0_wavecal, image.fiber1_wavecal, image.fiber2_wavecal = 1, 1, 0
        fiber_ids = IdentifyFibers.build_fiber_column(matched_ids, image, spectrum)
        assert np.allclose(fiber_ids, [0, 1, 0, 1, 0])
        spectrum = {'id': np.arange(6)}
        fiber_ids = IdentifyFibers.build_fiber_column(matched_ids, image, spectrum)
        assert np.allclose(fiber_ids, [0, 1, 0, 1, 0, 1])

    def test_build_fiber_column_single(self):
        # frames with one arc fiber lit.
        image = FakeImage()
        image.fiber0_wavecal, image.fiber1_wavecal, image.fiber2_wavecal = 0, 1, 0
        image.fiber0_lit, image.fiber1_lit, image.fiber2_lit = 0, 1, 1
        spectrum = {'id': np.arange(5)}
        matched_ids = [2]
        fiber_ids = IdentifyFibers.build_fiber_column(matched_ids, image, spectrum)
        assert np.allclose(fiber_ids, [1, 2, 1, 2, 1])
        image.fiber0_wavecal, image.fiber1_wavecal, image.fiber2_wavecal = 0, 0, 1
        fiber_ids = IdentifyFibers.build_fiber_column(matched_ids, image, spectrum)
        assert np.allclose(fiber_ids, [2, 1, 2, 1, 2])

    def test_build_ref_id_column_double(self):
        # frames with two arc fibers lit.
        image = FakeImage()
        image.fiber0_wavecal, image.fiber1_wavecal, image.fiber2_wavecal = 0, 1, 1
        image.fiber0_lit, image.fiber1_lit, image.fiber2_lit = 0, 1, 1
        spectrum = {'id': np.arange(5)}
        matched_ids = [2, 3]
        ref_ids = IdentifyFibers.build_ref_id_column(matched_ids, image, spectrum=spectrum, ref_id=50)
        assert np.allclose(ref_ids, [49, 49, 50, 50, 51])

    def test_build_ref_id_column_single(self):
        # frames with one arc fiber lit.
        image = FakeImage()
        image.fiber0_wavecal, image.fiber1_wavecal, image.fiber2_wavecal = 0, 0, 1
        image.fiber0_lit, image.fiber1_lit, image.fiber2_lit = 0, 1, 1
        spectrum = {'id': np.arange(5)}
        matched_ids = [3]
        ref_ids = IdentifyFibers.build_ref_id_column(matched_ids, image, spectrum=spectrum, ref_id=50)
        assert np.allclose(ref_ids, [49, 49, 50, 50, 51])
        matched_ids = [2]
        image.fiber0_wavecal, image.fiber1_wavecal, image.fiber2_wavecal = 0, 1, 0
        ref_ids = IdentifyFibers.build_ref_id_column(matched_ids, image, spectrum=spectrum, ref_id=50)
        assert np.allclose(ref_ids, [49, 49, 50, 50, 51])

    def test_calibration_type(self):
        assert IdentifyFibers(FakeContext()).calibration_type == 'TRACE'

    @mock.patch('echelle.fibers.IdentifyFibers.get_calibration_filename', return_value=None)
    @mock.patch('echelle.fibers.IdentifyFibers.on_missing_master_calibration')
    def test_do_stage_aborts_on_missing_cal(self, mock_warn, mock_cal):
        assert 'image' == IdentifyFibers(FakeContext()).do_stage(image='image')

    @mock.patch('echelle.fibers.IdentifyFibers.apply_master_calibration')
    @mock.patch('echelle.fibers.IdentifyFibers.get_calibration_filename', return_value='/path/')
    def test_do_stage(self, fake_cal, mock_apply_cal):
        IdentifyFibers(FakeContext()).do_stage(image='image')
        mock_apply_cal.assert_called_with('image', '/path/')

    @pytest.mark.integration
    @mock.patch('echelle.fibers.np.genfromtxt')
    def test_apply_master_calibration(self, mock_load):
        image = FakeImage()
        context = FakeContext()
        context.ref_id = 1
        image.fiber0_wavecal, image.fiber1_wavecal, image.fiber2_wavecal = 0, 1, 0
        image.fiber0_lit, image.fiber1_lit, image.fiber2_lit = 0, 1, 1
        image.set_header_val('read_noise', 0)
        spec = Table({'id': [0, 1, 2], 'flux': 100 * np.random.random((3, 30))})
        mock_load.return_value = spec['flux'].data[1].reshape(-1, 1)
        image.data_tables = {context.box_spectrum_name: spec}
        image = IdentifyFibers(context).apply_master_calibration(image, '')
        spec = image.data_tables[context.box_spectrum_name]
        assert np.allclose(spec['ref_id'], [0, 1, 1])
        assert np.allclose(spec['fiber'], [2, 1, 2])
    

@mock.patch('numpy.genfromtxt', return_value=np.array([np.ones(6), np.ones(6)]).T)
def test_construct_single_fiber_template(mock_template):
    template = IdentifyFibers.construct_single_fiber_template(None, 2)
    assert np.allclose(template, np.array([np.ones(5), np.zeros(5), np.ones(5)]))


@mock.patch('echelle.fibers.correlate2d')
def test_identify_matching_orders(fake_correlate):
    fake_correlate.return_value = np.array([[0, 1, 0],
                                            [0, 1, 2],
                                            [0, 0, 0]])
    assert np.allclose(IdentifyFibers.identify_matching_orders(None, None, 2), [1, 0])


def test_get_lit_fibers():
    image = FakeImage()
    image.fiber0_lit, image.fiber1_lit, image.fiber2_lit = 1, 1, 0
    image.fiber0_wavecal, image.fiber1_wavecal, image.fiber2_wavecal = 0, 1, 1
    assert np.allclose(lit_fibers(image), [0, 1])
    assert np.allclose(lit_wavecal_fibers(image), [1, 2])

