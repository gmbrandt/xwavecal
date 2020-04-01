import mock
import numpy as np
import pytest
from astropy.table import Table

from xwavecal.utils.fiber_utils import fiber_states_from_header, fibers_state_to_filename, \
                                          lit_fibers, lit_wavecal_fibers, wavecal_fibers_from_header
from xwavecal.fibers import IdentifyFibers, MakeFiberTemplate
from xwavecal.images import Image, DataProduct

from xwavecal.tests.utils import FakeImage, FakeContext


def test_creation_from_header():
    assert (1, 1, 0) == fiber_states_from_header('targ&ThAr&none')
    assert (0, 1, 1) == fiber_states_from_header('none&ThAr&targ')
    assert (0, 1, 0) == fiber_states_from_header('none&ThAr&none')


def test_wavecal_fibers_from_header():
    assert (0, 1, 0) == wavecal_fibers_from_header('targ&ThAr&none', 'thar')
    assert (0, 1, 1) == wavecal_fibers_from_header('none&ThAr&ThAr', 'thar')


def test_fiber_state_to_filename():
    image = Image(header={'fiber_state': 'tung&tung&none'})
    assert fibers_state_to_filename(image) == '110'


class TestIdentifyFibers:
    CONTEXT = FakeContext()

    def test_build_fiber_column(self):
        wavecal_fibers = np.array([1, 2])
        fibers = np.array([0, 1, 2])
        num_traces = 7
        for matched_ids in [[2], [3], [4]]:
            assert IdentifyFibers.build_fiber_column(matched_ids, fibers, wavecal_fibers, num_traces)[matched_ids[0]] == 1
            assert IdentifyFibers.build_fiber_column(matched_ids, fibers, wavecal_fibers, num_traces,
                                                     low_fiber_first=False)[matched_ids[0]] == 2

        assert np.allclose(IdentifyFibers.build_fiber_column([3], fibers, wavecal_fibers, num_traces),
                           [1, 2, 0, 1, 2, 0, 1])
        assert np.allclose(IdentifyFibers.build_fiber_column([3], fibers, wavecal_fibers, num_traces, low_fiber_first=False),
                           [2, 1, 0, 2, 1, 0, 2])

    def test_build_fiber_column_multi_match(self):
        wavecal_fibers = np.array([1, 2])
        fibers = np.array([0, 1, 2])
        num_traces = 7
        expected_fiber_designations = [1, 2, 0, 1, 2, 0, 1]

        assert np.allclose(IdentifyFibers.build_fiber_column([3, 4], fibers, wavecal_fibers, num_traces),
                           expected_fiber_designations)
        assert np.allclose(IdentifyFibers.build_fiber_column([4, 3], fibers, wavecal_fibers, num_traces),
                           expected_fiber_designations)

    def test_build_fiber_column_single_fiber(self):
        wavecal_fibers = np.array([1])
        fibers = np.array([1])
        num_traces = 7
        for matched_ids in [[2], [3], [4]]:
            assert IdentifyFibers.build_fiber_column(matched_ids, fibers, wavecal_fibers, num_traces)[matched_ids[0]] == 1
            assert IdentifyFibers.build_fiber_column(matched_ids, fibers, wavecal_fibers, num_traces,
                                                     low_fiber_first=False)[matched_ids[0]] == 1

        assert np.allclose(IdentifyFibers.build_fiber_column([3], fibers, wavecal_fibers, num_traces), np.ones(7))
        assert np.allclose(IdentifyFibers.build_fiber_column([3], fibers, wavecal_fibers, num_traces, low_fiber_first=False), np.ones(7))

    def test_build_ref_id_column(self):
        fiber_ids = np.array([1, 2, 1, 2, 1])
        ref_ids = IdentifyFibers.build_ref_id_column([2, 3], fiber_ids, ref_id=50)
        assert np.allclose(ref_ids, [49, 49, 50, 50, 51])
        fiber_ids = np.array([0, 1, 2, 0, 1])
        ref_ids = IdentifyFibers.build_ref_id_column([1, 2], fiber_ids, ref_id=50)
        assert np.allclose(ref_ids, [50, 50, 50, 51, 51])
        fiber_ids = np.array([0, 2, 1, 0, 2, 1])
        ref_ids = IdentifyFibers.build_ref_id_column([1, 2], fiber_ids, ref_id=50, low_fiber_first=False)
        assert np.allclose(ref_ids, [49, 50, 50, 50, 51, 51])
        fiber_ids = np.array([0, 0, 0, 0, 0, 0])
        ref_ids = IdentifyFibers.build_ref_id_column([1], fiber_ids, ref_id=50, low_fiber_first=False)
        assert np.allclose(ref_ids, [49, 50, 51, 52, 53, 54])

    def test_build_ref_id_column_reversed_match_order(self):
        fiber_ids = np.array([1, 2, 1, 2, 1])
        ref_ids = IdentifyFibers.build_ref_id_column([3, 2], fiber_ids, ref_id=50)
        assert np.allclose(ref_ids, [49, 49, 50, 50, 51])

    def test_calibration_type(self):
        assert IdentifyFibers(self.CONTEXT).calibration_type == 'FIBERS'

    @mock.patch('xwavecal.fibers.IdentifyFibers.get_calibration_filename', return_value=None)
    @mock.patch('xwavecal.fibers.IdentifyFibers.on_missing_master_calibration')
    def test_do_stage_aborts_on_missing_cal(self, mock_warn, mock_cal):
        assert 'image' == IdentifyFibers(self.CONTEXT).do_stage(image='image')

    @mock.patch('xwavecal.fibers.IdentifyFibers.get_calibration_filename', return_value='')
    @mock.patch('os.path.exists', return_value=True)
    def test_do_stage_aborts_on_no_wavecal_fibers(self, mock_warn, mock_exists):
        image = FakeImage()
        image.data_tables[self.CONTEXT.main_spectrum_name] = {'flux': [1, 2]}
        image.fiber0_wavecal, image.fiber1_wavecal, image.fiber2_wavecal = 0, 0, 0
        assert image == IdentifyFibers(self.CONTEXT).do_stage(image=image)

    @mock.patch('xwavecal.fibers.IdentifyFibers.get_calibration_filename', return_value='')
    @mock.patch('os.path.exists', return_value=True)
    def test_do_stage_aborts_on_length_zero_spectrum(self, mock_warn, mock_exists):
        image = FakeImage()
        image.data_tables[self.CONTEXT.main_spectrum_name] = {'flux': []}
        assert image == IdentifyFibers(FakeContext()).do_stage(image=image)

    @mock.patch('os.path.exists', return_value=True)
    @mock.patch('xwavecal.fibers.IdentifyFibers.apply_master_calibration')
    @mock.patch('xwavecal.fibers.IdentifyFibers.get_calibration_filename', return_value='/path/')
    def test_do_stage(self, fake_cal, mock_apply_cal, mock_os):
        IdentifyFibers().do_stage(image='image')
        mock_apply_cal.assert_called_with('image', '/path/')

    @pytest.mark.integration
    @mock.patch('xwavecal.images.fits.open')
    def test_apply_master_calibration(self, mock_load):
        image = FakeImage()
        context = FakeContext()
        context.ref_id = 1
        image.fiber0_wavecal, image.fiber1_wavecal, image.fiber2_wavecal = 0, 1, 0
        image.fiber0_lit, image.fiber1_lit, image.fiber2_lit = 0, 1, 1
        image.set_header_val('read_noise', 0)
        spec = Table({'id': [0, 1, 2], 'flux': 100 * np.random.random((3, 30)),
                      'fiber': [1, 1, 1], 'ref_id': [0, 0, 0]})
        mock_load.return_value = {'fibers': DataProduct(data=spec[1])}  # fake fiber template
        image.data_tables = {context.main_spectrum_name: spec}
        image = IdentifyFibers(context).apply_master_calibration(image, '')
        spec = image.data_tables[context.main_spectrum_name]
        assert np.allclose(spec['ref_id'], [0, 1, 1])
        assert np.allclose(spec['fiber'], [2, 1, 2])

    @pytest.mark.integration
    @mock.patch('xwavecal.images.fits.open')
    def test_apply_master_calibration_single_fiber(self, mock_load):
        image = FakeImage()
        context = FakeContext()
        context.ref_id = 1
        image.fiber0_wavecal, image.fiber1_wavecal, image.fiber2_wavecal = 0, 1, 0
        image.fiber0_lit, image.fiber1_lit, image.fiber2_lit = 0, 1, 0
        image.set_header_val('read_noise', 0)
        spec = Table({'id': [0, 1, 2], 'flux': 100 * np.random.random((3, 30)),
                      'fiber': [1, 1, 1], 'ref_id': [0, 0, 0]})
        mock_load.return_value = {'fibers': DataProduct(data=spec[1])}  # fake fiber template
        image.data_tables = {context.main_spectrum_name: spec}
        image = IdentifyFibers(context).apply_master_calibration(image, '')
        spec = image.data_tables[context.main_spectrum_name]
        assert np.allclose(spec['ref_id'], [0, 1, 2])
        assert np.allclose(spec['fiber'], [1, 1, 1])


class TestMakeFiberTemplate:
    def test_do_stage(self):
        context = FakeContext()
        context.template_trace_id = 2
        image = FakeImage()
        spec = Table({'id': [0, 1, 2, 3, 4], 'flux': 100 * np.random.random((5, 30))})
        image.data_tables = {context.main_spectrum_name: spec}
        image, template = MakeFiberTemplate(context).do_stage(image)
        assert template.get_header_val('type') == MakeFiberTemplate(context).calibration_type.lower()
        assert image.get_header_val('type') != MakeFiberTemplate(context).calibration_type.lower()
        assert np.allclose(spec[[0, 2, 4]]['flux'], template.data['flux'])


@mock.patch('xwavecal.images.fits.open')
def test_construct_single_fiber_template(mock_load):
    mock_load.return_value = {'fibers': DataProduct(data=Table({'flux': np.array([np.ones(5), np.ones(5)])}))}
    template = IdentifyFibers.construct_single_fiber_template(None, 'fibers', 2)
    assert np.allclose(template, np.array([np.ones(5), np.zeros(5), np.ones(5)]))
    mock_load.return_value = {'fibers': DataProduct(data=Table({'flux': np.zeros((3, 10))}))}
    template = IdentifyFibers.construct_single_fiber_template(None, 'fibers', 3)
    assert np.allclose(template.shape, (7, 10))
    template = IdentifyFibers.construct_single_fiber_template(None, 'fibers', 1)
    assert np.allclose(template.shape, (3, 10))


@mock.patch('xwavecal.fibers.correlate2d')
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

