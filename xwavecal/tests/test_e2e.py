from configparser import ConfigParser
import tempfile
import os
import mock
import numpy as np
from glob import glob
import pytest
from astropy.io import fits
from astropy.table import Table

from xwavecal.main import reduce_data, organize_config

@pytest.mark.e2e
@mock.patch('xwavecal.fibers.IdentifyFibers.get_calibration_filename',
            return_value='xwavecal/tests/data/nres_test_data/cpt_nres03_20190405_0014_fibers_011.fits')
def test_reduce_data(mock_arc_template):
    # TODO this test should be replaced with a proper pytest fixture which makes the tempfile
    #  once per session. Then we can make independent tests for lampflat and arc creation.
    with tempfile.TemporaryDirectory() as temp_directory:
        args = type('test_args', (), {'input_dir': 'xwavecal/tests/data/nres_test_data/',
                                      'output_dir': temp_directory, 'fpack': False})
        config = ConfigParser()
        config.read('xwavecal/tests/data/test_config.ini')
        config.set('reduction', 'database_path', '"' + os.path.join(temp_directory, 'test.db') + '"')
        data_paths = glob('xwavecal/tests/data/nres_test_data/*w00*.fits*')
        reduce_data(data_paths, args=args, config=config)
        # check the reduced trace, lampflat and blaze.
        # fetch the context for easy access to configuration information
        context = organize_config(config)[0]
        # check that the lampflat, traces and blaze exist
        trace_file = glob(os.path.join(temp_directory, '*trace*'))[0]
        assert os.path.exists(glob(os.path.join(temp_directory, '*lampflat*'))[0])
        assert os.path.exists(trace_file)
        assert os.path.exists(glob(os.path.join(temp_directory, '*blaze*'))[0])
        check_traces(trace_file, context)

        # reduce the arc lamp file and make a wavelength solution.
        data_paths = glob('xwavecal/tests/data/nres_test_data/*a00*.fits*')
        reduce_data(data_paths, args=args, config=config)
        # check the wavelength calibration
        wavecal_file = glob(os.path.join(temp_directory, '*wavecal*'))[0]
        assert os.path.exists(wavecal_file)
        check_wavelength_calibration(wavecal_file, context)


def check_wavelength_calibration(fname, context):
    with fits.open(fname) as wcs:
        # assert that both fibers have at least 1000 lines with residuals less than 0.004 Angstrom.
        lines = Table(wcs[context.emission_lines_table_name].data)
        both_fibers = [np.isclose(lines['fiber'], 1), np.isclose(lines['fiber'], 2)]
        for fiber in both_fibers:
            assert 1000 < np.count_nonzero(np.abs(lines[fiber]['wavelength'] - lines[fiber]['reference_wavelength']) < 0.004)
        # assert that both fibers have >25 overlaps with peaks>6 and >15 marked as good.
        overlaps = Table(wcs[context.overlap_table_name].data)
        both_fibers = [np.isclose(overlaps['fiber'], 1), np.isclose(overlaps['fiber'], 2)]
        for fiber in both_fibers:
            assert len(overlaps[fiber]['peaks'] >= 6) > 25
            assert len(overlaps[fiber]['good'] >= 6) > 15


def check_traces(fname, context):
    with fits.open(fname) as trace:
        assert len(Table(trace[context.trace_table_name].data)) == 135