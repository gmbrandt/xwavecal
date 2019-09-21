from configparser import ConfigParser
import tempfile
import os
from glob import glob
import pytest

from echelle.main import reduce_data
import logging as logger


@pytest.mark.e2e
def test_reduce_data():
    # TODO this test should be replaced with a proper pytest fixture which makes the tempfile
    #  once per session. Then we can make independent tests for lampflat and arc creation.
    with tempfile.TemporaryDirectory() as temp_directory:
        args = type('test_args', (), {'input_dir': 'echelle/tests/data/nres_test_data/',
                                      'output_dir': temp_directory, 'fpack': True})
        config = ConfigParser()
        config.read('echelle/tests/data/test_config.ini')
        config.set('reduction', 'database_path', '"' + os.path.join(temp_directory, 'test.db') + '"')
        data_paths = glob('echelle/tests/data/nres_test_data/*w00*.fits*')
        reduce_data(data_paths, args=args, config=config)
        # check that the lampflat, traces and blaze are in the database via query for match
        # check that the correct number of traces exists.

        #data_paths = glob('echelle/tests/data/nres_test_data/*a00*.fits*')
        #reduce_data(data_paths, args=args, config=config)
        assert True
