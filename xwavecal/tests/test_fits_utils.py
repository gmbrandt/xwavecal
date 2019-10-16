import mock
from astropy.io.fits import HDUList, PrimaryHDU
from tempfile import TemporaryDirectory

from xwavecal.utils.fits_utils import writeto
from xwavecal.utils.instrument_specific import parse_nres_region_keyword, parse_harps_region_keyword


class TmpDir(TemporaryDirectory):
    def __init__(self, suffix=None, prefix=None, dir=None):
        super(TmpDir, self).__init__(suffix=suffix, prefix=prefix, dir=dir)

    def __enter__(self):
        return 't3mp'


class TestWriteTo:
    @mock.patch('astropy.io.fits.HDUList.writeto')
    def test_writeto_no_fpack(self, mock_write):
        hdu_list = HDUList([PrimaryHDU()])
        writeto(hdu_list, 'path.fz', fpack=False, overwrite=True, output_verify='verify')
        mock_write.assert_called_with('path', overwrite=True, output_verify='verify')

    @mock.patch('shutil.move')
    @mock.patch('os.system')
    @mock.patch('astropy.io.fits.HDUList.writeto')
    @mock.patch('tempfile.TemporaryDirectory', new=TmpDir)
    def test_writeto_fpack(self, mock_write, mock_sys, mock_move):
        hdu_list = HDUList([PrimaryHDU()])
        writeto(hdu_list, 'path', fpack=True, overwrite=True, output_verify='verify', quant=64)
        mock_write.assert_called_with('t3mp/path', overwrite=True, output_verify='verify')
        mock_sys.assert_called_with('fpack -q 64 t3mp/path')
        mock_move.assert_called_with('t3mp/path.fz', 'path.fz')


def test_parse_region_keyword():
    assert (slice(0, 4096), slice(0, 4096)) == parse_nres_region_keyword('[1:4096,1:4096]')


def test_parse_harps_region_keyword():
    key = (4096, 4096, 50, 0, 50, 0, 4196, 4096)
    assert (slice(0, 4096), slice(50, 4146)) == parse_harps_region_keyword(key)
    key = (50, 0, 50, 0, 50, 0, 4196, 4096)
    assert (slice(0, 4096), slice(4146, 4196)) == parse_harps_region_keyword(key)
