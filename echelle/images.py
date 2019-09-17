from astropy.io import fits
from astropy.table import Table

from banzai.images import Image
from banzai_nres.utils.fiber_utils import fiber_states_from_header, wavecal_fibers_from_header
from banzai_nres.utils.fiber_utils import lit_wavecal_fibers, lit_fibers

import banzai_nres.settings as nres_settings  # import to override banzai settings
from banzai import settings
from banzai_nres.utils import fits_utils, db_utils
from banzai import dbs


class NRESImage(Image):
    def __init__(self, runtime_context, filename=None, data=None, data_tables=None,
                 header=None, extension_headers=None, bpm=None):
        super(NRESImage, self).__init__(runtime_context, filename=filename, data=data, data_tables=data_tables,
                                        header=header, extension_headers=extension_headers, bpm=bpm)
        self.trace = None
        self.rectified_2d_spectrum = None

        self.fiber0_lit, self.fiber1_lit, self.fiber2_lit = fiber_states_from_header(self.header)
        self.fiber0_wavecal, self.fiber1_wavecal, self.fiber2_wavecal = wavecal_fibers_from_header(self.header)

        self.wavelength_solution = {}

    def num_lit_fibers(self):
        return len(lit_fibers(self))

    def num_wavecal_fibers(self):
        return len(lit_wavecal_fibers(self))


class ImageBase(object):
    """
    Image like object with .write and .load methods. Used in trace and Blaze calibrations
    """
    def __init__(self, data=None, table_name=None, filepath=None, header=None, image=None, obstype='None'):
        if header is None:
            header = {}
        self.header = header
        self.filepath = filepath
        self.data = Table(data)
        self.table_name = table_name

        # banzai.images.Image attributes necessary for fits_utils.writeto
        self.obstype = obstype
        self.dateobs = getattr(image, 'dateobs', None)
        self.datecreated = getattr(image, 'datecreated', None)
        self.instrument = getattr(image, 'instrument', None)
        self.is_master = getattr(image, 'is_master', False)
        self.is_bad = getattr(image, 'is_bad', False)
        self.attributes = settings.CALIBRATION_SET_CRITERIA.get(self.obstype, {})
        for attribute in self.attributes:
            setattr(self, attribute, getattr(image, attribute, None))

    def write(self, runtime_context=None, update_db=True):
        hdu = fits.BinTableHDU(self.data, name=self.table_name, header=fits.Header(self.header))
        hdu_list = fits.HDUList([fits.PrimaryHDU(), hdu])
        self._update_filepath(runtime_context)
        fits_utils.writeto(hdu_list=hdu_list, filepath=self.filepath,
                           fpack=getattr(runtime_context, 'fpack', False),
                           overwrite=True, output_verify='fix+warn')
        if update_db:
            dbs.save_calibration_info(self.filepath, image=self,
                                      db_address=runtime_context.db_address)
            if runtime_context.post_to_archive:
                db_utils.post_to_archive(self.filepath)

    def _update_filepath(self, runtime_context):
        if getattr(runtime_context, 'fpack', False) and not self.filepath.endswith('.fz'):
            self.filepath += '.fz'

    @classmethod
    def load(cls, path, extension_name):
        hdu_list = fits.open(path)
        return cls(data=hdu_list[extension_name].data, header=hdu_list[extension_name].header)
