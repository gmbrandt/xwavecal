from astropy.io import fits
from astropy.table import Table

from xwavecal.utils.fiber_utils import fiber_states_from_header, wavecal_fibers_from_header
from xwavecal.utils.fiber_utils import lit_wavecal_fibers, lit_fibers
from xwavecal.utils import fits_utils

import logging
logger = logging.getLogger(__name__)


class DataProduct(object):
    """
    Base class with .write and .load methods for all data products.
    """
    def __init__(self, data=None, data_name=None, filepath=None, header=None,
                 translator=None):
        if header is None:
            header = {}
        self.header = header
        self.filepath = filepath
        self.data = Table(data)
        self.data_name = data_name
        self.translator = translator

    def _get_header_vals(self, keys):
        keys = [keys] if isinstance(keys, str) else keys
        out = tuple(self.header[key] for key in keys)
        return out if len(out) > 1 else out[0]

    def get_header_val(self, key):
        if self.translator is not None:
            keys = self.translator[key]
            return self.translator.type_keys.get(self._get_header_vals(keys), self._get_header_vals(keys))
        return self._get_header_vals(key)

    def set_header_val(self, key, value):
        if self.translator is not None:
            key = self.translator[key]
        self.header[key] = value

    def write(self, fpack=False):
        hdu = fits.BinTableHDU(self.data, name=self.data_name, header=fits.Header(self.header))
        hdu_list = fits.HDUList([fits.PrimaryHDU(), hdu])
        self._update_filepath(fpack)
        logger.info('Writing file to {filepath}'.format(filepath=self.filepath))
        fits_utils.writeto(hdu_list=hdu_list, filepath=self.filepath, fpack=fpack,
                           overwrite=True, output_verify='fix+warn')

    def _update_filepath(self, fpack=False):
        if fpack and not self.filepath.endswith('.fz'):
            self.filepath += '.fz'

    @classmethod
    def load(cls, path, extension_name, translator=None):
        hdu_list = fits.open(path)
        return cls(data=hdu_list[extension_name].data, header=hdu_list[extension_name].header,
                   filepath=path, translator=translator)


class Image(DataProduct):
    """
    Class for two-dimensional data frames.
    """
    def __init__(self, filepath=None, data=None, header=None, data_tables=None, translator=None, trace=None,
                 data_name=None, ivar=None, wavelength_solution=None):
        super(Image, self).__init__(filepath=filepath, data=data, header=header, translator=translator, data_name=data_name)
        if data_tables is None:
            data_tables = {}
        if wavelength_solution is None:
            wavelength_solution = {}
        self.data_tables = data_tables
        self.data = data
        self.ivar = ivar
        self.filepath = filepath
        self.trace = trace
        self.rectified_2d_spectrum = None
        self.rectified_ivar = None
        self.fiber0_lit, self.fiber1_lit, self.fiber2_lit = fiber_states_from_header(self.get_header_val('fiber_state'))
        self.fiber0_wavecal, self.fiber1_wavecal, self.fiber2_wavecal = wavecal_fibers_from_header(self.get_header_val('fiber_state'))
        self.wavelength_solution = wavelength_solution

    def write(self, fpack=False):
        if self.data_name is None:
            hdus = [fits.PrimaryHDU(data=self.data, header=fits.Header(self.header))]
        else:
            hdus = [fits.PrimaryHDU(), fits.ImageHDU(data=self.data, header=fits.Header(self.header), name=self.data_name)]
        if self.ivar is not None:
            hdus.append(fits.ImageHDU(data=self.ivar, header=None, name='ivar'))

        table_hdus = [fits.BinTableHDU(table, name=name, header=table.meta.get('header'))
                      for name, table in self.data_tables.items()]
        hdu_list = fits.HDUList([*hdus, *table_hdus])
        self._update_filepath(fpack)
        logger.info('Writing file to {filepath}'.format(filepath=self.filepath))
        fits_utils.writeto(hdu_list=hdu_list, filepath=self.filepath, fpack=fpack,
                           overwrite=True, output_verify='fix+warn')

    def num_lit_fibers(self):
        return len(lit_fibers(self))

    def num_wavecal_fibers(self):
        return len(lit_wavecal_fibers(self))

    @classmethod
    def load(cls, path, extension_name, translator=None):
        hdu_list = fits.open(path)
        data_tables = {hdu.name: Table(hdu.data) for hdu in hdu_list if type(hdu) is fits.BinTableHDU}
        return cls(data=hdu_list[extension_name].data, header=hdu_list[extension_name].header,
                   filepath=path, translator=translator, data_tables=data_tables)


class SplitHeaderImage(Image):
    """
    Functionally the same as xwavecal.images.Image , except this assumes that the [0] extension
    of the raw data .fits file has information in the header integral to reduction. The [0].header is
    appended onto the primary_extension (e.g. [1], specified in the config.ini) header for the image.
    Example: This class is suitable for HARPS data.
    """
    def __init__(self, filepath=None, data=None, header=None, data_tables=None, translator=None, trace=None,
                 data_name=None, ivar=None):
        super(SplitHeaderImage, self).__init__(filepath=filepath, data=data, header=header, translator=translator,
                                               data_name=data_name, data_tables=data_tables, trace=trace, ivar=ivar)

    @classmethod
    def load(cls, path, extension_name, translator=None):
        hdu_list = fits.open(path)
        data_tables = {hdu.name: Table(hdu.data) for hdu in hdu_list if type(hdu) is fits.BinTableHDU}
        header = hdu_list[extension_name].header
        header.extend(hdu_list[0].header)
        return cls(data=hdu_list[extension_name].data, header=header, filepath=path, translator=translator,
                   data_tables=data_tables)