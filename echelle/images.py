from astropy.io import fits
from astropy.table import Table

from echelle.utils.fiber_utils import fiber_states_from_header, wavecal_fibers_from_header
from echelle.utils.fiber_utils import lit_wavecal_fibers, lit_fibers
from echelle.utils import fits_utils

import logging as logger


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

    def get_header_val(self, key):
        if self.translator is not None:
            key = self.translator[key]
            return self.translator.type_keys.get(self.header[key], self.header[key])
        return self.header[key]

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
    # TODO THIS IS NRES SPECIFIC because of fiber_states_from_header
    def __init__(self, filepath=None, data=None, header=None, data_tables=None, translator=None, trace=None,
                 data_name=None, ivar=None):
        super(Image, self).__init__(filepath=filepath, data=data, header=header, translator=translator, data_name=data_name)
        if data_tables is None:
            data_tables = {}
        self.data_tables = data_tables
        self.data = data
        self.ivar = ivar
        self.filepath = filepath
        self.trace = trace
        self.rectified_2d_spectrum = None
        self.fiber0_lit, self.fiber1_lit, self.fiber2_lit = fiber_states_from_header(self.header)
        self.fiber0_wavecal, self.fiber1_wavecal, self.fiber2_wavecal = wavecal_fibers_from_header(self.header)

        self.wavelength_solution = {}

    def write(self, fpack=False):
        if self.data_name is None:
            hdus = [fits.PrimaryHDU(data=self.data, header=fits.Header(self.header))]
        else:
            hdus = [fits.PrimaryHDU(), fits.ImageHDU(data=self.data, header=fits.Header(self.header), name=self.data_name)]
        if self.ivar is not None:
            hdus.append(fits.ImageHDU(data=self.ivar, header=None, name='ivar'))

        table_hdus = [fits.BinTableHDU(table, name=name) for name, table in self.data_tables.items()]
        hdu_list = fits.HDUList([*hdus, *table_hdus])
        self._update_filepath(fpack)
        logger.info('Writing file to {filepath}'.format(filepath=self.filepath))
        fits_utils.writeto(hdu_list=hdu_list, filepath=self.filepath, fpack=fpack,
                           overwrite=True, output_verify='fix+warn')

    def num_lit_fibers(self):
        return len(lit_fibers(self))

    def num_wavecal_fibers(self):
        return len(lit_wavecal_fibers(self))