import numpy as np
import abc
import logging as logger
import copy
from astropy.table import Table

from banzai_nres.utils import extract_utils
import banzai_nres.settings as nres_settings

from banzai.stages import Stage
from banzai.images import DataTable


class Extract(Stage):
    def __init__(self, runtime_context=None):
        super(Extract, self).__init__(runtime_context=runtime_context)

    @abc.abstractmethod
    def extract(self, rectified_2d_spectrum):
        """
         :param rectified_2d_spectrum: Dictionary
               Dictionary where the keys are the trace id's from trace.get_id(),
               where rectified_2d_spectrum[trace_id]['flux'] is a 2d float array (flux for the trace_id order).
               If half extraction window was 10, then rectified_2d_spectrum[trace_id]['flux']
               is 21 rows by 4096 columns (for a 4096 pixel wide image). One would sum this 2d
               spectrum along-columns to get a box extracted spectrum.
        :return:
        """
        pass

    @staticmethod
    def extract_order(twod_spectrum, weights=None):
        """
        :param twod_spectrum: 2d float array.
                              array of flux where the center of the order is the central row.
        :param weights: 2d array
                        weights, the same shape as twod_spectrum
        :return: 1d array.
                 weighted sum along the columns of twod_spectrum.
        """
        if weights is None:
            return np.sum(twod_spectrum, axis=0)
        else:
            return np.sum(twod_spectrum * weights, axis=0)

    def do_stage(self, image):
        return image


class BoxExtract(Extract):
    def __init__(self, runtime_context=None):
        super(Extract, self).__init__(runtime_context=runtime_context)
        self.extraction_half_window = nres_settings.BOX_EXTRACTION_HALF_WINDOW
        self.max_extraction_half_window = nres_settings.MAX_EXTRACTION_HALF_WINDOW

    def extract(self, rectified_2d_spectrum):
        extracted_spectrum_per_order = {'id': [], 'flux': [], 'pixel': []}
        for order_id in list(rectified_2d_spectrum.keys()):
            flux = self.extract_order(rectified_2d_spectrum[order_id]['flux'])
            extracted_spectrum_per_order['flux'].append(flux)
            extracted_spectrum_per_order['pixel'].append(np.arange(len(flux)))
            extracted_spectrum_per_order['id'].append(order_id)
        return Table(extracted_spectrum_per_order)

    def do_stage(self, image):
        logger.info('Box extracting spectrum', image=image)
        rectified_2d_spectrum = self._trim_rectified_2d_spectrum(image.rectified_2d_spectrum)
        spectrum = self.extract(rectified_2d_spectrum)
        table_name = nres_settings.BOX_SPECTRUM_NAME
        image.data_tables[table_name] = DataTable(data_table=spectrum,
                                                  name=table_name)
        return image

    def _trim_rectified_2d_spectrum(self, rectified_2d_spectrum):
        """
        :param rectified_2d_spectrum: Dictionary
               Dictionary where the keys are the trace id's from trace.get_id(),
               where rectified_2d_spectrum[trace_id]['flux'] is a 2d float array (flux for the trace_id order).
               If half extraction window was 10, then rectified_2d_spectrum[trace_id]['flux']
               is 21 rows by 4096 columns (for a 4096 pixel wide image). One would sum this 2d
               spectrum along-columns to get a box extracted spectrum.
        :return rectified_2d_spectrum: Dictionary
                 Same as input but trimmed so that each order's 2d spectrum only has 2 * extraction_half_window + 1 rows.
        NOTE: The output spectra per order have the center of the trace at the center of the spectrum. E.g. if
        extraction_half_window is 10, then the 2d spectra have 21 rows and the trace center (peak flux) lies at
        index 10 (indexing from 0).
        """
        trimmed_rectified_spectrum = copy.deepcopy(rectified_2d_spectrum)
        if self.extraction_half_window >= self.max_extraction_half_window:
            # short circuit
            logger.warning('Box extraction window was chosen to be >= the max extraction window defined in settings.py.'
                           ' Defaulting to the max extraction window.')
            return rectified_2d_spectrum
        trim = self.max_extraction_half_window - self.extraction_half_window
        for order_id in list(rectified_2d_spectrum.keys()):
            for data_type in list(rectified_2d_spectrum[order_id].keys()):
                trimmed_rectified_spectrum[order_id][data_type] = rectified_2d_spectrum[order_id][data_type][trim:-trim]
        return trimmed_rectified_spectrum


class RectifyTwodSpectrum(Stage):
    def __init__(self, runtime_context=None):
        super(RectifyTwodSpectrum, self).__init__(runtime_context=runtime_context)
        self.max_extraction_half_window = nres_settings.MAX_EXTRACTION_HALF_WINDOW

    def do_stage(self, image):
        logger.info('Rectifying the 2d spectrum', image=image)
        if image.trace is None:
            logger.error('Image has empty trace attribute. Aborting extraction.', image=image)
            image.is_bad = True
            image.rectified_2d_spectrum = {}
            return image
        rectified_2d_spectrum = extract_utils.rectify_orders(image.data, image.trace,
                                                             half_window=self.max_extraction_half_window)
        image.rectified_2d_spectrum = rectified_2d_spectrum
        return image
