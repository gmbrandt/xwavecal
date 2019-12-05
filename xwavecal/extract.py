import numpy as np
import copy
from astropy.table import Table

from xwavecal.utils import extract_utils
from xwavecal.stages import Stage


class BoxExtract(Stage):
    def __init__(self, runtime_context=None):
        super(BoxExtract, self).__init__(runtime_context=runtime_context)
        self.extraction_half_window = runtime_context.box_extraction_half_window
        self.max_extraction_half_window = runtime_context.max_extraction_half_window
        self.table_name = runtime_context.box_spectrum_name

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

    def extract(self, rectified_2d_spectrum, rectified_ivar):
        """
         :param rectified_2d_spectrum: Dictionary
                Dictionary where the keys are the trace id's from trace.get_id(),
                where rectified_2d_spectrum[trace_id]['val'] is a 2d float array (flux for the trace_id order).
                If half extraction window was 10, then rectified_2d_spectrum[trace_id]['flux']
                is 21 rows by 4096 columns (for a 4096 pixel wide image). One would sum this 2d
                spectrum along-columns to get a box extracted spectrum.
         :param rectified_ivar: Dictionary
                Same as rectified_2d_spectrum but ['val'] is the inverse variance.
        :return:
        """
        extracted_spectrum_per_order = {'id': [], 'flux': [], 'stderr': [], 'pixel': [],
                                        'fiber': [], 'ref_id': [], 'wavelength': []}
        for order_id in list(rectified_2d_spectrum.keys()):
            weights = self._weights(rectified_2d_spectrum[order_id]['val'], rectified_ivar[order_id]['val'])
            flux = self.extract_order(rectified_2d_spectrum[order_id]['val'], weights)
            stdvar = self.extract_order(np.power(rectified_ivar[order_id]['val'], -1,
                                                 where=~np.isclose(rectified_ivar[order_id]['val'], 0),
                                                 out=np.zeros_like(rectified_ivar[order_id]['val'])),
                                        safe_pow(weights, 2))
            extracted_spectrum_per_order['flux'].append(flux)
            extracted_spectrum_per_order['stderr'].append(np.sqrt(stdvar))
            extracted_spectrum_per_order['pixel'].append(np.arange(len(flux)))
            extracted_spectrum_per_order['id'].append(order_id)

            extracted_spectrum_per_order['wavelength'].append(np.ones_like(flux, dtype=float) * np.nan)
            extracted_spectrum_per_order['ref_id'].append(np.nan)
            extracted_spectrum_per_order['fiber'].append(np.nan)
        return Table(extracted_spectrum_per_order)

    def _weights(self, order_rect_spectrum, order_rect_ivar):
        return None

    def do_stage(self, image):
        self.logger.info('Extracting spectrum')
        rectified_2d_spectrum = self._trim_rectified_2d_spectrum(image.rectified_2d_spectrum)
        rectified_ivar = self._trim_rectified_2d_spectrum(image.rectified_ivar)
        spectrum = self.extract(rectified_2d_spectrum, rectified_ivar)
        image.data_tables[self.table_name] = Table(spectrum)
        return image

    def _trim_rectified_2d_spectrum(self, rectified_2d_spectrum):
        """
        :param rectified_2d_spectrum: Dictionary
               Dictionary where the keys are the trace id's from trace.get_id(),
               where rectified_2d_spectrum[trace_id]['val'] is a 2d float array (flux for the trace_id order).
               If half extraction window was 10, then rectified_2d_spectrum[trace_id]['val']
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
            self.logger.warning('Extraction window was chosen to be >= the max extraction window defined in settings.py.'
                           ' Defaulting to the max extraction window.')
            return rectified_2d_spectrum
        trim = self.max_extraction_half_window - self.extraction_half_window
        for order_id in list(rectified_2d_spectrum.keys()):
            for data_type in list(rectified_2d_spectrum[order_id].keys()):
                trimmed_rectified_spectrum[order_id][data_type] = rectified_2d_spectrum[order_id][data_type][trim:-trim]
        return trimmed_rectified_spectrum


class IVarExtract(BoxExtract):
    """
    Extraction with each pixel weighted by its inverse variance.
    """
    def __init__(self, runtime_context=None):
        super(IVarExtract, self).__init__(runtime_context=runtime_context)
        self.extraction_half_window = runtime_context.sne_extraction_half_window
        self.max_extraction_half_window = runtime_context.max_extraction_half_window
        self.table_name = runtime_context.ivar_spectrum_name

    def _weights(self, order_rect_spectrum, order_rect_ivar):
        normalization = self.extract_order(order_rect_ivar)
        return np.divide(order_rect_ivar, normalization, out=np.zeros_like(order_rect_ivar),
                         where=~np.isclose(normalization, 0))


class BlazeCorrectedExtract(IVarExtract):
    """
    Same as IVarExtract, meant to run after dividing the 2d spectrum by the blaze AND re-rectifying the spectrum.
    """
    def __init__(self, runtime_context=None):
        super(BlazeCorrectedExtract, self).__init__(runtime_context=runtime_context)
        self.extraction_half_window = runtime_context.sne_extraction_half_window
        self.max_extraction_half_window = runtime_context.max_extraction_half_window
        self.table_name = runtime_context.blaze_corrected_spectrum_name


class RectifyTwodSpectrum(Stage):
    def __init__(self, runtime_context=None):
        super(RectifyTwodSpectrum, self).__init__(runtime_context=runtime_context)
        self.max_extraction_half_window = runtime_context.max_extraction_half_window

    def do_stage(self, image):
        self.logger.info('Rectifying the 2d spectrum')
        if image.trace is None:
            self.logger.error('Image has empty trace attribute. Aborting extraction.')
            image.is_bad = True
            image.rectified_2d_spectrum = {}
            image.rectified_ivar = {}
            return image
        rectified_2d_spectrum = extract_utils.rectify_orders(image.data, image.trace,
                                                             half_window=self.max_extraction_half_window)
        image.rectified_2d_spectrum = rectified_2d_spectrum
        image.ivar = np.full_like(image.data, np.nan, dtype=np.double) if image.ivar is None else image.ivar
        image.rectified_ivar = extract_utils.rectify_orders(image.ivar, image.trace,
                                                            half_window=self.max_extraction_half_window,
                                                            nullify_mapped_values=False)
        return image


def safe_pow(arr, power):
    if arr is None:
        return arr
    return np.power(arr, power)
