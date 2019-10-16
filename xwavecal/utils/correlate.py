from scipy import signal
import numpy as np


def correlate2d(array, template, max_lag=100):
    """
    wrapper for signal.correlate2d in which you can specify a max offset for the correlation.
    works by padding the input array with zeros and enforcing 'valid' mode in signal.correlate 2d so
    that the template cannot slide off the array.
    Both dimensions of array must be larger than template.
    -----
    :param array: 2d ndarray to correlate with the template.
    :param template: 2d ndarray.
    :param max_lag: The maximum correlation offset in either dimension.
    :return 2d cross correlation signal as a function of offset.
    """
    top_pad = int(template.shape[0] / 2)
    temp_array = np.pad(array, pad_width=((top_pad, top_pad), (max_lag, max_lag)),
                        mode='constant', constant_values=0)
    return signal.correlate2d(temp_array, template, mode='valid')
