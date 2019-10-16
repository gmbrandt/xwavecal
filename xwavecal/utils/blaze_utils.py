import numpy as np


def normalize_orders(image_data, trace, half_window=10, n=100):
    """
    :param image_data: numpy.ndarray
    :param trace: xwavecal.utils.trace_utils.Trace
    :param half_window: int
                        the number of pixels above an below a diffraction order to try and normalize
    :param n: int
              the n brightest pixels of those near the trace will have their median computed, then
              all the pixels near the trace will be divided by that median.
    :return: normalization_factor: numpy.ndarray
             The array that normalizes each diffraction order of image_data to 1.
    """
    maxy = np.max(image_data.shape[0])
    x, y = np.meshgrid(np.arange(image_data.shape[1]), np.arange(image_data.shape[0]))
    normalization_factor = np.ones_like(image_data).astype(float) * np.inf
    for single_order in trace.data['centers']:
        roi = slice(int(max(0, np.min(single_order) - half_window)),
                    int(min(np.max(single_order) + half_window, maxy)))
        near_trace = np.where(np.isclose(y[roi] - single_order, 0, atol=half_window))
        median_of_n_brightest = np.median(np.sort(image_data[roi][near_trace])[::-1][:n])
        normalization_factor[roi][near_trace] = median_of_n_brightest
    return normalization_factor
