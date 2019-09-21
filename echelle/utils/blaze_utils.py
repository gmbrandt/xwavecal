import numpy as np


def normalize_orders(image_data, trace, minval=1, half_window=10, n=100):
    """
    :param image_data: numpy.ndarray
    :param trace: echelle.utils.trace_utils.Trace
    :param minval: float
                   minimum value allowed. E.g. 10 * read_noise.
    :param buffer: int
    :return: numpy.ndarray
             image_data where each order has had it's brightest pixels normalized to 1.
    """
    maxy = np.max(image_data.shape[0])
    x, y = np.meshgrid(np.arange(image_data.shape[1]), np.arange(image_data.shape[0]))
    normalization_factor = np.ones_like(image_data).astype(float)
    for single_order in trace.data['centers']:
        roi = slice(int(max(0, np.min(single_order) - half_window)),
                    int(min(np.max(single_order) + half_window, maxy)))
        near_trace = np.where(np.isclose(y[roi] - single_order, 0, atol=half_window))
        median_of_n_brightest = np.median(np.sort(image_data[roi][near_trace])[::-1][:n])
        normalization_factor[roi][near_trace] = median_of_n_brightest
    image_data[image_data < minval] = minval
    image_data = image_data / normalization_factor
    return image_data
