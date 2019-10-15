import numpy as np


def median_subtract_channels_y(data, num_channels):
    """
    :param data: array_like
                 Input array. Must be 2D
    :param num_channels: The number of readout channels along axis=0 (rows).
    :return: ndarray
             Input array with each horizontal slice of size data.shape[0]/num_channels
             subtracted by its median.

    Examples
    --------
    >>> import numpy as np
    >>> a = (np.arange(3) * np.ones((3, 3))).T
    >>> print(a)
    array([[0.,  0.,  0.],
           [1.,  1.,  1.],
           [2.,  2.,  2.]])
    >>> print(median_subtract_channels_y(a, 3))
    array([[0.,  0.,  0.],
           [0.,  0.,  0.],
           [0.,  0.,  0.]])
    """
    reshaped = data.reshape(num_channels, data.shape[1], -1)
    medians = np.array([np.median(readout_channel) for readout_channel in reshaped])
    return np.subtract(reshaped, medians.reshape(num_channels, 1, 1)).reshape(data.shape)
