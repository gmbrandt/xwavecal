import numpy as np

from xwavecal.utils.correlate import correlate2d


def test_correlate2d():
    arr = np.ones((5, 5))
    arr[1:4, 1:4] = 2
    sig = correlate2d(arr, 2 * np.ones((3, 3)), max_lag=1)
    assert sig[2, 2] == np.max(sig)