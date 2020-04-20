import numpy as np

from scipy import signal, optimize, interpolate
import warnings


def normalize_by_brightest(two_d_spectrum, n=10):
    """
    :param two_d_spectrum: ndarray 2 dimensional, where two_d_spectrum[i] is the ith one d spectrum.
    :return: two_d_spectrum where each row has been divided by norm,
             where norm is the median of the flux from the brightest n pixels of that row.
    NOTE: norm can be no smaller than 1.
    """
    temp = np.copy(two_d_spectrum)
    n = min(two_d_spectrum.shape[1], n)
    for i in range(two_d_spectrum.shape[0]):
        brightest_pixels = np.sort(temp[i])[::-1][:n]
        temp[i] /= max(np.median(brightest_pixels), 1)
    return temp


def n_largest_elements(arr, n):
    """
    :param arr: ndarray
    :param n: int
    :return: returns the indices for the n largest elements of arr.
    """
    return arr.argsort()[-n:][::-1]


def minmax(arr):
    return np.min(arr), np.max(arr)


def overlap_region(a, b):
    """
    :param a: ndarray or list
    :param b: ndarray or list
    :return: (left, right):
             left: int. smallest value of the intersection of a and b if they were continuous arrays.
             right: int. largest value of the intersection of a and b if they were continuous arrays.
             returns None if there is no overlap
    Examples:
            overlap_region([1, 2.1], [1.9, 3]) returns (1.9, 2.1)
            overlap_region([1, 2.1], [2.2, 3]) returns None
    """
    left, right = np.max([np.min(a), np.min(b)]), np.min([np.max(a), np.max(b)])
    return (left, right) if right > left else None


def find_nearest(array_a, array_b):
    """
    # TODO consider replacing with np.searchsorted.
    :param array_a: ndarray
    :param array_b: ndarray. Must be Sorted.
    :return: a list of elements of B which are closest to each element in A, searched for in order of A.
    e.g. if A = [1, 2, 2, 2] and B = [0.99, 1.99, 4]
    then we would return [0.99, 1.99, 1.99, 1.99]
    Believe it or not, this is 1000 times faster (for a list of length 10000) then a python implementation
    of bisection-searching array_b for every element in array_a.
    """
    find_func = interpolate.interp1d(array_b, array_b, kind='nearest', bounds_error=False,
                                     fill_value=(np.min(array_b), np.max(array_b)))
    return find_func(array_a)


def find_peaks(y, x, yerr, height=0, distance=1, prominence=None, window=6, peak_width=2):
    """
    :param y: ndarray. float. Amplitude as a function of x
    :param x: ndarray. dtype=int. x coordinates for y.
    :param yerr: ndarray. float. 1 sigma errors for y.
    :param height: float. peaks with amplitudes below height will be excluded. scalar or ndarray of length y.
    :param distance: int. minimum number of grid points (distance) between two identified peaks.
    :param prominence: float. The vertical distance between the peak and its lowest contour line.
                       Recommend prominence = 0.5 * np.abs(y)
    :param peak_width: estimate for the full width at half max of the peaks.
    :param window: int. A guassian will be fit around the peak but within a window of width.
    :return: ndarray, ndarray, ndarray
             peaks, peak_errs, peak_indices
             peaks is the peak centroid in the given x coordinates
             peak_errs is the 1 sigma error on each centroid in peaks.
             peak_indices the indices of the closest x coordinate to each location in peak_locations.

    Note on prominence:
        Strategy to compute a peak’s prominence:
        - Extend a horizontal line from the current peak to the left and right until the line either
          reaches the window border (see wlen) or intersects the signal again at the slope of a
          higher peak. An intersection with a peak of the same height is ignored.
        - On each side find the minimal signal value within the interval defined above.
          These points are the peak’s bases.
        - The higher one of the two bases marks the peak’s lowest contour line. The prominence can
          then be calculated as the vertical difference between the peaks height itself and
          its lowest contour line.

    Note: might do better to clip y for all points below height, instead of using signal.find_peaks height argument.
    Note: If the gaussian fit a peak position which differs by more than 3 from the initial guess,
          then we return the initial guess.
    """
    peak_indices = signal.find_peaks(y, distance=distance, prominence=prominence, height=height)[0]  # identify peaks.
    peaks = fit_peaks(x, y, yerr, window // 2, peak_indices, std=peak_width)  # fit a Gaussian to each peak.
    # remove fits that were wildly different from init guess, or have non finite errors:
    bad_fits = np.where(np.logical_or(~np.isclose(peaks[:, 0], x[peak_indices], atol=5), ~np.isfinite(peaks[:, 1])))
    peaks[bad_fits, 0], peaks[bad_fits, 1] = x[peak_indices[bad_fits]], 1
    return peaks[:, 0].flatten(), peaks[:, 1].flatten(), peak_indices


def fit_peaks(x, y, yerr, half_width, index_guess, std):
    peaks = []
    amplitude_guess = y[index_guess]
    for i, ampl in zip(index_guess, amplitude_guess):
        sl = slice(i - half_width, i + half_width + 1)
        center, err = x[i], 1
        if len(x[sl]) >= 3:  # num params in the model
            try:
                params, pcov = optimize.curve_fit(gaussian, x[sl], y[sl], [ampl, center, std],
                                                  sigma=yerr[sl], absolute_sigma=True)
                center, err = params[1], np.sqrt(np.diag(pcov))[1]
            except RuntimeError as e:
                warnings.warn(str(e))

        peaks.append([center, err])
    return np.array(peaks)


def gaussian(x, ampl, mu, std):
    return ampl * np.exp(-np.power(x - mu, 2.) / (2 * np.power(std, 2.)))


def brute_local_min(fun, x0, args, rrange=(0.1, 10), filtw=501, step=10, finish='Nelder-Mead', **kwargs):
    """
    :param fun: function you wish to find the deepest local minimum for. fun(x, args) must be able to
                accept an ndarray x and return an ndarray of fun evaluated at all the points in x.
    :param x0: initial guess
    :param args: extra args to pass to fun
    :param rrange: (low, high). The brute force minimization over the grid low*x0 and high*x0 with spacing step.
    :param filtw: the filter window for the median filter of fun() evaluating over the grid of points
    :param step: interval size in the brute force grid
    :param finish: Finish method to pass to scipy.optimize.minimize to refine the brute force minimum.
                   Note: use a method like Newton-Raphson that prefers local minima.
    :param kwargs: dict. Placeholder which allows scipy.optimize to pass standard extraneous arguments.
    :return: OptimizeResult with OptimizeResult.x = np.array([minimum]).
    Note: this is intended to be used with optimize.minimize with optimize.minimize(..., method=brute_local_min)
    """
    eval_points = np.arange(x0 * rrange[0], x0 * rrange[1], step=step)
    fun_vals = np.array(fun(eval_points, *args))
    fun_vals /= signal.medfilt(fun_vals, filtw)
    minimum_pt = eval_points[fun_vals == np.min(fun_vals[filtw:-filtw])][:1]
    if finish is not None:
        minimum_pt = optimize.minimize(fun, x0=minimum_pt, args=args, method=finish).x
    return optimize.OptimizeResult(x=minimum_pt)
