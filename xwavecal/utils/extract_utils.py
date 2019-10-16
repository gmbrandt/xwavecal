import numpy as np
from scipy.ndimage import map_coordinates


def rectify_orders(image_data, trace, half_window=10, debug=False, nullify_mapped_values=True):
    """
    :param image_data: xwavecal.images.Image.data
    :param trace: Trace object
    :param half_window: half of the size of the extraction window about the center of the traces. E.g. 10 means 21 pixel
                        window will be extracted.
    :param debug: boolean for returning an internal copy of the image which is used in extract.
    :param nullify_mapped_values: bool.
                                  True if the flux from a pixel can only be used in one order's spectrum.
                                  False if a pixel's flux can be double (or triple etc..) counted.
    :return rectified_orders: Dictionary
    Dictionary where the keys are the trace id's from trace.get_id(),
    where rectified_2d_spectrum[trace_id]['val'] is a 2d float array (flux for the trace_id order).
    If half extraction window was 10, then rectified_2d_spectrum[trace_id]['val']
    is 21 rows by 4096 columns (for a 4096 pixel wide image). One would sum this 2d
    spectrum along-columns to get a box extracted spectrum.
    """
    rectified_orders = {}
    num_orders = trace.num_traces_found()
    image_copy = np.copy(image_data)
    x_coordinates, y_coordinates = np.meshgrid(np.arange(image_data.shape[1]), np.arange(image_data.shape[0]))
    image_coordinates = {'x': x_coordinates, 'y': y_coordinates}
    for index in np.arange(num_orders):
        rectified_orders[trace.get_id(index)], image_copy = rectify_order(image_copy, image_coordinates,
                                                                          trace.get_centers(index),
                                                                          half_window=half_window,
                                                                          nullify_mapped_values=nullify_mapped_values)
    if debug:
        return rectified_orders, image_copy
    else:
        return rectified_orders


def rectify_order(image_data, image_coordinates, single_order_centers, half_window=10, nullify_mapped_values=True):
    """
    :param image_data: xwavecal.images.Image.data
    :param image_coordinates: dictionary with x and y coordinates for each pixel.
    :param single_order_centers: the y centers of the center of the rectification window
    :param half_window: the half width of the window about which one is rectifying the trace/order
    :param nullify_mapped_values: True if the flux from a pixel can only be used in one order's spectrum.
                                  False if a pixel's flux can be double (or triple etc..) counted.
    :return: An array of 2*half_window + 1 rows by width(image_data) about the order. The center of the order
             is right at the index of half_window, e.g. rectified_order[half_window] gives the flux down the
             center of the order.
    """
    spectrum_shape = (2*half_window + 1, image_data.shape[1])
    rectified_order = {'val': np.zeros(spectrum_shape),
                       'y': np.zeros(spectrum_shape),
                       'x': np.zeros(spectrum_shape)}
    x_coords = np.arange(image_data.shape[1])
    for offset, row in zip(np.arange(-half_window, half_window + 1), np.arange(spectrum_shape[0])):
        mapped_y_values = map_coordinates(image_coordinates['y'], [single_order_centers + offset, x_coords],
                                          order=0, mode='constant', cval=0, prefilter=False)
        rectified_order['val'][row] = image_data[(mapped_y_values, x_coords)]
        rectified_order['y'][row] = mapped_y_values
        rectified_order['x'][row] = x_coords
        if nullify_mapped_values:
            image_data[(mapped_y_values, x_coords)] = 0
    return rectified_order, image_data
