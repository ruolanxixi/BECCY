# Load modules
import numpy as np
import matplotlib as mpl


###############################################################################

def truncate_colormap(cmap_in, minval=0.0, maxval=1.0, n=100):
    """Truncate colormap.

    Truncate colormap to specific range between [0.0, 1.0].

    Parameters
    ----------
    cmap_in: colormap object
        Input colormap object
    minval: float
        Lower truncation value [0.0, 1.0] [-]
    maxval: float
        Upper truncation value [0.0, 1.0] [-]
    n: int
        Number of levels [-]

    Returns
    -------
    cmap_trun: colormap object
        Truncated colormap object

    Notes
    -----
    Author: Christian Steger (christian.steger@env.ethz.ch)"""

    # Check input arguments
    if minval >= maxval:
        raise TypeError("'minval' must be smaller than 'maxval'")

    cmap_trun = mpl.colors.LinearSegmentedColormap.from_list(
        "trunc({n},{a:.2f},{b:.2f})".format(
            n=cmap_in.name, a=minval, b=maxval),
        cmap_in(np.linspace(minval, maxval, n)))

    return cmap_trun

###############################################################################


def spat_agg_1d(data, agg_num, operation="sum"):
    """Aggregate 1-dimensional array.

    Parameters
    ----------
    data: ndarray
        Array (1-dimensional) with data [arbitrary]
    agg_num: int
        Aggregation number [-]
    operation: str
        Aggregation operation (sum or mean)

    Returns
    -------
    data_agg: ndarray
        Array (2-dimensional) with aggregated data [arbitrary]

    Notes
    -----
    Author: Christian Steger (christian.steger@env.ethz.ch)"""

    # Check input arguments
    if data.ndim != 1:
        raise TypeError("'data' must be 1-dimensional")
    if data.size % agg_num != 0:
        raise TypeError("argument 'agg_numb' does not fit with shape of data")
    if operation not in ("sum", "mean"):
        raise TypeError("unknown operation")
    if data.dtype == bool:
        print("Cast boolean array to type 'np.int32'")
        data = data.astype(np.int32)

    # Perform aggregation
    if operation == "sum":
        data_agg = np.sum(data.reshape(int(data.size / agg_num), agg_num),
                          axis=1)
    else:
        data_agg = np.mean(data.reshape(int(data.size / agg_num), agg_num),
                           axis=1)

    return data_agg


###############################################################################

def spat_agg_2d(data, agg_num_0, agg_num_1, operation="sum"):
    """Aggregate 2-dimensional array.

    Parameters
    ----------
    data: ndarray
        Array (2-dimensional) with data [arbitrary]
    agg_num_0: int
        Aggregation number along first dimension [-]
    agg_num_1: int
        Aggregation number along second dimension [-]
    operation: str
        Aggregation operation (sum or mean)

    Returns
    -------
    data_agg: ndarray
        Array (2-dimensional) with aggregated data [arbitrary]

    Notes
    -----
    Author: Christian Steger (christian.steger@env.ethz.ch)"""

    # Check input arguments
    if data.ndim != 2:
        raise TypeError("'data' must be 2-dimensional")
    if (data.shape[0] % agg_num_0 != 0) or (data.shape[1] % agg_num_1 != 0):
        raise TypeError("argument 'agg_numb' does not fit with shape of data")
    if operation not in ("sum", "mean"):
        raise TypeError("unknown operation")
    if data.dtype == bool:
        print("Cast boolean array to type 'np.int32'")
        data = data.astype(np.int32)

    # Perform aggregation
    y = np.arange(0, data.shape[0], agg_num_0)
    temp = np.add.reduceat(data, y, axis=0, dtype=data.dtype)
    x = np.arange(0, data.shape[1], agg_num_1)
    data_agg = np.add.reduceat(temp, x, axis=1, dtype=data.dtype)
    if operation == "mean":
        if np.issubdtype(data.dtype, np.integer):
            print("Cast integer array to type 'np.float32'")
            data_agg = data_agg.astype(np.float32)
        data_agg /= float(agg_num_0 * agg_num_1)

    return data_agg


###############################################################################

def gridcoord(x_cent, y_cent):
    """Compute edge coordinates.

    Compute edge coordinates from grid cell centre coordinates

    Parameters
    ----------
    x_cent: array_like
        Array (1-dimensional) with x-coordinates of grid centres [arbitrary]
    y_cent: array_like
        Array (1-dimensional) with y-coordinates of grid centres [arbitrary]

    Returns
    -------
    x_edge: array_like
        Array (1-dimensional) with x-coordinates of grid edges [arbitrary]
    y_edge: array_like
        Array (1-dimensional) with y-coordinates of grid edges [arbitrary]

    Notes
    -----
    Author: Christian Steger (christian.steger@env.ethz.ch)"""

    # Check input arguments
    if len(x_cent.shape) != 1 or len(y_cent.shape) != 1:
        raise TypeError("number of dimensions of input arrays is not 1")
    if (np.any(np.diff(np.sign(np.diff(x_cent))) != 0) or
            np.any(np.diff(np.sign(np.diff(y_cent))) != 0)):
        sys.exit("input arrays are not monotonically in- or decreasing")

    # Compute grid spacing if not provided
    dx = np.diff(x_cent).mean()
    dy = np.diff(y_cent).mean()

    # Compute grid coordinates
    x_edge = np.hstack((x_cent[0] - (dx / 2.),
                        x_cent[:-1] + np.diff(x_cent) / 2.,
                        x_cent[-1] + (dx / 2.))).astype(x_cent.dtype)
    y_edge = np.hstack((y_cent[0] - (dy / 2.),
                        y_cent[:-1] + np.diff(y_cent) / 2.,
                        y_cent[-1] + (dy / 2.))).astype(y_cent.dtype)

    return x_edge, y_edge
