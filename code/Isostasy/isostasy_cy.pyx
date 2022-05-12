#cython: boundscheck=False, wraparound=False, cdivision=True, language_level=3

# Load modules
import numpy as np
from libc.math cimport sqrt
from libc.math cimport sin, cos, acos
from libc.math cimport M_PI
from cython.parallel import prange
from scipy.special.cython_special cimport kei
from libc.stdlib cimport abs

# -----------------------------------------------------------------------------

# To do
# - allow computation of geodesic instead of great circle distance ->
#   (more accurate) -> use geographiclib with Cython bindings:
#   https://github.com/megaserg/geographiclib-cython-bindings

###############################################################################
# Compute deflection due to isostatic adjustment
###############################################################################

def deflection_xy(double[:] x, double[:] y, double[:, :] q,
                  double alpha, double d, imp="default"):
    """Compute vertical deflection for planar grid.

    Compute vertical deflection due to isostatic adjustment for
    planar grid.

    Parameters
    ----------
    x : ndarray of double
        Array (one-dimensional) with x-coordinates [m]
    y : ndarray of double
        Array (one-dimensional) with y-coordinates [m]
    q : ndarray of double
        Array (two-dimensional) with point loads [N]
    alpha : double
        Scalar with two-dimensional flexural parameter [m]
    d : double
        Scalar with flexural rigidity [m2 kg s-2]
    imp : str
        Implementation of algorithm (default or fast)

    Returns
    -------
    w : ndarray of double
        Array (two-dimensional) with vertical deflection [m]

    Sources
    -------
    - Wickert (2016): Open-source modular solutions for flexural isostasy:
      gFlex v1.0
    - Jha et al. (2017): Toolbox for Analysis of Flexural Isostasy (TAFI)—
      A MATLAB toolbox for modeling flexural deformation of the lithosphere

    Notes
    -----
    Author: Christian Steger (christian.steger@env.ethz.ch)"""

    cdef int len_0 = y.shape[0]
    cdef int len_1 = x.shape[0]
    cdef int i, j, k, m
    cdef float r
    cdef int ind
    cdef double[:, :] w = np.zeros((len_0, len_1), dtype=np.float64)
    cdef double[:, :] term = np.empty((len_0, len_1), dtype=np.float64)

    # Compute deflection due to point loads (default version)
    if imp == "default":
        for i in prange(len_0, nogil=True, schedule="static"):
            for j in range(len_1):
                for k in range(len_0):
                    for m in range(len_1):
                        r = sqrt((x[j] - x[m]) ** 2 + (y[i] - y[k]) ** 2)
                        w[i, j] += q[k, m] * alpha ** 2 / (2.0 * M_PI * d) \
                                   * kei(r / alpha)

    # Compute deflection due to point loads (fast version)
    else:
        for i in range(len_0):
            for k in prange(len_0, nogil=True, schedule="static"):
                for m in range(len_1):
                    r = sqrt((x[0] - x[m]) ** 2 + (y[i] - y[k]) ** 2)
                    term[k, m] = alpha ** 2 / (2.0 * M_PI * d) * kei(r / alpha)
            for j in prange(len_1, nogil=True, schedule="static"):
                for k in range(len_0):
                    for m in range(len_1):
                        ind = abs(j - m)
                        w[i, j] += q[k, m] * term[k, ind]

    return np.asarray(w)


# -----------------------------------------------------------------------------

def deflection_lonlat(double[:] lon, double[:] lat, double[:, :] q,
                      double alpha, double d, imp="default"):
    """Compute vertical deflection for spherical grid.

    Compute vertical deflection due to isostatic adjustment for
    spherical grid.

    Parameters
    ----------
    lon : ndarray of double
        Array (one-dimensional) with geographic longitude [rad]
    lat : ndarray of double
        Array (one-dimensional) with geographic latitude [rad]
    q : ndarray of double
        Array (two-dimensional) with point loads [N]
    alpha : double
        Scalar with two-dimensional flexural parameter [m]
    d : double
        Scalar with flexural rigidity [m2 kg s-2]
    imp : str
        Implementation of algorithm (default or fast)

    Returns
    -------
    w : ndarray of double
        Array (two-dimensional) with vertical deflection [m]

    Sources
    -------
    - Wickert (2016): Open-source modular solutions for flexural isostasy:
      gFlex v1.0
    - Jha et al. (2017): Toolbox for Analysis of Flexural Isostasy (TAFI)—
      A MATLAB toolbox for modeling flexural deformation of the lithosphere

    Notes
    -----
    Author: Christian Steger (christian.steger@env.ethz.ch)"""

    cdef int len_0 = lat.shape[0]
    cdef int len_1 = lon.shape[0]
    cdef int i, j, k, m
    cdef float r
    cdef int ind
    cdef double[:, :] w = np.zeros((len_0, len_1), dtype=np.float64)
    cdef double[:, :] term = np.empty((len_0, len_1), dtype=np.float64)

    # Compute deflection due to point loads (default version)
    if imp == "default":
        for i in prange(len_0, nogil=True, schedule="static"):
            for j in range(len_1):
                for k in range(len_0):
                    for m in range(len_1):
                        r = great_circ_dist(lat[i], lon[j], lat[k], lon[m])
                        w[i, j] += q[k, m] * alpha ** 2 / (2.0 * M_PI * d) \
                                   * kei(r / alpha)

    # Compute deflection due to point loads (fast version)
    else:
        for i in range(len_0):
            for k in prange(len_0, nogil=True, schedule="static"):
                for m in range(len_1):
                    r = great_circ_dist(lat[i], lon[0], lat[k], lon[m])
                    term[k, m] = alpha ** 2 / (2.0 * M_PI * d) * kei(r / alpha)
            for j in prange(len_1, nogil=True, schedule="static"):
                for k in range(len_0):
                    for m in range(len_1):
                        ind = abs(j - m)
                        w[i, j] += q[k, m] * term[k, ind]

    return np.asarray(w)


###############################################################################
# Auxiliary functions
###############################################################################

cdef inline double great_circ_dist(double lat_1, double lon_1,
                                   double lat_2, double lon_2) nogil:
    """Compute great-circle distance between two points on the sphere.

    Parameters
    ----------
    lat_1 : double
        Scalar with geographic latitude of first point [rad]
    lon_1 : double
        Scalar with geographic longitude of first point [rad]
    lat_2 : double
        Scalar with geographic latitude of second point [rad]
    lon_2 : double
        Scalar with geographic longitude of second point [rad]

    Returns
    -------
    dist : double
        Scalar with great-circle distance [m]

    Sources
    -------
    - https://en.wikipedia.org/wiki/Great-circle_distance"""

    cdef double rad_earth, term

    rad_earth = 6370997.0  # earth radius (PROJ default value) [m]
    term = sin(lat_1) * sin(lat_2) + cos(lat_1) * cos(lat_2) \
           * cos(lon_2 - lon_1)
    if term > 1.0:
        term = 1.0
    return rad_earth * acos(term)
