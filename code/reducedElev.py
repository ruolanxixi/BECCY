# Description: Generate reduced elevation topography
# Author: Ruolan Xiang, November 2021

# Load modules
import xarray as xr
from netCDF4 import Dataset
from pyproj import CRS, Transformer
from scipy.ndimage.interpolation import rotate
import matplotlib.pyplot as plt
import numpy as np
import math
from mpl_toolkits.axisartist.axislines import SubplotZero
import cartopy.crs as ccrs
import cartopy.feature as cfeature

###############################################################################
# Settings
###############################################################################

# Get files
path = "/Users/rxiang/Desktop/topo/"
src_file = 'GLOBE_H10.nc'
tar_file = 'GLOBE_H10_reduced.nc'
# src_file = 'GLOBE_H_filt_tukey_0.75_3.0_it4.nc'
# tar_file = 'GLOBE_H_filt_tukey_0.75_3.0_it4_reduced1.nc'
# src_file = 'GLOBE_H_filt_lanczos_window.nc'
# tar_file = 'GLOBE_H_filt_lanczos_window_reduced1.nc'

# open file
ds = xr.open_dataset(path + src_file)
HSURF = ds["altitude"].values
lat = ds["lat"].values
lon = ds["lon"].values
ds.close()

HSURF_reduced = HSURF

lon_0, lat_0 = lon[612], lat[2012]  # lon[1000], lat[2580]
crs_wgs84 = CRS.from_epsg(4326)
crs_aeqd = CRS.from_proj4("+proj=aeqd +lat_0=" + str(lat_0) + " +lon_0="
                          + str(lon_0) + " +datum=WGS84 +units=m")
transformer = Transformer.from_crs(crs_wgs84, crs_aeqd, always_xy=True)
lon_2d, lat_2d = np.meshgrid(lon, lat)
x, y = transformer.transform(lon_2d, lat_2d)
i = np.arrange(0, len(x[:, 0]))
j = np.arrange(len(x[0, :]), 0)


def compute_dis(rlon, rlat):
    distance = np.sqrt(rlon ** 2 + rlat ** 2) / 1000.0
    return distance


azi = np.empty(np.shape(x))
azi[i, j] = math.atan2(y[i, j], x[i, j])

azi_deg = np.rad2deg(azi)
azi45 = azi + np.pi / 4


def compute_frac_azi(rlon, rlat):
    if -np.pi / 2 <= azi45[rlon, rlat] <= np.pi / 2:
        fraction = np.cos(azi45[rlon, rlat]) ** 2
    else:
        fraction = 0
    return fraction


def compute_frac_dis(rlon, rlat):
    if compute_dis(rlon, rlat) <= 1800:
        fraction = np.sin(compute_dis(rlon, rlat) / 1800 * np.pi) ** 2
    else:
        fraction = 0
    return fraction


reduce = np.empty(np.shape(x))
weight_dis = np.empty(np.shape(x))

# compute reduce factor
reduce[i, j] = compute_frac_azi(i, j) * compute_frac_dis(x[i][j], y[i][j]) * 0.75

# compute redueced topography
HSURF_reduced[i, j] = max(500, HSURF[i, j] * (1 - reduce[i, j]))


tar_nc = Dataset(path + tar_file, 'a')
tar_nc['altitude'][:] = HSURF_reduced[:]
tar_nc.close()
