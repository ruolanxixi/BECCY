# Description: Generate reduced elevation topography

# Load modules
import xarray as xr
from netCDF4 import Dataset
from pyproj import CRS, Transformer
import numpy as np
import math
from pathlib import Path

###############################################################################
# Settings
###############################################################################

# Get files
path = "/Users/kaktus/Documents/ETH/BECCY/myscripts/topo/"
#files = [path + i for i in ("MERIT_N60-N30_E060-E090.nc")]
files = [path + i for i in ("MERIT_N60-N30_E060-E090.nc", "MERIT_N60-N30_E090-E120.nc",
                            "MERIT_N30-N00_E060-E090.nc", "MERIT_N30-N00_E090-E120.nc")]
ds = xr.open_mfdataset(files)
elev = ds["Elevation"].values
lat = ds["lat"].values
lon = ds["lon"].values
ds.close()

elev_reduced = elev
ny, nx = elev_reduced.shape

lon_0, lat_0 = 95.10000, 33.23000  # lon[612], lat[2012]  # lon = 95.10 lat = 33.23
crs_wgs84 = CRS.from_epsg(4326)
crs_aeqd = CRS.from_proj4("+proj=aeqd +lat_0=" + str(lat_0) + " +lon_0="
                          + str(lon_0) + " +datum=WGS84 +units=m")
transformer = Transformer.from_crs(crs_wgs84, crs_aeqd, always_xy=True)
lon_2d, lat_2d = np.meshgrid(lon, lat)
x, y = transformer.transform(lon_2d, lat_2d)

azi = np.arctan2(y, x)
azi_deg = np.rad2deg(azi)
azi45 = azi + np.pi / 4

frac_azi = np.cos(azi45) ** 2
frac_azi_mask = np.ma.masked_outside(frac_azi, -np.pi / 2, np.pi / 2)

distance = np.sqrt(x ** 2 + y ** 2) / 1000.0
frac_dis = np.sin(distance / 1800 * np.pi) ** 2
frac_dis_mask = np.ma.masked_greater(frac_dis, 1800)

mask = np.ma.getmask(frac_azi_mask*frac_dis_mask)  # no change = True
# compute reduce factor
reduce_factor = frac_azi * frac_dis * 0.75
elev_reduced = elev * (1 - reduce_factor*~mask)

# Set where the elevation was higher than 500 m and has been changed lower than 500 m to 500 m
elev_diff = elev - elev_reduced
mask = np.logical_and(elev_diff > 0.1, elev > 500, elev_reduced < 500)
elev_reduced = np.ma.masked_array(elev_reduced, mask=mask)

elev_reduced.filled(fill_value=500)

for p in Path(f'{path}').rglob('Reduced*'):
    print(p.name)
    ds = Dataset(path + p.name, 'a')
    if p.name == 'Reduced_MERIT_N60-N30_E060-E090.nc':
        ds['Elevation'][:] = elev_reduced[0:ny / 2, 0:nx / 2]
    elif p.name == 'Reduced_MERIT_N60-N30_E090-E180.nc':
        ds['Elevation'][:] = elev_reduced[0:ny / 2, nx / 2:]
    elif p.name == 'Reduced_MERIT_N30-N00_E060-E90.nc':
        ds['Elevation'][:] = elev_reduced[ny / 2:, 0:nx / 2]
    else:
        ds['Elevation'][:] = elev_reduced[ny / 2:, nx / 2:]
    ds.close()
