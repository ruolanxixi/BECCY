# Description: Generate reduced elevation topography

# Load modules
import xarray as xr
from netCDF4 import Dataset
from pyproj import CRS, Transformer
import matplotlib.pyplot as plt
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfeature

###############################################################################
# Settings
###############################################################################

# Get files
path = "/Users/kaktus/Documents/ETH/BECCY/myscripts/topo/"
# for files in ("MERIT_N60-N30_E060-E090.nc", "MERIT_N60-N30_E090-E120.nc",
#              "MERIT_N30-N00_E060-E090.nc", "MERIT_N30-N00_E090-E120.nc"):
files = "MERIT_N30-N00_E090-E120.nc"
# files = "GLOBE_H10.nc"
ds = xr.open_dataset(path + files)
# elev = ds["altitude"].values
elev = ds["Elevation"].values
lat = ds["lat"].values
lon = ds["lon"].values
ds.close()
del ds

elev_reduced = elev
ny, nx = elev_reduced.shape

lon_0, lat_0 = 95.10000, 33.23000  # lon[612], lat[2012]  # lon = 95.10 lat = 33.23
crs_wgs84 = CRS.from_epsg(4326)
crs_aeqd = CRS.from_proj4("+proj=aeqd +lat_0=" + str(lat_0) + " +lon_0="
                          + str(lon_0) + " +datum=WGS84 +units=m")
transformer = Transformer.from_crs(crs_wgs84, crs_aeqd, always_xy=True)
lon_2d, lat_2d = np.meshgrid(lon, lat)
x, y = transformer.transform(lon_2d, lat_2d)
del lon_2d, lat_2d

azi = np.arctan2(y, x)  # degree in rad
# azi_deg = np.rad2deg(azi)
azi45 = azi + np.pi / 4
del azi

frac_azi = np.cos(azi45) ** 2
frac_azi_mask = np.ma.masked_inside(azi45, -np.pi / 2, np.pi / 2)
del azi45

distance = np.sqrt(x ** 2 + y ** 2) / 1000.0
frac_dis = np.sin(distance / 1800 * np.pi) ** 2
frac_dis_mask = np.ma.masked_less_equal(distance, 1800)
del x, y, distance

mask = np.ma.getmask(frac_azi_mask) * np.ma.getmask(frac_dis_mask)  # change = True

# compute reduce factor
reduce_factor = frac_azi * frac_dis * 0.75
elev_reduced = elev * (1 - reduce_factor * mask)
del reduce_factor, mask, frac_azi_mask, frac_dis_mask, frac_azi, frac_dis

# Set where the elevation was higher than 500 m and has been changed lower than 500 m to 500 m
elev_diff = elev - elev_reduced
mask2 = np.ma.getmask(np.ma.masked_greater(elev_diff, 1)) * \
       np.ma.getmask(np.ma.masked_greater(elev, 500)) * \
       np.ma.getmask(np.ma.masked_less(elev_reduced, 500))

elev_reduced_modi = np.ma.masked_array(elev_reduced, mask=mask2)
del elev, elev_diff, mask, elev_reduced

elev_reduced_modi = elev_reduced_modi.filled(fill_value=500)

ds = Dataset(path + f"Reduced_{files}", 'a')
ds['Elevation'][:] = elev_reduced_modi
# ds['altitude'][:] = elev_reduced_modi
ds.close()
del ds
del elev_reduced_modi

###############################################################################
# Plot
###############################################################################
# fig = plt.figure()
# ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
#
# X = lat
# Y = lon
# Z = reduce_factor
#
# ax.plot(projection=ccrs.PlateCarree())
#
# # ax.set_title("Modern topography")
#
# ax.set_extent([90, 115, 15, 40], crs=ccrs.PlateCarree())
# ax.add_feature(cfeature.LAND)
# ax.add_feature(cfeature.OCEAN, zorder=100, edgecolor='k')
# ax.add_feature(cfeature.COASTLINE)
# ax.add_feature(cfeature.BORDERS, linestyle=':')
# ax.add_feature(cfeature.LAKES, alpha=0.5)
# ax.add_feature(cfeature.RIVERS)
#
# cs = ax.contourf(Y, X, Z, transform=ccrs.PlateCarree(), levels=20, cmap='YlOrBr')
#
# ax.gridlines(draw_labels=True, linewidth=1, color='grey', alpha=0.5, linestyle='--')
#
# cb = fig.colorbar(cs, orientation='horizontal', shrink=0.6, ax=ax, pad=0.1)
# # cb.set_ticks([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
# fig.show()

