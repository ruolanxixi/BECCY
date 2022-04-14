# Description: Plot modified topography
#
# Authors: Ruolan Xiang, Christian R. Steger, IAC ETH Zurich

# Load modules
import xarray as xr
from pyproj import CRS, Transformer
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

mpl.style.use("classic")

###############################################################################
# Settings
###############################################################################

# Terrain reduction settings
lat_0, lon_0 = 33.23000, 95.10000  # reference location [degree]
rad_red = 1800.0  # reduction radius [km]
alpha_0, alpha_1 = -135.0, 45.0  # measured anti-clockwise from East [degree]
fac_amp = 0.75  # amplitude of terrain reduction
topo_min = 500.0  # minimal allowed elevation for terrain reduction

# Map projection
crs_wgs84 = CRS.from_epsg(4326)
crs_aeqd = CRS.from_proj4("+proj=aeqd +lat_0=" + str(lat_0) + " +lon_0="
                          + str(lon_0) + " +datum=WGS84 +units=m")

# DEM tiles
tiles_dem = ("MERIT_N60-N30_E060-E090.nc", "MERIT_N60-N30_E090-E120.nc",
             "MERIT_N30-N00_E060-E090.nc", "MERIT_N30-N00_E090-E120.nc")

# Paths
# path_dem = "/Users/kaktus/Documents/ETH/BECCY/myscripts/topo/"
path_dem = "/Users/csteger/Dropbox/IAC/Data/DEMs/MERIT/Tiles/"
path_out = "/Users/csteger/Desktop/Temp/"

###############################################################################
# Process MERIT data
###############################################################################

# fig = plt.figure()
# ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
#
# X = lat
# Y = lon
# Z = elev
#
# ax.plot(projection=ccrs.PlateCarree())
#
# # ax.set_title("Modern topography")
#
# ax.set_extent([90, 115, 15, 40], crs=ccrs.PlateCarree())
# ax.add_feature(cfeature.LAND)
# #ax.add_feature(cfeature.OCEAN, zorder=100, edgecolor='k')
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
