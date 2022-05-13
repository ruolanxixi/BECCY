# Description: Compute reduced topography part II (plot modified topography)
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
import matplotlib.ticker as mticker
import matplotlib.gridspec as gridspec
import sys
from cmcrameri import cm

mpl.style.use("classic")

# Load required functions
sys.path.append("/Users/csteger/Downloads/BECCY/code/")
from auxiliary import truncate_colormap, spat_agg_1d, spat_agg_2d

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
path_dem = "/Users/csteger/Dropbox/IAC/Data/DEMs/MERIT/Tiles/"
path_in_out = "/Users/csteger/Dropbox/IAC/Data/Model/BECCY/MERIT_tiles_red/"
path_plot = "/Users/csteger/Dropbox/IAC/Plots/BECCY/Topo_reduce/"

###############################################################################
# Process MERIT data
###############################################################################

# Aggregation number
agg_num = 10

# Load unmodified topography and aggregate
ds = xr.open_mfdataset([path_dem + i for i in tiles_dem])
ds = ds.sel(lon=slice(90.0, 110.0 - 0.0001), lat=slice(41.0 - 0.0001, 15.0))
topo = spat_agg_2d(ds["Elevation"].values, agg_num, agg_num, "mean")
lon = spat_agg_1d(ds["lon"].values, agg_num, "mean")
lat = spat_agg_1d(ds["lat"].values, agg_num, "mean")
ds.close()

# Load modified topography
ds = xr.open_mfdataset([path_in_out + i for i in tiles_dem])
ds = ds.sel(lon=slice(90.0, 110.0 - 0.0001), lat=slice(41.0 - 0.0001, 15.0))
topo_red = spat_agg_2d(ds["Elevation"].values, agg_num, agg_num, "mean")
fac_red = spat_agg_2d(ds["fac_red"].values, agg_num, agg_num, "mean")
ds.close()

# Plot reduction factor
levels = np.arange(0., 0.95, 0.05)
cmap = cm.lajolla
norm = mpl.colors.BoundaryNorm(levels, ncolors=cmap.N)
fig = plt.figure()
ax = plt.axes(projection=ccrs.PlateCarree())
plt.pcolormesh(lon, lat, fac_red, shading="auto", cmap=cmap, norm=norm)
ax.add_feature(cfeature.COASTLINE)
ax.add_feature(cfeature.BORDERS, linestyle=':')
ax.set_extent([90.1, 109.9, 15.1, 39.9], crs=ccrs.PlateCarree())
gl = ax.gridlines(draw_labels=True, linewidth=1, color="grey", alpha=0.5,
                  linestyle="--")
gl.xlocator = mticker.FixedLocator([90, 95, 100, 105, 110])
gl.top_labels = False
gl.right_labels = False
plt.colorbar()
plt.title("Topography reduction factor [-]", fontsize=12, fontweight="bold")
fig.savefig(path_plot + "reduction_factor.png", dpi=300, bbox_inches="tight")
plt.close(fig)

# Plot topography (unmodified, reduced and difference)
map_ext = [90.1, 109.9, 15.1, 39.9]
fig = plt.figure(figsize=(11.0, 13.0))
gs = gridspec.GridSpec(3, 3, left=0.1, bottom=0.1, right=0.9,
                       top=0.9, hspace=0.08, wspace=0.04,
                       width_ratios=[1, 1, 0.043], height_ratios=[1, 1, 0.043])
# -----------------------------------------------------------------------------
levels = np.arange(0., 6500.0, 500.0)
ticks = np.arange(0., 6500.0, 500.0)
cmap = truncate_colormap(cm.bukavu, 0.55, 1.0)
norm = mpl.colors.BoundaryNorm(levels, ncolors=cmap.N, extend="max")
ax = plt.subplot(gs[0, 0], projection=ccrs.PlateCarree())
plt.pcolormesh(lon, lat, topo, shading="auto", cmap=cmap, norm=norm)
ax.add_feature(cfeature.COASTLINE)
ax.add_feature(cfeature.BORDERS, linestyle=':')
ax.add_feature(cfeature.RIVERS, alpha=0.5)
ax.set_aspect("auto")
plt.axis(map_ext)
plt.title("Modern topography [m]", fontsize=12, fontweight="bold")
# -----------------------------------------------------------------------------
ax = plt.subplot(gs[0, 1], projection=ccrs.PlateCarree())
plt.pcolormesh(lon, lat, topo_red, shading="auto", cmap=cmap, norm=norm)
ax.add_feature(cfeature.COASTLINE)
ax.add_feature(cfeature.BORDERS, linestyle=':')
ax.set_aspect("auto")
plt.axis(map_ext)
plt.title("Reduced topography [m]", fontsize=12, fontweight="bold")
# -----------------------------------------------------------------------------
ax = plt.subplot(gs[0, 2])
cb = mpl.colorbar.ColorbarBase(ax, cmap=cmap, norm=norm, ticks=ticks,
                               orientation="vertical")
# cb.ax.tick_params(labelsize=10)
# plt.title("[m]", fontsize=10, y=1.018, loc="left")
# -----------------------------------------------------------------------------
levels = np.arange(0., 2750.0, 250.0)
ticks = np.arange(0., 2750.0, 500.0)
cmap = cm.lajolla
norm = mpl.colors.BoundaryNorm(levels, ncolors=cmap.N, extend="max")
ax = plt.subplot(gs[1, 0], projection=ccrs.PlateCarree())
plt.pcolormesh(lon, lat, topo - topo_red, shading="auto", cmap=cmap, norm=norm)
ax.add_feature(cfeature.COASTLINE)
ax.add_feature(cfeature.BORDERS, linestyle=':')
ax.set_aspect("auto")
plt.axis(map_ext)
plt.title("Reduction in topography [m]", fontsize=12, fontweight="bold")
# -----------------------------------------------------------------------------
ax = plt.subplot(gs[2, 0])
cb = mpl.colorbar.ColorbarBase(ax, cmap=cmap, norm=norm, ticks=ticks,
                               orientation="horizontal")
# -----------------------------------------------------------------------------
fig.savefig(path_plot + "reduced_topography.png", dpi=300, bbox_inches="tight")
plt.close(fig)
