# Description: Find elevation threshold for permanent snow/ice coverage
#
# Source of data (beside MERIT DEM):
# - GlobCover 2009
#   - http://due.esrin.esa.int/files/Globcover2009_V2.3_Global_.zip
#
# Author: Christian R. Steger, IAC ETH Zurich

# Load modules
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.gridspec as gridspec
from matplotlib import colors
import cartopy.crs as ccrs
from osgeo import gdal
from cmcrameri import cm
import xarray as xr
import subprocess
from skimage.measure import label

mpl.style.use("classic")

###############################################################################
# Settings
###############################################################################

# Paths
path_globcov = "/Users/csteger/Dropbox/IAC/Data/Observations/Globcover2009/" \
               + "Globcover2009_V2/"
path_dem = "/Users/csteger/Dropbox/IAC/Data/DEMs/MERIT/Tiles/"
path_temp = "/Users/csteger/Dropbox/IAC/Temp/BECCY_glaciation/"
path_plot = "/Users/csteger/Dropbox/IAC/Plots/BECCY/"

###############################################################################
# Remap MERIT data to GlobCover 2009 grid
###############################################################################

# Define domain boundaries
dom_bound = (85.0, 105.0, 25.0, 35.0)  # lon (min/max), lat (min/max)

# Get relevant MERIT domain for remapping
tiles_dem = ("MERIT_N60-N30_E060-E090.nc", "MERIT_N60-N30_E090-E120.nc",
             "MERIT_N30-N00_E060-E090.nc", "MERIT_N30-N00_E090-E120.nc")
ds = xr.open_mfdataset([path_dem + i for i in tiles_dem])
ds = ds.sel(lon=slice(dom_bound[0] - 0.1, dom_bound[1] + 0.1),
            lat=slice(dom_bound[3] + 0.1, dom_bound[2] - 0.1))
ds.to_netcdf(path_temp + "MERIT_subdom.nc")

# Create remapping text file
xfirst, yfirst = dom_bound[0], dom_bound[3]
xinc, yinc = (1.0 / 360.0), -(1.0 / 360.0)
xsize = int((dom_bound[1] - dom_bound[0]) / np.abs(xinc))
ysize = int((dom_bound[3] - dom_bound[2]) / np.abs(yinc))
file = open(path_temp + "grid_out.txt", "w")
file.write("#\n")
file.write("# gridID 1\n")
file.write("#\n")
file.write("gridtype  = lonlat\n")
file.write("xsize     = " + str(xsize) + "\n")
file.write("ysize     = " + str(ysize) + "\n")
file.write("xfirst    = " + str(xfirst) + "\n")
file.write("xinc      = " + str(xinc) + "\n")
file.write("yfirst    = " + str(yfirst) + "\n")
file.write("yinc      = " + str(yinc) + "\n")
file.close()

# Remap conservatively with CDO
cmd = "cdo remapcon," + path_temp + "grid_out.txt"
sf = path_temp + "MERIT_subdom.nc"
tf = path_temp + "MERIT_GlobCover_2009.nc"
# subprocess.call(cmd + " " + sf + " " + tf, shell=True)

###############################################################################
# Load data
###############################################################################

# Define relevant domain size
dom_bound = (90.0, 105.0, 25.0, 35.0)  # lon (min/max), lat (min/max)

# Load GlobCover 2009 data (land classes and elevation)
ds = xr.open_dataset(tf)
ds = ds.sel(lon=slice(dom_bound[0], dom_bound[1]),
            lat=slice(dom_bound[3], dom_bound[2]))
topo_gc2009 = ds["Elevation"].values
ds.close()
ds = gdal.Open(path_globcov + "GLOBCOVER_L4_200901_200912_V2.3.tif")
lon_gc2009 = np.linspace(-180.0,
                         -180.0 + ds.RasterXSize * ds.GetGeoTransform()[1],
                         ds.RasterXSize + 1)[:-1]
lat_gc2009 = np.linspace(90.0,
                         90.0 + ds.RasterYSize * ds.GetGeoTransform()[5],
                         ds.RasterYSize + 1)[:-1]
slic = (slice(np.argmin(np.abs(lat_gc2009 - dom_bound[3])),
              np.argmin(np.abs(lat_gc2009 - dom_bound[2]))),
        slice(np.argmin(np.abs(lon_gc2009 - dom_bound[0])),
              np.argmin(np.abs(lon_gc2009 - dom_bound[1]))))
class_gc2009 = ds.GetRasterBand(1).ReadAsArray()[slic]
lon_gc2009, lat_gc2009 = lon_gc2009[slic[1]], lat_gc2009[slic[0]]

###############################################################################
# Analyse data
###############################################################################

# Compute mean elevation of snow/ice line for contiguous snow/ice pixel regions
labels, num = label((class_gc2009 == 220).astype(np.int32), connectivity=1,
                    return_num=True)
# 220: Permanent snow and ice
# 1-connectivity -> only check 4 directions
# labels = 0 -> background
num_grad = np.zeros(num, dtype=np.int32)
grad = np.zeros(num, dtype=np.float32)
print(" compute gradients in x-direction ".center(60, "-"))
for i in range(labels.shape[0]):
    for j in range(labels.shape[1] - 1):
        if np.diff(labels[i, j:(j + 2)])[0] != 0:
            ind = labels[i, j:(j + 2)].max() - 1
            num_grad[ind] += 1
            grad[ind] += topo_gc2009[i, j:(j + 2)].mean()
    if i % 500 == 0:
        print("Iteration i = " + str(i) + " completed")
print(" compute gradients in y-direction ".center(60, "-"))
for i in range(labels.shape[1]):
    for j in range(labels.shape[0] - 1):
        if np.diff(labels[j:(j + 2), i])[0] != 0:
            ind = labels[j:(j + 2), i].max() - 1
            num_grad[ind] += 1
            grad[ind] += topo_gc2009[j:(j + 2), i].mean()
    if i % 500 == 0:
        print("Iteration i = " + str(i) + " completed")
print("Minimum/maximum of 'num_grad': " + str(num_grad.min()) + ", "
      + str(num_grad.max()))
grad_mean = grad / num_grad

# Assign mean elevation of snow/ice line to regions
line_elev_mean = np.empty(labels.shape, dtype=np.float32)
line_elev_mean.fill(np.nan)
for i in range(labels.shape[0]):
    for j in range(labels.shape[1]):
        if labels[i, j] > 0:
            line_elev_mean[i, j] = grad_mean[labels[i, j] - 1]

###############################################################################
# Plot data
###############################################################################

# Colormaps
levels_topo = np.arange(0.0, 6500.0, 500.0)
cmap_topo = cm.grayC_r
norm_topo = mpl.colors.BoundaryNorm(levels_topo, ncolors=cmap_topo.N,
                                    extend="max")
levels = np.arange(3000.0, 6200.0, 200.)
cmap = cm.roma_r
norm = mpl.colors.BoundaryNorm(levels, ncolors=cmap.N, extend="both")

# Plot
font = {"family": "normal", "weight": "normal", "size": 10.5}
mpl.rc("font", **font)
map_ext = dom_bound
fig = plt.figure(figsize=(10.0, 13.5))
gs = gridspec.GridSpec(7, 1, left=0.1, bottom=0.1, right=0.9,
                       top=0.9, hspace=0.15, wspace=0.04,
                       height_ratios
                       =[0.9, 0.005, 0.030, 0.020, 0.030, 0.025, 0.4])
# -----------------------------------------------------------------------------
ax = plt.subplot(gs[0], projection=ccrs.PlateCarree())
plt.pcolormesh(lon_gc2009, lat_gc2009, topo_gc2009, cmap=cmap_topo,
               norm=norm_topo, shading="auto")
data_plot = np.ma.masked_where(np.isnan(line_elev_mean), line_elev_mean)
plt.pcolormesh(lon_gc2009, lat_gc2009, data_plot, cmap=cmap, norm=norm,
               shading="auto")
ax.set_aspect("auto")
plt.axis(map_ext)
gl = ax.gridlines(draw_labels=True, linestyle="-", linewidth=0.0)
gl.top_labels, gl.right_labels = False, False
plt.text(0.03, 0.96, "(a)", fontweight="bold",
         horizontalalignment="center", verticalalignment="center",
         transform=ax.transAxes)
# -----------------------------------------------------------------------------
ax = plt.subplot(gs[2])
mpl.colorbar.ColorbarBase(ax, cmap=cmap, norm=norm, orientation="horizontal")
plt.xlabel("Mean permanent snow/ice line [m]")
# -----------------------------------------------------------------------------
ax = plt.subplot(gs[4])
mpl.colorbar.ColorbarBase(ax, cmap=cmap_topo, norm=norm_topo,
                          orientation="horizontal")
plt.xlabel("Terrain elevation [m]")
# -----------------------------------------------------------------------------
ax = plt.subplot(gs[6])
data = line_elev_mean  # [mask_gc2009]
data = data[np.isfinite(data)]
l0 = plt.hist(data, bins=30, density=True, color="blue", alpha=0.7, zorder=2)
for i in (5.0, 50.0, 95.0):  # percentiles
    perc_value = np.percentile(data, i)
    plt.vlines(perc_value, ymin=0.0, ymax=1.0, color="black",
               lw=2.0, ls="-", zorder=3)
    txt = "~" + str(int(np.round(perc_value / 10.0, decimals=0) * 10)) + " m"
    plt.text(perc_value + 35.0, 0.0012, txt, fontweight="bold")
plt.axis([2500.0, 7000.0, 0.0, 0.0013])
plt.xlabel("Elevation [m]")
plt.ylabel("Density [-]")
plt.text(0.03, 0.94, "(b)", fontweight="bold",
         horizontalalignment="center", verticalalignment="center",
         transform=ax.transAxes)
# -----------------------------------------------------------------------------
fig.savefig(path_plot + "Snow_ice_threshold_elevation.png", dpi=300,
            bbox_inches="tight")
plt.close(fig)
