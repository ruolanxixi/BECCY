# Description: Find elevation threshold for permanent snow/ice coverage
#
# Source of data (beside MERIT DEM):
# - GlobCover 2009
#   - http://due.esrin.esa.int/files/Globcover2009_V2.3_Global_.zip
# - GAMDAM glacier inventory
#   - Central Asia:	   https://doi.org/10.1594/PANGAEA.891415
#   - South Asia East: https://doi.org/10.1594/PANGAEA.891417
#
# Authors: Christian R. Steger, IAC ETH Zurich

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
import math
from skimage.measure import label
from pyproj import CRS, Transformer
# from scipy.spatial import cKDTree

mpl.style.use("classic")

# Paths to folders
root_IAC = os.getenv("HOME") + "/Dropbox/IAC/"
path_globcov = root_IAC + "Data/Observations/Globcover2009/Globcover2009_V2/"
path_GAMDAM = root_IAC + "Data/Shapefiles/GAMDAM/Area_altitude_distribution/"
path_dem = "/Users/csteger/Dropbox/IAC/Data/DEMs/MERIT/Tiles/"
path_temp = root_IAC + "Temp/BECCY_glaciation/"
path_out = "/Users/csteger/Desktop/"

###############################################################################
# Settings
###############################################################################

# Reference location
lat_0, lon_0 = 33.23000, 95.10000  # reference location [degree]

# Map projection
crs_wgs84 = CRS.from_epsg(4326)
crs_aeqd = CRS.from_proj4("+proj=aeqd +lat_0=" + str(lat_0) + " +lon_0="
                          + str(lon_0) + " +datum=WGS84 +units=m")

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

# -----------------------------------------------------------------------------

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

# -----------------------------------------------------------------------------

# Load glacier data from GAMDAM inventory
usecols = (2, 3, 5)  # longitude/latitude [deg], elevation [m a.s.l.] (median)
data_GAMDAM = np.vstack((
    np.genfromtxt(path_GAMDAM + "GAMDAM_area-altitude_Central_Asia.tab",
                  skip_header=150, usecols=usecols),
    np.genfromtxt(path_GAMDAM + "GAMDAM_area-altitude_South_Asia_East.tab",
                  skip_header=182, usecols=usecols)))

# Select relevant data for domain
mask_sel = np.zeros(data_GAMDAM.shape[0], dtype=bool)
for i in range(data_GAMDAM.shape[0]):
    if (dom_bound[0] <= data_GAMDAM[i, 0] <= dom_bound[1]) \
            and (dom_bound[2] <= data_GAMDAM[i, 1] <= dom_bound[3]):
        mask_sel[i] = True
data_GAMDAM = data_GAMDAM[mask_sel, :]
print("Number of glaciers: " + str(data_GAMDAM.shape[0]))

###############################################################################
# Analyse and plot data
###############################################################################

# Compute mean elevation of snow/ice line for contiguous snow/ice pixel regions
labels, num = label((class_gc2009 == 220).astype(np.int32), connectivity=1,
                    return_num=True)
# 220: Permanent snow and ice
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

# # Masks for south-eastern region
# # ---------------------------------------------------------------------------
# lon, lat = np.meshgrid(lon_gc2009, lat_gc2009)
# trans_ellps2cart = Transformer.from_crs(crs_wgs84, crs_aeqd, always_xy=True)
# x, y = trans_ellps2cart.transform(lon, lat)
# azim = np.rad2deg(np.arctan2(y, x))  # [degree]
# mask_gc2009 = (azim >= -135.0) & (azim <= 45.0)
# # ---------------------------------------------------------------------------
# trans_ellps2cart = Transformer.from_crs(crs_wgs84, crs_aeqd, always_xy=True)
# x, y = trans_ellps2cart.transform(data_GAMDAM[:, 0], data_GAMDAM[:, 1])
# azim = np.rad2deg(np.arctan2(y, x))  # [degree]
# mask_GAMDAM = (azim >= -135.0) & (azim <= 45.0)
# # ---------------------------------------------------------------------------
# trans_cart2ellps = Transformer.from_crs(crs_aeqd, crs_wgs84, always_xy=True)
# num = 500
# rad = np.concatenate([np.linspace(1000.0, 0.0, num),
#                       np.linspace(0.0, 1000.0, num)])  # [km]
# azim = np.concatenate([np.repeat(-135.0, num), np.repeat(45.0, num)])
# x = rad * 1000.0 * np.cos(np.deg2rad(azim))
# y = rad * 1000.0 * np.sin(np.deg2rad(azim))
# lon_bound, lat_bound = trans_cart2ellps.transform(x, y)
# # ---------------------------------------------------------------------------

# Colormaps
levels_topo = np.arange(0.0, 6500.0, 500.0)
cmap_topo = cm.grayC_r
norm_topo = mpl.colors.BoundaryNorm(levels_topo, ncolors=cmap_topo.N,
                                    extend="max")
levels = np.arange(3000.0, 6200.0, 200.)
cmap = cm.roma_r
norm = mpl.colors.BoundaryNorm(levels, ncolors=cmap.N, extend="both")

# Plot maps
map_ext = [90.0, 104.8, 26.7, 35.0]
fig = plt.figure(figsize=(10.0, 17.5))
gs = gridspec.GridSpec(7, 1, left=0.1, bottom=0.1, right=0.9,
                       top=0.9, hspace=0.15, wspace=0.04,
                       height_ratios=[1, 0.01, 1, 0.003, 0.025, 0.01, 0.025])
# -----------------------------------------------------------------------------
ax = plt.subplot(gs[0], projection=ccrs.PlateCarree())
plt.pcolormesh(lon_gc2009, lat_gc2009, topo_gc2009, cmap=cmap_topo,
               norm=norm_topo, shading="auto")
data_plot = np.ma.masked_where(np.isnan(line_elev_mean), line_elev_mean)
plt.pcolormesh(lon_gc2009, lat_gc2009, data_plot, cmap=cmap, norm=norm,
               shading="auto")
# plt.plot(lon_bound, lat_bound, lw=1.5, color="black", linestyle="--")
# plt.scatter(lon_0, lat_0, color="black", s=50)
ax.set_aspect("auto")
plt.axis(map_ext)
gl = ax.gridlines(draw_labels=True, linestyle="-", linewidth=0.0)
gl.top_labels, gl.right_labels = False, False
plt.title("Permanent snow/ice line (GlobCover 2009)", fontsize=12,
          fontweight="bold", y=1.01)
# -----------------------------------------------------------------------------
ax = plt.subplot(gs[2], projection=ccrs.PlateCarree())
plt.pcolormesh(lon_gc2009, lat_gc2009, topo_gc2009, cmap=cmap_topo,
               norm=norm_topo, shading="auto")
plt.scatter(data_GAMDAM[:, 0], data_GAMDAM[:, 1], c=data_GAMDAM[:, 2],
            cmap=cmap, norm=norm, s=15)
# plt.plot(lon_bound, lat_bound, lw=1.5, color="black", linestyle="--")
# plt.scatter(lon_0, lat_0, color="black", s=50)
ax.set_aspect("auto")
plt.axis(map_ext)
gl = ax.gridlines(draw_labels=True, linestyle="-", linewidth=0.0)
gl.top_labels, gl.right_labels = False, False
plt.title("Median glacier elevation (GAMDAM)", fontsize=12,
          fontweight="bold", y=1.01)
# -----------------------------------------------------------------------------
ax = plt.subplot(gs[4])
cb = mpl.colorbar.ColorbarBase(ax, cmap=cmap, norm=norm,
                               orientation="horizontal")
cb.ax.tick_params(labelsize=10)
plt.xlabel("Permanent snow/ice line or median glacier elevation[m]",
           fontsize=10)
# -----------------------------------------------------------------------------
ax = plt.subplot(gs[6])
cb = mpl.colorbar.ColorbarBase(ax, cmap=cmap_topo, norm=norm_topo,
                               orientation="horizontal")
cb.ax.tick_params(labelsize=10)
plt.xlabel("Terrain elevation [m]", fontsize=10)
# -----------------------------------------------------------------------------
fig.savefig(path_out + "Snow_ice_threshold_elevation_map.png",
            bbox_inches="tight", dpi=300)
plt.close(fig)

# Plot histogram
fig = plt.figure(figsize=(12, 6))
data = line_elev_mean  # [mask_gc2009]
data = data[np.isfinite(data)]
l0 = plt.hist(data, bins=30, density=True, color="blue", alpha=0.7, zorder=2)
# for i in (10.0, 50.0, 90.0):  # percentiles
for i in (5.0, 50.0, 95.0):  # percentiles
    perc_value = np.percentile(data, i)
    plt.vlines(perc_value, ymin=0.0, ymax=1.0, color="black",
               lw=2.0, ls="-", zorder=3)
    txt = "~" + str(int(np.round(perc_value / 10.0, decimals=0) * 10)) + " m"
    plt.text(perc_value + 35.0, 0.0017, txt, fontsize=11, fontweight="bold")
data = data_GAMDAM[:, 2]  # [mask_GAMDAM]
data = data[np.isfinite(data)]
l1 = plt.hist(data, bins=30, density=True, color="darkgrey", alpha=0.3,
              zorder=1)
plt.axis([2500.0, 7000.0, 0.0, 0.0018])
plt.xlabel("Elevation [m]")
plt.ylabel("Density [-]")
plt.title("Distribution of permanent snow/ice line and glacier median "
          + "elevation", fontsize=13, fontweight="bold", y=1.01)
plt.legend([l0[-1], l1[-1]], ["GlobCover 2009", "GAMDAM"], loc="upper right",
           fontsize=12, frameon=False)
fig.savefig(path_out + "Snow_ice_threshold_elevation_histogram.pdf",
            bbox_inches="tight")
plt.close(fig)

# Notes
# - assume that median elevation represents the equilibrium line altitude
#   -> Braithwaite and Raper (2009), Braithwaite (2015)

###############################################################################
# Interpolate equilibrium line altitude (ELA) by inverse distance weighting
###############################################################################

# # Output grid
# lon_out = np.linspace(85.0, 110.0, 250)
# lat_out = np.linspace(20.0, 40.0, 200)
# lon_out, lat_out = np.meshgrid(lon_out, lat_out)
#
# # Compute map projection coordinates (orthographic projection)
# crs_wgs84 = CRS.from_epsg(4326)
# crs_ortho = CRS.from_proj4("+proj=ortho +lat_0=" + str(lat.mean())
#                            + " +lon_0=" + str(lon.mean()) + " +x_0=0 +y_0=0")
# transformer = Transformer.from_crs(crs_wgs84, crs_ortho, always_xy=True)
# lon_2d, lat_2d = np.meshgrid(lon, lat)
# x, y = transformer.transform(lon_2d, lat_2d)
# x_pts, y_pts = transformer.transform(data_glac[:, 0], data_glac[:, 1])
# x_out, y_out = transformer.transform(lon_out, lat_out)
#
# # Build k-d tree
# glac_pts = np.hstack((x_pts[:, np.newaxis], y_pts[:, np.newaxis]))
# tree = cKDTree(glac_pts)
#
# # Interpolate with Inverse distance weighting (IDW)
# rad_search = 200000.0  # search radius [m]
# elev_ip = np.empty(x_out.shape, dtype=np.float32)
# elev_ip.fill(np.nan)
# for i in range(x_out.shape[0]):
#     for j in range(x_out.shape[1]):
#         pts_cc = np.array([x_out[i, j], y_out[i, j]],
#                           dtype=np.float32).reshape(1, 2)
#         ind = tree.query_ball_point(pts_cc, r=rad_search)[0]
#         if len(ind) > 0:
#             dist = np.sqrt(((glac_pts[ind, :] - pts_cc) ** 2).sum(axis=1))
#             elev_ip[i, j] = np.average(data_glac[ind, 2],
#                                        weights=(1.0 / dist ** 0.0))
