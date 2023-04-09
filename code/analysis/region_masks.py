# Description: Compute Eastern Tibet / Hengduan Mountains regions masks for
#              different products
#
# Author: Christian R. Steger, April 2023

# Load modules
import os
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.gridspec as gridspec
from shapely.geometry import Polygon, LineString
from shapely.geometry import shape
from shapely.ops import unary_union
from shapely.ops import transform
from shapely.ops import split
import cartopy.crs as ccrs
from cartopy.io import shapereader
import cartopy.feature as cfeature
import fiona
from descartes import PolygonPatch
from pyproj import CRS, Transformer
import utilities
from cmcrameri import cm
from netCDF4 import Dataset
import pickle

mpl.style.use("classic")

# Path to folders
path_out = "/project/pr133/csteger/Data/Model/BECCY/region_masks/"
# output directory

###############################################################################
# Define regions in rotated latitude/longitude coordinates
###############################################################################

# Load CTR04 topography
ds = xr.open_dataset("/project/pr133/rxiang/data/extpar/"
                     + "extpar_BECCY_4.4km_merit_unmod_topo.nc")
ds = ds.isel(rlon=slice(30, 680), rlat=slice(30, 680))
rlon = ds["rlon"].values
rlat = ds["rlat"].values
topo = ds["HSURF"].values
lsm = ds["FR_LAND"].values
crs_rot_pole = ccrs.RotatedPole(
    pole_longitude=ds["rotated_pole"].grid_north_pole_longitude,
    pole_latitude=ds["rotated_pole"].grid_north_pole_latitude)
ds.close()

# Define regions (polygons)
regions = {}
box = (-24.84, -6.67 - 0.9, -3.00, 12.57)  # (x_min, y_min, x_max, y_max)
print("Boundary with (right): %.2f" % (rlon[-1] - box[2]) + " deg")
print("Boundary with (top): %.2f" % (rlat[-1] - box[3]) + " deg")
regions["ET"] = utilities.grid.polygon_rectangular(box, spacing=0.01)
# Eastern Tibet
box = (-19.0, -6.67, -9.0, 6.0)
regions["HM"] = utilities.grid.polygon_rectangular(box, spacing=0.01)
# Hengduan Mountains

# Get borders of certain countries
countries = ("India", "Myanmar")
file_shp = shapereader.natural_earth("10m", "cultural", "admin_0_countries")
ds = fiona.open(file_shp)
geom_names = [i["properties"]["NAME"] for i in ds]
poly_count = unary_union([shape(ds[geom_names.index(i)]["geometry"])
                          for i in countries])  # merge all polygons
crs_count = CRS.from_string(ds.crs["init"])
ds.close()

# Load CTR04 rainy season precipitation
ds = xr.open_dataset("/project/pr133/rxiang/data/cosmo/EAS04_ctrl/mon/"
                     + "TOT_PREC/2001-2005.TOT_PREC.nc")
ds = ds.sel(time=(ds["time.month"] >= 5) & (ds["time.month"] <= 9))
precip = ds["TOT_PREC"].values.mean(axis=0)  # [mm day-1]
ds.close()

# Transform country polygon to rotated latitude/longitude coordinates
project = Transformer.from_crs(crs_count, CRS.from_user_input(crs_rot_pole),
                               always_xy=True).transform
poly_count_rot = transform(project, poly_count)

# Intersect polygons
regions["HMU"] = regions["HM"].intersection(poly_count_rot)
# Hengduan Mountains Upstream
regions["HMC"] = regions["HM"].difference(regions["HMU"])
# Hengduan Mountains Centre
line = LineString([(-28.24, -3.6), (-2.28, -3.6)])
geom_split = split(regions["HMU"], line)
regions["HMUS"] = geom_split[0]  # South
regions["HMUN"] = geom_split[1]  # North

# Save region polygons to disk
path = path_out + "region_polygons/"
if not os.path.isdir(path_out + "region_polygons/"):
    os.mkdir(path)
for i in regions.keys():
    with open(path + i + "_rot_coord.poly", "wb") as file:
        pickle.dump(regions[i], file, pickle.HIGHEST_PROTOCOL)

# # Test load of polygon
# regions = {}
# for i in ["ET", "HM", "HMU", "HMC", "HMUS", "HMUN"]:
#     with open(path + i + "_rot_coord.poly", "rb") as file:
#         regions[i] = pickle.load(file)

###############################################################################
# Plot regions
###############################################################################

# Colormap
levels = np.arange(0.0, 6500.0, 500.0)
ticks = np.arange(0.0, 6500.0, 500.0)
cmap = utilities.plot.truncate_colormap(cm.bukavu, 0.5, 1.0)
norm = mpl.colors.BoundaryNorm(levels, ncolors=cmap.N, extend="max")
# levels = np.arange(0.0, 21.0, 1.0)
# ticks = np.arange(0.0, 22.0, 2.0)
# cmap = cm.lapaz_r
# norm = mpl.colors.BoundaryNorm(levels, ncolors=cmap.N, extend="max")

# Map plot
fig = plt.figure(figsize=(9, 8), dpi=200)
gs = gridspec.GridSpec(1, 2, left=0.1, bottom=0.1, right=0.9, top=0.9,
                       hspace=0.05, wspace=0.07, width_ratios=[1, 0.05])
# -----------------------------------------------------------------------------
ax = plt.subplot(gs[0], projection=crs_rot_pole)
ax.set_facecolor(cm.bukavu(0.4))
data_plot = np.ma.masked_where(lsm < 0.5, topo)
# data_plot = precip
plt.pcolormesh(rlon, rlat, data_plot, cmap=cmap, norm=norm,
               transform=crs_rot_pole)
poly_plot = PolygonPatch(regions["ET"], facecolor="none",
                         edgecolor="black", lw=2.5, ls="--",
                         transform=crs_rot_pole, zorder=3)
ax.add_patch(poly_plot)
poly_plot = PolygonPatch(regions["HM"], facecolor="none",
                         edgecolor="black", lw=2.5, ls="--",
                         transform=crs_rot_pole, zorder=3)
ax.add_patch(poly_plot)
for i in ["HMC", "HMUS", "HMUN"]:
    poly_plot = PolygonPatch(regions[i], facecolor="none",
                             edgecolor="orangered", lw=3.0, ls="-",
                             transform=crs_rot_pole, zorder=2)
    ax.add_patch(poly_plot)
# -----------------------------------------------------------------------------
ax.coastlines(resolution="10m")
bord_10m = cfeature.NaturalEarthFeature("cultural",
                                        "admin_0_boundary_lines_land",
                                        "10m",
                                        edgecolor="black",
                                        facecolor="none", zorder=1)
ax.add_feature(bord_10m)
ax.set_aspect("auto")
gl = ax.gridlines(crs=ccrs.PlateCarree(),
                  xlocs=np.arange(0.0, 180.0, 5.0),
                  ylocs=np.arange(0.0, 90.0, 5.0),
                  linewidth=1, color="darkgrey", alpha=1.0, linestyle=":",
                  draw_labels=True, dms=True, x_inline=False, y_inline=False)
gl.top_labels = False
gl.right_labels = False
# -----------------------------------------------------------------------------
ax = plt.subplot(gs[1])
cb = mpl.colorbar.ColorbarBase(ax, cmap=cmap, norm=norm,
                               ticks=ticks, orientation="vertical")
# -----------------------------------------------------------------------------
# plt.show()
fig.savefig(path_out + "Region_masks.png", dpi=200, bbox_inches="tight")
plt.close(fig)

# Full names of regions
regions_fn = {"ET": "Eastern Tibet",
              "HM": "Hengduan Mountains",
              "HMU": "Hengduan Mountains upstream",
              "HMC": "Hengduan Mountains centre",
              "HMUS": "Hengduan Mountains upstream south",
              "HMUN": "Hengduan Mountains upstream north"}

###############################################################################
# Compute region masks for different products
###############################################################################

# List of products (with file and geographic reference)
products = {
    # -------------------------------------------------------------------------
    "APHRODITE": {"file": "/project/pr133/rxiang/data/obs/pr/APHRO/"
                          + "APHRO_2001-2005_ymonmean.nc",
                  "geo_ref": CRS.from_epsg(4326)},
    "GPCC": {"file": "/project/pr133/rxiang/data/obs/pr/GPCC/mo/"
                     + "GPCC.2001-2005.025.mon.nc",
             "geo_ref": CRS.from_epsg(4326)},
    "IMERG": {"file": "/project/pr133/rxiang/data/obs/pr/IMERG/day_old/"
                      + "IMERG.ydaymean.2001-2005.mon.nc4",
              "geo_ref": CRS.from_epsg(4326)},
    "ERA5": {"file": "/project/pr133/rxiang/data/era5/pr/mo/"
                     + "era5.mo.2001-2005.mon.nc",
             "geo_ref": CRS.from_epsg(4326)},
    "PBCOR": {"file": "/project/pr133/csteger/Data/Observations/"
                      + "PBCOR_V1.0/CHELSA_V12.nc",
              "geo_ref": CRS.from_epsg(4326)},
    # -------------------------------------------------------------------------
    "CRU": {"file": "/project/pr133/rxiang/data/obs/tmp/cru/mo/"
                    + "cru.2001-2005.05.mon.nc",
            "geo_ref": CRS.from_epsg(4326)},
    # -------------------------------------------------------------------------
    "CTRL04": {"file": "/project/pr133/rxiang/data/cosmo/EAS04_ctrl/mon/"
                       + "TOT_PREC/2001-2005.TOT_PREC.nc",
               "geo_ref": CRS.from_user_input(crs_rot_pole)},
    "CTRL11": {"file": "/project/pr133/rxiang/data/cosmo/EAS11_ctrl/mon/"
                       + "TOT_PREC/2001-2005.TOT_PREC.nc",
               "geo_ref": CRS.from_user_input(crs_rot_pole)},
    # -------------------------------------------------------------------------
    }

# Loop through products
for i in products.keys():

    # Load spatial coordinates
    ds = xr.open_dataset(products[i]["file"])
    if "rlon" in list(ds.coords):
        x = ds["rlon"].values
        y = ds["rlat"].values
    elif "lon" in list(ds.coords):
        x = ds["lon"].values
        y = ds["lat"].values
    elif "longitude" in list(ds.coords):
        x = ds["longitude"].values
        y = ds["latitude"].values
    else:
        raise ValueError("Unknown spatial coordinates")

    if i == "IMERG":
        x[x < 0.0] += 360.0  # make latitude axis monotonically increasing

    # Transform region polygons and compute (binary) masks
    project = Transformer.from_crs(CRS.from_user_input(crs_rot_pole),
                                   products[i]["geo_ref"],
                                   always_xy=True).transform
    x_edge, y_edge = np.meshgrid(*utilities.grid.coord_edges(x, y, atol=1e-04))
    region_masks = {}
    for j in regions.keys():
        region_trans = transform(project, regions[j])
        area_frac = utilities.grid.polygon_inters_exact(
            x_edge, y_edge, region_trans,
            agg_cells=np.array([10, 5, 2]))
        region_masks[j] = (area_frac > 0.5).astype(np.int32)

    # Save to NetCDF file
    ncfile = Dataset(filename=path_out + i + "_region_masks.nc", mode="w",
                     format="NETCDF4")
    ncfile.createDimension(dimname="y", size=area_frac.shape[0])
    ncfile.createDimension(dimname="x", size=area_frac.shape[1])
    for j in region_masks.keys():
        nc_data = ncfile.createVariable(varname=j, datatype="i",
                                        dimensions=("y", "x"))
        nc_data.units = "-"
        nc_data.long_name = regions_fn[j]
        nc_data[:] = region_masks[j]
    ncfile.close()

    print("File " + i + "_region_masks.nc created")

###############################################################################
# Create nice overview plot with regions and precipitation patterns
###############################################################################

# Load IMERG data for rainy season (May -> September)
ds = xr.open_dataset("/project/pr133/rxiang/data/obs/pr/IMERG/day_old/"
                     + "IMERG.ydaymean.2001-2005.mon.nc4")
ds = ds.isel(lon=slice(250, 750), lat=slice(100, 500))
ds = ds.sel(time=(ds["time.month"] >= 5) & (ds["time.month"] <= 9))
prec = ds["pr"].values.mean(axis=0)  # [mm day-1]
lon = ds["lon"].values
lat = ds["lat"].values
ds.close()

# # Load ERA5 data for rainy season (May -> September)
# ds = xr.open_dataset("/project/pr133/rxiang/data/era5/pr/mo/"
#                      + "era5.mo.2001-2005.mon.nc")
# ds = ds.isel(longitude=slice(200, 500))
# ds = ds.sel(time=(ds["time.month"] >= 5) & (ds["time.month"] <= 9))
# prec = ds["tp"].values.mean(axis=0) * 1000.0  # [mm/day-1]
# lon = ds["longitude"].values
# lat = ds["latitude"].values
# ds.close()

# Colormap
levels = np.arange(0.0, 21.0, 1.0)
ticks = np.arange(0.0, 22.0, 2.0)
cmap = cm.lapaz_r
norm = mpl.colors.BoundaryNorm(levels, ncolors=cmap.N, extend="max")

# Regional labels settings
labels = {"ET": {"pos": (86.5, 37.0), "color": "black"},
          "HM": {"pos": (94.8, 31.7), "color": "black"},
          "HMC": {"pos": (101.5, 32.5), "color": "red"},
          "HMUS": {"pos": (96.2, 21.8), "color": "red"},
          "HMUN": {"pos": (95.4, 26.2), "color": "red"}}

# Map plot
fig = plt.figure(figsize=(9, 8), dpi=200)
gs = gridspec.GridSpec(1, 2, left=0.1, bottom=0.1, right=0.9, top=0.9,
                       hspace=0.05, wspace=0.07, width_ratios=[1, 0.05])
# -----------------------------------------------------------------------------
ax = plt.subplot(gs[0], projection=crs_rot_pole)
ax.set_facecolor(cm.bukavu(0.4))
plt.pcolormesh(lon, lat, prec, cmap=cmap, norm=norm,
               transform=ccrs.PlateCarree())
poly_plot = PolygonPatch(regions["ET"], facecolor="none",
                         edgecolor="black", lw=2.5, ls="--",
                         transform=crs_rot_pole, zorder=3)
ax.add_patch(poly_plot)
poly_plot = PolygonPatch(regions["HM"], facecolor="none",
                         edgecolor="black", lw=2.5, ls="--",
                         transform=crs_rot_pole, zorder=3)
ax.add_patch(poly_plot)
for i in ["HMC", "HMUS", "HMUN"]:
    poly_plot = PolygonPatch(regions[i], facecolor="none",
                             edgecolor="orangered", lw=3.0, ls="-",
                             transform=crs_rot_pole, zorder=2)
    ax.add_patch(poly_plot)
for i in ["ET", "HM", "HMC", "HMUS", "HMUN"]:
    plt.text(*labels[i]["pos"], i, color=labels[i]["color"], fontsize=11,
             fontweight="bold", transform=ccrs.PlateCarree())
plt.title("IMERG precipitation (2001 - 2005; rainy season) [mm day-1]", y=1.01,
          fontsize=14, fontweight="bold")
# -----------------------------------------------------------------------------
ax.coastlines(resolution="10m")
bord_10m = cfeature.NaturalEarthFeature("cultural",
                                        "admin_0_boundary_lines_land",
                                        "10m",
                                        edgecolor="black",
                                        facecolor="none", zorder=1)
ax.add_feature(bord_10m)
ax.set_aspect("auto")
gl = ax.gridlines(crs=ccrs.PlateCarree(),
                  xlocs=np.arange(0.0, 180.0, 5.0),
                  ylocs=np.arange(0.0, 90.0, 5.0),
                  linewidth=1, color="darkgrey", alpha=1.0, linestyle=":",
                  draw_labels=True, dms=True, x_inline=False, y_inline=False)
gl.top_labels = False
gl.right_labels = False
ax.set_extent([-29.5, 0.5, -11.5, 15.5], crs=crs_rot_pole)
# -----------------------------------------------------------------------------
ax = plt.subplot(gs[1])
cb = mpl.colorbar.ColorbarBase(ax, cmap=cmap, norm=norm,
                               ticks=ticks, orientation="vertical")
# -----------------------------------------------------------------------------
# plt.show()
fig.savefig(path_out + "Region_masks_labels.png", dpi=200, bbox_inches="tight")
plt.close(fig)
