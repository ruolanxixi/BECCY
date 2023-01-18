# Description: Compute envelop topography part I (generate envelope topography
#              with convex hull and increased surface curvature)
#
# Author: Christian R. Steger, IAC ETH Zurich

# Load modules
import os
import sys
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.gridspec as gridspec
import matplotlib.ticker as mticker
from scipy.interpolate import LinearNDInterpolator
from pyproj import CRS, Transformer
import time
from cmcrameri import cm
import cartopy.crs as ccrs
from scipy.spatial import ConvexHull
from netCDF4 import Dataset

mpl.style.use("classic")

# Change latex fonts
mpl.rcParams["mathtext.fontset"] = "custom"
# custom mathtext font (set default to Bitstream Vera Sans)
mpl.rcParams["mathtext.default"] = "rm"
mpl.rcParams["mathtext.rm"] = "Bitstream Vera Sans"

# Load required functions
sys.path.append("/Users/csteger/Downloads/BECCY/code/")
from auxiliary import truncate_colormap
from auxiliary import spat_agg_2d, spat_agg_1d

###############################################################################
# Settings
###############################################################################

# Paths
path_dem = "/Users/csteger/Dropbox/IAC/Data/DEMs/MERIT/Tiles/"
path_plot = "/Users/csteger/Dropbox/IAC/Plots/BECCY/Topo_envelop/"
path_out = "/Users/csteger/Dropbox/IAC/Data/Model/BECCY/env_topo_raw/"

# Constants
rad_earth = 6370997.0  # default PROJ sphere radius [m]

# Settings
fac_curv = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
file_netcdf = "MERIT_envelope_topo.nc"

###############################################################################
# Compute and save envelope topographies for different surface curvatures
###############################################################################

# Load MERIT DEM
files_dem = ("MERIT_N60-N30_E060-E090.nc", "MERIT_N30-N00_E060-E090.nc",
             "MERIT_N60-N30_E090-E120.nc", "MERIT_N30-N00_E090-E120.nc")
ds = xr.open_mfdataset([path_dem + i for i in files_dem])
# ds = ds.sel(lon=slice(96.0, 97.5), lat=slice(28.5, 27.5))  # (test; small)
# ds = ds.sel(lon=slice(96, 100), lat=slice(30, 27))  # (test; medium)
# ds = ds.sel(lon=slice(94.5, 107.3), lat=slice(34.0, 20.2))  # (old)
ds = ds.sel(lon=slice(94.4, 107.4), lat=slice(34.1, 20.1))  # (new; +0.1 deg)
topo = ds["Elevation"].values  # 32-bit float
lon = ds["lon"].values
lat = ds["lat"].values
ds.close()
print("Size of DEM: " + str(topo.shape))

# # Artificial test data
# lon = np.linspace(95.0, 105.0, 20)
# lat = np.linspace(33.0, 21.0, 20)
# topo = np.zeros((len(lat), len(lon)), dtype=np.float32)

# Pre-process MERIT data
mask_ocean = np.isnan(topo)
print("Number of ocean grid cells: " + str(mask_ocean.sum()))
topo[mask_ocean] = 0.0  # set ocean grid cells to 0.0 m
print("Minimal DEM elevation: " + str(topo.min()) + " m")

# -----------------------------------------------------------------------------
# Generate and save envelope topography
# -----------------------------------------------------------------------------

# Create NetCDF file
topo_fill_val = -32767
ncfile = Dataset(filename=path_out + file_netcdf, mode="w")
ncfile.createDimension(dimname="lat", size=topo.shape[0])
nc_lat = ncfile.createVariable(varname="lat", datatype="f", dimensions="lat")
nc_lat.units = "degrees_north"
nc_lat[:] = lat
ncfile.createDimension(dimname="lon", size=topo.shape[1])
nc_lat = ncfile.createVariable(varname="lon", datatype="f", dimensions="lon")
nc_lat.units = "degrees_east"
nc_lat[:] = lon
nc_data = ncfile.createVariable(varname="Elevation", datatype="i2",
                                dimensions=("lat", "lon"),
                                fill_value=topo_fill_val)
nc_data.units = "m"
nc_data.long_name = "Unmodified topography"
nc_data[:] = np.where(mask_ocean, topo_fill_val, topo).astype(np.int16)
ncfile.close()

# Loop through different curvature factors
for i in fac_curv:

    print((" Surface curvature factor %.1f" % i + " ").center(79, "#"))

    # Construct coordinates centred at prime meridian / equator
    dx = (2.0 * np.pi * rad_earth * np.cos(np.deg2rad(lat.mean()))) \
        / 360.0 * np.abs(np.diff(lon).mean())  # [m]
    dy = (2.0 * np.pi * rad_earth) / 360.0 * np.abs(np.diff(lat).mean())  # [m]
    lon_rang = (360.0 / (2.0 * np.pi * rad_earth)) * dx * topo.shape[1] * i
    lat_rang = (360.0 / (2.0 * np.pi * rad_earth)) * dy * topo.shape[0] * i
    print("Longitude/latitude range [deg]: %.3f" % (lon[-1] - lon[0]) + ", "
          + "%.3f" % (lat[0] - lat[-1]))
    print("Longitude/latitude range (equator / prime meridian cent.)"
          + " [deg]: %.3f" % lon_rang + ", " + "%.3f" % lat_rang)
    if (lon_rang > 179.0) or (lat_rang > 179.0):
        raise ValueError("Curvature factor too large for domain")
    lon_cen_pme = np.linspace(-(lon_rang / 2.0), (lon_rang / 2.0),
                              topo.shape[1])
    lat_cen_pme = np.linspace((lat_rang / 2.0), -(lat_rang / 2.0),
                              topo.shape[0])

    # Compute coordinates in map projection
    time_beg = time.time()
    crs_latlon = CRS.from_dict({"proj": "latlon", "ellps": "sphere",
                                "R": rad_earth})
    crs_ecef = CRS.from_dict({"proj": "geocent", "ellps": "sphere",
                              "R": rad_earth})
    lon_2d, lat_2d = np.meshgrid(lon_cen_pme, lat_cen_pme)
    transformer = Transformer.from_crs(crs_latlon, crs_ecef, always_xy=True)
    x_ecef, y_ecef, z_ecef = transformer.transform(lon_2d, lat_2d, topo)
    print("Computing ECEF coordinates: " + "%.2f" % (time.time() - time_beg)
          + " sec")

    # Add Earth centre to points
    time_beg = time.time()
    pts = np.hstack((x_ecef.ravel()[:, np.newaxis],
                     y_ecef.ravel()[:, np.newaxis],
                     z_ecef.ravel()[:, np.newaxis]))
    pts_earth_cen = np.array([0.0, 0.0, 0.0]).reshape(1, 3)
    pts = np.vstack((pts, pts_earth_cen))
    del x_ecef, y_ecef, z_ecef
    print("Rearrange data: " + "%.2f" % (time.time() - time_beg) + " sec")

    # Compute convex hull
    if pts.shape[0] >= 2147483647:
        raise ValueError("Number of points will exceed 32 bit integer size of "
                         + "convex hull simplices")
    time_beg = time.time()
    convex_hull = ConvexHull(pts)
    print("Computing convex hull: " + "%.2f" % (time.time() - time_beg)
          + " sec")
    print("Number of triangles: " + str(convex_hull.nsimplex))

    # Remove triangles facing sideways
    mask_side = np.ones(topo.shape, dtype=bool)
    mask_side[1:-1, 1:-1] = False
    mask_side = np.append(mask_side.ravel(), np.array([True]))
    ind_del = []
    for j in range(convex_hull.nsimplex):
        ind = convex_hull.simplices[j, :]
        if np.all(mask_side[ind]):
            ind_del.append(j)
    simplices_rel = np.delete(convex_hull.simplices, ind_del, axis=0)
    print("Number of relevant triangles: " + str(simplices_rel.shape[0]))
    del mask_side

    # # 3D test plot (only affordable for small DEM size!)
    # fig = plt.figure()
    # ax = fig.add_subplot(projection="3d")
    # ax.scatter(pts[:, 0], pts[:, 1], pts[:, 2])
    # ax.scatter(pts[-1, 0], pts[-1, 1], pts[-1, 2])
    # for j in range(convex_hull.nsimplex):
    #     ind_cl = np.append(convex_hull.simplices[j, :],
    #                        convex_hull.simplices[j, 0])
    #     plt.plot(pts[ind_cl, 0], pts[ind_cl, 1], pts[ind_cl, 2],
    #              color="red")
    # for j in range(simplices_rel.shape[0]):
    #     ind_cl = np.append(simplices_rel[j, :], simplices_rel[j, 0])
    #     plt.plot(pts[ind_cl, 0], pts[ind_cl, 1], pts[ind_cl, 2],
    #              color="green")
    # ax.set_xlabel("x-axis")
    # ax.set_ylabel("y-axis")
    # ax.set_zlabel("z-axis")

    # Get relevant nodes
    ind_nodes = convex_hull.vertices
    ind_del = np.where(pts[ind_nodes, :].sum(axis=1) == 0.0)[0][0]
    ind_nodes = np.delete(ind_nodes, ind_del)
    del pts, convex_hull

    # Rasterise
    time_beg = time.time()
    pts_nodes = np.hstack((lon_2d.ravel()[ind_nodes][:, np.newaxis],
                           lat_2d.ravel()[ind_nodes][:, np.newaxis]))
    interp_loc = LinearNDInterpolator(pts_nodes, topo.ravel()[ind_nodes])
    # -> assume that lon/lat-coordinates represent a planar grid (map
    #    projection also introduce slight inaccuracy for larger domain)
    topo_env = interp_loc(lon_2d, lat_2d).astype(np.float32)
    del pts_nodes, lon_2d, lat_2d, interp_loc
    print("Rasterise triangle mesh data : %.2f" % (time.time() - time_beg)
          + " sec")
    print("Maximal elevation of ocean grid cell: %.1f"
          % np.nanmax(topo_env[mask_ocean]) + " m")

    # Append to NetCDF file
    ncfile = Dataset(filename=path_out + file_netcdf, mode="a")
    nc_data = ncfile.createVariable(varname="Elev_curv_%.1f" % i,
                                    datatype="i2", dimensions=("lat", "lon"),
                                    fill_value=topo_fill_val)
    nc_data.units = "m"
    nc_data.long_name = "Envelope topography (surface curvature " \
                        + "factor: %.1f" % i + ")"
    nc_data[:] = np.where(np.isnan(topo_env), topo_fill_val, topo_env) \
        .astype(np.int16)
    ncfile.close()

    np.savez(path_out + "aux_curv_%.1f" % i + ".npz", simplices=simplices_rel,
             indices_nodes=ind_nodes)

    del ind_nodes, simplices_rel

###############################################################################
# Plot selected envelope topographies
###############################################################################

# Select experiments
fac_curv_sel = [1.0, 3.0, 5.0, 7.0]

# Load topographies
topo_env = {}
ds = xr.open_dataset(path_out + file_netcdf)
for i in fac_curv_sel:
    topo_env[i] = ds["Elev_curv_%.1f" % i].values  # 32-bit float
ds.close()

# Spatially aggregate data for plotting
agg_num = 10
slic = (slice(0, (topo.shape[0] // agg_num) * agg_num),
        slice(0, (topo.shape[1] // agg_num) * agg_num))
topo_agg = spat_agg_2d(topo[slic], agg_num, agg_num, operation="mean")
mask_ocean_agg = spat_agg_2d(mask_ocean[slic], agg_num, agg_num,
                             operation="sum")
lon_agg = spat_agg_1d(lon[slic[1]], agg_num, operation="mean")
lat_agg = spat_agg_1d(lat[slic[0]], agg_num, operation="mean")
topo_env_agg = {}
for i in fac_curv_sel:
    topo_env_agg[i] = spat_agg_2d(topo_env[i][slic], agg_num, agg_num,
                                  operation="mean")

# Colormap
levels_topo = np.arange(0., 5750.0, 250.0)
ticks_topo = np.arange(0., 6000.0, 1000.0)
cmap_topo = truncate_colormap(cm.bukavu, 0.55, 1.0)
norm_topo = mpl.colors.BoundaryNorm(levels_topo, ncolors=cmap_topo.N,
                                    extend="max")
levels_diff = np.arange(0., 3250.0, 250.0)
ticks_diff = np.arange(0., 3500.0, 500.0)
cmap_diff = cm.davos_r
norm_diff = mpl.colors.BoundaryNorm(levels_diff, ncolors=cmap_diff.N,
                                    extend="both")

# -----------------------------------------------------------------------------
# Map plot
# -----------------------------------------------------------------------------

lon_2d, lat_2d = np.meshgrid(lon, lat)

# Plot
fig = plt.figure(figsize=(11.6, 16.5))  # 4 rows
# fig = plt.figure(figsize=(11.3, 19.5))  # 5 rows
gs = gridspec.GridSpec(len(fac_curv_sel) + 1, 3, left=0.1, bottom=0.1,
                       right=0.9, top=0.9, hspace=0.05, wspace=0.05,
                       height_ratios=[1.0] * len(fac_curv_sel) + [0.06])
# -----------------------------------------------------------------------------
count = 0
for i in fac_curv_sel:
    ax = plt.subplot(gs[count, 0], projection=ccrs.PlateCarree())
    if topo.size <= (2000 * 2000):
        plt.pcolormesh(lon, lat, topo, cmap=cmap_topo, norm=norm_topo,
                       shading="auto")
    else:
        plt.pcolormesh(lon_agg, lat_agg, topo_agg, cmap=cmap_topo,
                       norm=norm_topo, shading="auto")
    data_aux = np.load(path_out + "aux_curv_%.1f" % i + ".npz")
    if i == 1.0:  # set manually
        for j in range(data_aux["simplices"].shape[0]):
            ind = data_aux["simplices"][j, :]
            ind_cl = np.append(ind, ind[0])
            rang_0 = np.ptp(ind // lon_2d.shape[1])
            rang_1 = np.ptp(ind % lon_2d.shape[1])
            if (rang_0 > 1) or (rang_1 > 1):
                plt.plot(lon_2d.ravel()[ind_cl], lat_2d.ravel()[ind_cl],
                         "black", lw=1.0)
                # -> do not plot triangles that are formed by directly
                #    neighbouring grid cells (occur over ocean) -> too
                #    expensive to plot
    ind_nodes = data_aux["indices_nodes"]
    data_aux.close()
    if i in (1.0, 3.0, 4.0):  # set manually
        plt.scatter(lon_2d.ravel()[ind_nodes], lat_2d.ravel()[ind_nodes], s=1,
                    color="black")
    ax.set_aspect("auto")
    plt.axis([lon.min(), lon.max(), lat.min(), lat.max()])
    plt.text(-0.05, 0.5, "Curvature factor: %.1f" % i, fontsize=11,
             fontweight="bold",
             horizontalalignment="center", verticalalignment="center",
             rotation=90, transform=ax.transAxes)
    if count == 0:
        plt.title("Present-day topography [m]", fontsize=11, fontweight="bold",
                  y=1.01)
    count += 1
# -----------------------------------------------------------------------------
ax = plt.subplot(gs[-1, 0])
mpl.colorbar.ColorbarBase(ax, cmap=cmap_topo, norm=norm_topo, ticks=ticks_topo,
                          orientation="horizontal")
# -----------------------------------------------------------------------------
cmap_oc = mpl.colors.ListedColormap(["grey"])
bounds_red = [0.5, 1.5]
norm_red = mpl.colors.BoundaryNorm(bounds_red, cmap_oc.N)
count = 0
for i in fac_curv_sel:
    ax = plt.subplot(gs[count, 1], projection=ccrs.PlateCarree())
    if topo.size <= (2000 * 2000):
        plt.pcolormesh(lon, lat, topo_env[i], cmap=cmap_topo, norm=norm_topo,
                       shading="auto")
        data_plot = np.ma.masked_where(mask_ocean == 0, np.ones_like(topo))
        plt.pcolormesh(lon, lat, data_plot, cmap=cmap_oc, norm=norm_red,
                       shading="auto")
    else:
        plt.pcolormesh(lon_agg, lat_agg, topo_env_agg[i], cmap=cmap_topo,
                       norm=norm_topo, shading="auto")
        data_plot = np.ma.masked_where(mask_ocean_agg == 0,
                                       np.ones_like(topo_agg))
        plt.pcolormesh(lon_agg, lat_agg, data_plot, cmap=cmap_oc,
                       norm=norm_red, shading="auto")
    ax.set_aspect("auto")
    plt.axis([lon.min(), lon.max(), lat.min(), lat.max()])
    if count == 0:
        plt.title("Raw envelope topography [m]", fontsize=11,
                  fontweight="bold", y=1.01)
    count += 1
# -----------------------------------------------------------------------------
count = 0
for i in fac_curv_sel:
    ax = plt.subplot(gs[count, 2], projection=ccrs.PlateCarree())
    if topo.size <= (2000 * 2000):
        plt.pcolormesh(lon, lat, (topo_env[i] - topo), cmap=cmap_diff,
                       norm=norm_diff, shading="auto")
    else:
        plt.pcolormesh(lon_agg, lat_agg, (topo_env_agg[i] - topo_agg),
                       cmap=cmap_diff, norm=norm_diff, shading="auto")
    ax.set_aspect("auto")
    plt.axis([lon.min(), lon.max(), lat.min(), lat.max()])
    if count == 0:
        plt.title("Elevation difference [m]", fontsize=11, fontweight="bold",
                  y=1.01)
    count += 1
# -----------------------------------------------------------------------------
ax = plt.subplot(gs[-1, 2])
mpl.colorbar.ColorbarBase(ax, cmap=cmap_diff, norm=norm_diff, ticks=ticks_diff,
                          orientation="horizontal")
# -----------------------------------------------------------------------------
fig.savefig(path_plot + "Envelope_convex_hull_map.png", dpi=300,
            bbox_inches="tight")
plt.close(fig)

del lon_2d, lat_2d

# -----------------------------------------------------------------------------
# Cross-section plot
# -----------------------------------------------------------------------------

# Latitude transect values
lat_trans = [29.3, 27.5, 26.5, 25.5]  # [degree]
ind_lat = [np.argmin(np.abs(lat - i)) for i in lat_trans]

# Plot
fig = plt.figure(figsize=(16.0, 13.5))
gs = gridspec.GridSpec(len(ind_lat) + 2, 3, left=0.1, bottom=0.1, right=0.9,
                       top=0.9, hspace=0.05, wspace=0.05,
                       height_ratios=[1.0] * len(ind_lat) + [0.2, 1.4],
                       width_ratios=[1.0, 0.5, 1.0])
# -----------------------------------------------------------------------------
count_trans = 0
# cols = ["peru", "dodgerblue", "forestgreen", "orangered", "darkviolet"]
cols = cm.roma([0.0, 0.2, 0.4, 0.6, 0.8])
for i in ind_lat:
    ax = plt.subplot(gs[count_trans, :])
    plt.fill_between(lon, 0.0, topo[i, :], color="lightgray")
    plt.plot(lon, topo[i, :], color="black", lw=1.0)
    count_exp = 0
    for j in fac_curv_sel:
        plt.plot(lon, topo_env[j][i, :], color=cols[count_exp, :], lw=1.5,
                 label="%.1f" % j)
        count_exp += 1
    if count_trans != len(ind_lat) - 1:
        plt.xticks(np.arange(90.0, 111.0, 1.0), "")
    else:
        plt.xticks(np.arange(90.0, 111.0, 1.0))
        plt.xlabel("Longitude [$^{\circ}$]")
    if count_trans == 0:
        plt.legend(loc="upper right", frameon=False, fontsize=11,
                   title="Curvature factor:", title_fontsize=11,
                   ncol=2)
    t = plt.text(0.07, 0.90, "Latitude: %.2f" % lat[i] + "$^{\circ}$N",
                 fontsize=11,
                 horizontalalignment="center", verticalalignment="center",
                 transform=ax.transAxes, fontweight="bold")
    t.set_bbox(dict(facecolor="white", alpha=0.7, edgecolor="white"))
    plt.yticks(range(1000, 9000, 1000))
    plt.ylabel("Elevation [m]")
    plt.axis([lon.min(), lon.max(), 0.0, 7500.0])
    count_trans += 1
# -----------------------------------------------------------------------------
ax = plt.subplot(gs[-1, 1], projection=ccrs.PlateCarree())
if topo.size <= (2000 * 2000):
    plt.pcolormesh(lon, lat, topo, cmap=cmap_topo, norm=norm_topo,
                   shading="auto", rasterized=True)
else:
    plt.pcolormesh(lon_agg, lat_agg, topo_agg, cmap=cmap_topo, norm=norm_topo,
                   shading="auto", rasterized=True)
for i in ind_lat:
    plt.plot([lon[0], lon[-1]], [lat[i], lat[i]], color="black", lw=2.0)
gl = ax.gridlines(crs=ccrs.PlateCarree(), linewidth=1, color="black",
                  alpha=0.5, linestyle=":", draw_labels=True)
gl.xlocator = mticker.FixedLocator(range(95, 110, 3))
gl.ylocator = mticker.FixedLocator(range(20, 34, 2))
gl.top_labels = False
gl.right_labels = False
ax.set_aspect("auto")
# -----------------------------------------------------------------------------
fig.savefig(path_plot + "Envelope_convex_hull_cross.pdf", dpi=300,
            bbox_inches="tight")
plt.close(fig)

###############################################################################
# Check mean elevation increase of envelope topographies for subdomain
###############################################################################

# Load topographies
topo_env_mean = []
ds = xr.open_dataset(path_out + file_netcdf)
ds = ds.sel(lon=slice(96.0, 106.0), lat=slice(31.0, 22.0))
if np.any(np.isnan(ds["Elevation"].values)):
    raise ValueError("nan-value in data")
topo_mean = ds["Elevation"].values.mean()
for i in fac_curv:
    if np.any(np.isnan(ds["Elev_curv_%.1f" % i].values)):
        raise ValueError("nan-value in data")
    topo_env_mean.append(ds["Elev_curv_%.1f" % i].values.mean())
ds.close()

# Plot
fig = plt.figure(figsize=(6, 8))
data_plot = ((topo_env_mean / topo_mean) - 1.0) * 100.0
plt.plot(fac_curv, data_plot, color="dodgerblue", lw=1.5)
plt.scatter(fac_curv, data_plot, s=30, color="dodgerblue")
plt.xlabel("Curvature factor [-]")
plt.ylabel("Increase in mean elevation of envelope topography [%]", labelpad=8)
fig.savefig(path_plot + "Curv_fac_vs_mean_elev.png", dpi=300,
            bbox_inches="tight")
plt.close(fig)
