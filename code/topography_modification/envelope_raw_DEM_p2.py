# Description: Compute envelop topography part II (correct for isostatic
#              adjustment and fill terrain depressions)
#
# Source of data (beside MERIT DEM):
# - HydroBASINS
#   - https://www.hydrosheds.org/products/hydrobasins
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
from pyproj import CRS, Transformer
import time
from cmcrameri import cm
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import fiona
from shapely.geometry import shape
from descartes import PolygonPatch
import richdem as rd
from skimage.measure import label
from scipy import interpolate
from netCDF4 import Dataset

mpl.style.use("classic")

# Load required functions
sys.path.append("/Users/csteger/Downloads/BECCY/code/")
from auxiliary import truncate_colormap
from auxiliary import gridcoord
from auxiliary import spat_agg_2d, spat_agg_1d
sys.path.append("/Users/csteger/Downloads/BECCY/code/Isostasy/")
from isostasy_cy import deflection_lonlat

###############################################################################
# Settings
###############################################################################

# Paths
path_dem = "/Users/csteger/Dropbox/IAC/Data/DEMs/MERIT/Tiles/"
path_basins = "/Users/csteger/Dropbox/IAC/Data/Shapefiles/" \
              + "HydroBASINS/hybas_as_lev01-12_v1c/"
path_plot = "/Users/csteger/Dropbox/IAC/Plots/BECCY/Topo_envelop/"
path_in_out = "/Users/csteger/Dropbox/IAC/Data/Model/BECCY/"

# Miscellaneous
dom_proc = {"lat": slice(36.5 - 0.0001, 18.0),
            "lon": slice(89.0, 110.0 - 0.0001)}  # MERIT domain to process
agg_num_plot = 10  # spatial aggregation for plotting
agg_num_iso = 25  # spatial aggregation for isostatic adjustment (~2.3 km)
env_topo_sel = "Elev_curv_5.0"  # select envelope topography
# (old: "Elev_curv_8.0")

# Region with envelope topography
env_cen = (26.50, 100.80)  # centre of circle (latitude/longitude) [deg]
env_rad = 500.0 * 1000.0  # radius of circle [m]
env_bound = 100.0 * 1000.0  # boundary zone width [m]

# Constants for isostatic adjustment
rho_m = 3500.0  # mantle density  [kg m-3] (gFlex: 3300.0, old: 3400.0)
rho_nsr = 2300.0  # density of near-surface rock [kg m-3] (old: 2500.0)
rho_fill = 0.0  # infill material density (density of air: ~1.2) [kg m-3]
g = 9.78  # acceleration due to gravity [m s-2] (old: 9.81)
Te = 30000.0  # Elastic thickness [m] (old/gFlex: 35000.0)
E = 100E9  # Young's modulus [Pa] (old/gFlex: 65E9)
nu = 0.27  # Poisson's ratio [-] (old/gFlex: 0.25)

# Further constants
rad_earth = 6370997.0  # default PROJ sphere radius [m]
topo_fill_val = -32767  # MERIT 16-bit integer fill value

# -----------------------------------------------------------------------------
# References for constants
# -----------------------------------------------------------------------------

# - g: Chen et al. (2013)
# - Te: Chen et al. (2014), Wang et al. (2015)
# - E: Chen et al. (2013), Chen et al. (2014), Ou et al. (2021)
# - nu: Wang et al. (2017)

# Chen et al. (2013): Variations of the effective elastic thickness over China
#                     and surroundings and their relation to the lithosphere
#                     dynamics (http://dx.doi.org/10.1016/j.epsl.2012.12.022)
# Chen et al. (2014): Elastic thickness, mechanical anisotropy and deformation
#                     of the southeastern Tibetan Plateau
#                     (http://dx.doi.org/10.1016/j.tecto.2014.09.007)
# Ou et al. (2021): Contrasting exhumation histories and relief development
#                   within the Three Rivers Region (south-east Tibet)
#                   (https://doi.org/10.5194/se-12-563-2021)
# Wang et al. (2015): Flexural bending of southern Tibet in a retro foreland
#                     setting (https://doi.org/10.1038/srep12076)
# Wang et al. (2017): Crustal thickness and Poisson’s ratio in southwest China
#                     based on data from dense seismic arrays
#                     (https://doi.org/10.1002/2017JB013978)

###############################################################################
# Load/initialise data and make overview plot
###############################################################################

# Load raw MERIT DEM
files_dem = ("MERIT_N60-N30_E060-E090.nc", "MERIT_N30-N00_E060-E090.nc",
             "MERIT_N60-N30_E090-E120.nc", "MERIT_N30-N00_E090-E120.nc")
ds = xr.open_mfdataset([path_dem + i for i in files_dem])
ds = ds.sel(lon=dom_proc["lon"], lat=dom_proc["lat"])
topo_raw = ds["Elevation"].values  # 32-bit float
mask_ocean_raw = np.isnan(topo_raw)
lon_raw, lat_raw = ds["lon"].values, ds["lat"].values
ds.close()
if sum([topo_raw.shape[0] % agg_num_plot, topo_raw.shape[1] % agg_num_plot,
        topo_raw.shape[0] % agg_num_iso, topo_raw.shape[1] % agg_num_iso]) \
        != 0:
    raise ValueError("Dimensions lengths of DEM must be exact multiples "
                     + "of 'agg_num_plot' and 'agg_num_iso'")
print("Size of DEM: " + str(topo_raw.shape) + ", %.2f"
      % (topo_raw.nbytes / (10 ** 9)) + " GB")
print("Number of ocean grid cells: " + str(mask_ocean_raw.sum()))

# Load modified envelope MERIT data
file = "MERIT_envelope_topo.nc"
ds = xr.open_dataset(path_in_out + "env_topo_raw/" + file)
topo_env = ds[env_topo_sel].values.astype(np.float32)  # 32-bit float
# -> does not contain ocean grid cells but can contain nan-values close to
#    the edges
lon_env, lat_env = ds["lon"].values, ds["lat"].values
ind_lon = np.where(lon_raw == lon_env[0])[0][0]
ind_lat = np.where(lat_raw == lat_env[0])[0][0]
ds.close()

# Load river basins shapefiles
ds = fiona.open(path_basins + "hybas_as_lev03_v1c.shp")

# Spatial aggregation of data
data_agg_plot = {
    "lat": spat_agg_1d(lat_raw, agg_num_plot, operation="mean"),
    "lon": spat_agg_1d(lon_raw, agg_num_plot, operation="mean"),
    "topo_raw": spat_agg_2d(topo_raw, agg_num_plot, agg_num_plot,
                            operation="mean")
}
data_agg_iso = {
    "lat": spat_agg_1d(lat_raw, agg_num_iso, operation="mean"),
    "lon": spat_agg_1d(lon_raw, agg_num_iso, operation="mean"),
    "topo_ref": spat_agg_2d(topo_raw, agg_num_iso, agg_num_iso,
                            operation="mean")
}

# Colormap
levels_topo = np.arange(0., 5750.0, 250.0)
ticks_topo = np.arange(0., 6000.0, 1000.0)
cmap_topo = truncate_colormap(cm.bukavu, 0.55, 1.0)
norm_topo = mpl.colors.BoundaryNorm(levels_topo, ncolors=cmap_topo.N,
                                    extend="max")

# Map instances
crs_wgs84 = CRS.from_epsg(4326)
crs_aeqd = CRS.from_dict({"proj": "aeqd", "lat_0": env_cen[0],
                          "lon_0": env_cen[1], "datum": "WGS84", "units": "m"})

# Overview plot to define envelope region
fig = plt.figure(figsize=(11, 11))
ax = plt.axes(projection=ccrs.PlateCarree())
plt.pcolormesh(data_agg_plot["lon"], data_agg_plot["lat"],
               data_agg_plot["topo_raw"], cmap=cmap_topo,
               norm=norm_topo, shading="auto", zorder=1)
gl = ax.gridlines(crs=ccrs.PlateCarree(), linewidth=1, color="black",
                  alpha=0.5, linestyle=":", draw_labels=True)
gl.xlocator = mticker.FixedLocator(range(90, 110, 2))
gl.ylocator = mticker.FixedLocator(range(20, 40, 2))
gl.top_labels = False
gl.right_labels = False
ax.add_feature(cfeature.COASTLINE, zorder=3)
ax.add_feature(cfeature.RIVERS, zorder=3)
# -----------------------------------------------------------------------------
# for i in range(len(ds)):
for i in (16, 14, 11, 9, 7, 5):
    shp_geom = shape(ds[i]["geometry"])
    poly = PolygonPatch(shp_geom, facecolor="none", edgecolor="black",
                        zorder=2)
    ax.add_patch(poly)
    # cen_x, cen_y = shp_geom.centroid.xy
    # cen_x, cen_y = cen_x[0], cen_y[0]
    # plt.text(cen_x, cen_y, str(i), fontsize=12, fontweight="bold")
# -----------------------------------------------------------------------------
plt.scatter(env_cen[1], env_cen[0], s=50, color="red")
ls = ["-", "--"]
phi = np.linspace(0.0, 2.0 * np.pi, 1000)
count = 0
transformer = Transformer.from_crs(crs_aeqd, crs_wgs84, always_xy=True)
for i in (env_rad, env_rad + env_bound):
    x_circ, y_circ = i * np.cos(phi), i * np.sin(phi)
    lon_circ, lat_circ = transformer.transform(x_circ, y_circ)
    poly = plt.Polygon(list(zip(lon_circ, lat_circ)), facecolor="none",
                       edgecolor="red", linewidth=1.8, linestyle=ls[count])
    ax.add_patch(poly)
    count += 1
plt.axis([lon_raw.min(), lon_raw.max(), lat_raw.min(), lat_raw.max()])
fig.savefig(path_plot + "Envelope_region.png", dpi=300, bbox_inches="tight")
plt.close(fig)

# Create NetCDF file
file_netcdf = "MERIT_envelope_topo.nc"
ncfile = Dataset(filename=path_in_out + file_netcdf, mode="w")
ncfile.createDimension(dimname="lat", size=topo_raw.shape[0])
nc_lat = ncfile.createVariable(varname="lat", datatype="f", dimensions="lat")
nc_lat.units = "degrees_north"
nc_lat[:] = lat_raw
ncfile.createDimension(dimname="lon", size=topo_raw.shape[1])
nc_lat = ncfile.createVariable(varname="lon", datatype="f", dimensions="lon")
nc_lat.units = "degrees_east"
nc_lat[:] = lon_raw
nc_data = ncfile.createVariable(varname="Elevation", datatype="i2",
                                dimensions=("lat", "lon"),
                                fill_value=topo_fill_val)
nc_data.units = "m"
nc_data.long_name = "Unmodified topography"
nc_data[:] = np.where(mask_ocean_raw, topo_fill_val, topo_raw).astype(np.int16)
ncfile.close()

###############################################################################
# Compute envelope topography
###############################################################################

# Generate array with distance to centre
time_beg = time.time()
transformer = Transformer.from_crs(crs_wgs84, crs_aeqd, always_xy=True)
lon_2d, lat_2d = np.meshgrid(lon_env, lat_env)
x, y = transformer.transform(lon_2d, lat_2d)
dist = np.sqrt(x ** 2 + y ** 2)  # distance from centre [m]
del lon_2d, lat_2d, x, y
print("Compute distance to centre: %.1f" % (time.time() - time_beg) + " sec")

# # Test plot
# slic_agg = (slice(0, (topo_env.shape[0] // agg_num_plot) * agg_num_plot),
#             slice(0, (topo_env.shape[1] // agg_num_plot) * agg_num_plot))
# dist_agg = spat_agg_2d(dist[slic_agg], agg_num_plot, agg_num_plot,
#                        operation="mean")
# lon_agg = spat_agg_1d(lon_env[slic_agg[1]], agg_num_plot, operation="mean")
# lat_agg = spat_agg_1d(lat_env[slic_agg[0]], agg_num_plot, operation="mean")
# plt.figure()
# plt.pcolormesh(lon_agg, lat_agg, dist_agg, shading="auto")
# plt.colorbar()
# plt.contour(lon_agg, lat_agg, dist_agg,
#             levels=[env_rad, env_rad + env_bound], colors="black")

# Spatial aggregation of data
slic_env = (slice(ind_lat, ind_lat + topo_env.shape[0]),
            slice(ind_lon, ind_lon + topo_env.shape[1]))
mask_env = np.zeros(topo_raw.shape, dtype=bool)
mask_env[slic_env] = (dist <= env_rad)
data_agg_plot["mask_env"] = spat_agg_2d(mask_env.astype(np.float32),
                                        agg_num_plot, agg_num_plot,
                                        operation="mean")
del mask_env

# Create weights
mask = (dist > env_rad) & (dist < (env_rad + env_bound))
dist_scal = (dist[mask] - env_rad) / env_bound  # range [0.0, 1.0]
weights = np.ones_like(dist)
weights[dist <= env_rad] = 0.0
weights[mask] = (np.sin(dist_scal * np.pi - np.pi / 2.0) + 1.0) / 2.0
del mask, dist_scal, dist
if (np.any(weights[0, :] != 1.0) or np.any(weights[-1, :] != 1.0)
        or np.any(weights[:, 0] != 1.0) or np.any(weights[:, -1] != 1.0)):
    raise ValueError("envelope topography does not cover selected domain")

# # Test plot
# weights_agg = spat_agg_2d(weights[slic_agg], agg_num_plot, agg_num_plot,
#                           operation="mean")
# plt.figure()
# plt.pcolormesh(lon_agg, lat_agg, weights_agg, shading="auto")
# plt.colorbar()

# Deal with potential NaN-values in envelope topography
if np.any(np.isnan(topo_env)):
    print("NaN-values found in envelope topography (number: "
          + str(np.isnan(topo_env).sum()) + ")")
    if np.any(np.isnan(topo_env[(1.0 - weights) > 0.0])):
        raise ValueError("NaN-values are located in utilized region")
    else:
        print("Set NaN-values of envelope topography to 0.0")
        topo_env[np.isnan(topo_env)] = 0.0

# Merge topographies
mask_mod = np.zeros(topo_raw.shape, dtype=bool)  # modified grid cells
topo_raw[slic_env] = topo_raw[slic_env] * weights + topo_env * (1.0 - weights)
mask_mod[slic_env] = (weights < 1.0)
del weights, topo_env
if not np.all(np.isnan(topo_raw) == mask_ocean_raw):
    raise ValueError("inconsistency between topography and ocean mask")
data_agg_plot["topo_env"] = spat_agg_2d(topo_raw, agg_num_plot, agg_num_plot,
                                        operation="mean")
data_agg_plot["mask_mod_0"] = spat_agg_2d(mask_mod.astype(np.float32),
                                          agg_num_plot, agg_num_plot,
                                          operation="mean")
data_agg_iso["topo_mod"] = spat_agg_2d(topo_raw, agg_num_iso, agg_num_iso,
                                       operation="mean")

# Append to NetCDF file
ncfile = Dataset(filename=path_in_out + file_netcdf, mode="a")
nc_data = ncfile.createVariable(varname="Elevation_env", datatype="i2",
                                dimensions=("lat", "lon"),
                                fill_value=topo_fill_val)
nc_data.units = "m"
nc_data.long_name = "Envelope topography"
nc_data[:] = np.where(mask_ocean_raw, topo_fill_val, topo_raw).astype(np.int16)
ncfile.close()

###############################################################################
# Adjust for isostasy and fill depressions iteratively
###############################################################################

# -----------------------------------------------------------------------------
# Constants for isostatic adjustment
# -----------------------------------------------------------------------------

# Compute parameters
d = E * Te ** 3 / (12.0 * (1.0 - nu ** 2))  # Flexural rigidity [m2 kg s-2]
alpha = (d / ((rho_m - rho_fill) * g)) ** 0.25  # 2D flexural parameter [m]

# Compute surface area of spatially aggregated grid cells
lon_grid, lat_grid = gridcoord(np.deg2rad(data_agg_iso["lon"]),
                               np.deg2rad(data_agg_iso["lat"]))  # [rad]
d_lon = np.diff(np.deg2rad(data_agg_iso["lon"])).mean()  # [rad]
gc_area = rad_earth ** 2 * d_lon * \
          (np.sin(lat_grid[:-1]) - np.sin(lat_grid[1:]))  # [m2]
gc_area = np.repeat(gc_area[:, np.newaxis], len(data_agg_iso["lon"]), axis=1)

# -----------------------------------------------------------------------------
# Iteration
# -----------------------------------------------------------------------------


def print_save(txt, file_handle):
    print(txt)
    file_handle.write(txt + "\n")


iter_con = True
iter_num = 1
file_txt = open(path_in_out + file_netcdf[:-3] + ".txt", mode="w")
while iter_con:

    print_save((" Iteration " + str(iter_num) + " ").center(79, "#"), file_txt)

    # -------------------------------------------------------------------------
    # Isostatic adjustment
    # -------------------------------------------------------------------------
    print_save(" Isostatic adjustment ".center(79, "-"), file_txt)

    # Compute point loadings (spatially aggregated)
    elev_diff = (data_agg_iso["topo_mod"] - data_agg_iso["topo_ref"])
    print_save("Minimal elevation difference: %.2f" % np.nanmin(elev_diff)
               + " m", file_txt)
    elev_diff[np.isnan(elev_diff)] = 0.0
    elev_diff = elev_diff.clip(min=0.0)
    q = (elev_diff * gc_area) * rho_nsr * g  # [N]

    # Compute deflection
    time_beg = time.time()
    w = deflection_lonlat(np.deg2rad(data_agg_iso["lon"]),
                          np.deg2rad(data_agg_iso["lat"]),
                          q, alpha, d, "fast")
    print("Isostatic adjustment computation: %.1f" % (time.time() - time_beg)
          + " sec")
    txt = "Maximal absolute deflections at boundaries (spatial agg.): " \
          + "%.2f" % np.abs(w[:, -1]).max() + " m (east), " \
          + "%.2f" % np.abs(w[-1, :]).max() + " m (south), " \
          + "%.2f" % np.abs(w[:, 0]).max() + " m (west), " \
          + "%.2f" % np.abs(w[0, :]).max() + " m (north)"
    print_save(txt, file_txt)

    # Interpolate from spatially aggregated to raw grid
    f_ip = interpolate.RectBivariateSpline(data_agg_iso["lat"][::-1],
                                           data_agg_iso["lon"],
                                           np.flipud(w),  kx=1, ky=1)
    w_ip = np.flipud(f_ip(lat_raw[::-1], lon_raw).astype(np.float32))
    txt = "Maximal absolute deflections at boundaries: " \
          + "%.2f" % np.abs(w_ip[:, -1]).max() + " m (east), " \
          + "%.2f" % np.abs(w_ip[-1, :]).max() + " m (south), " \
          + "%.2f" % np.abs(w_ip[:, 0]).max() + " m (west), "\
          + "%.2f" % np.abs(w_ip[0, :]).max() + " m (north)"
    print_save(txt, file_txt)

    # Adjust topography
    topo_raw += w_ip
    data_agg_plot["deflection_" + str(iter_num)] \
        = spat_agg_2d(w_ip, agg_num_plot, agg_num_plot, operation="mean")
    print_save("Spatial mean adjustment: %.1f" % w_ip.mean() + " m", file_txt)

    # Update grid cell modification mask
    mask_mod = (mask_mod | (np.abs(w_ip) > 10.0))
    data_agg_plot["mask_mod_" + str(iter_num)] = \
        spat_agg_2d(mask_mod.astype(np.float32), agg_num_plot, agg_num_plot,
                    operation="mean")
    del w_ip

    # -------------------------------------------------------------------------
    # Fill depressions
    # -------------------------------------------------------------------------
    print_save(" Fill depressions ".center(79, "-"), file_txt)

    # Find depressions
    no_data = -9999.0
    topo_raw[mask_ocean_raw] = no_data
    time_beg = time.time()
    topo_fill = np.array(rd.FillDepressions(rd.rdarray(topo_raw,
                                                       no_data=no_data),
                                            in_place=False))
    # -> outflow possible both at boundaries and "no data regions"
    # -> "no_data=np.nan" produces erroneous output
    topo_raw[mask_ocean_raw] = np.nan
    topo_fill[mask_ocean_raw] = np.nan
    print("Depression filling: %.1f" % (time.time() - time_beg) + " sec")

    # Select relevant depressions
    time_beg = time.time()
    mask_dep = ((topo_fill - topo_raw) > 0.0).astype(np.int8)
    labels, num = label(mask_dep, return_num=True)
    labels_sel = np.unique(labels[mask_mod])
    if labels_sel[0] == 0:
        labels_sel = labels_sel[1:]
    print_save("Number of labels (total/relevant): " + str(num) + "/"
               + str(len(labels_sel)), file_txt)
    mask_dep_sel = np.isin(labels, labels_sel)
    del labels
    print("Select relevant depressions: %.1f" % (time.time() - time_beg)
          + " sec")

    # Adjust topography
    dep_fill = np.zeros(topo_raw.shape, dtype=np.float32)
    dep_fill[mask_dep_sel] = (topo_fill - topo_raw)[mask_dep_sel]
    print_save("Spatial mean adjustment: %.1f" % dep_fill.mean() + " m",
               file_txt)
    dep_depth_max = dep_fill.max()
    print_save("Maximal depression depth: %.1f" % dep_depth_max + " m",
               file_txt)
    data_agg_plot["dep_fill_" + str(iter_num)] \
        = spat_agg_2d(dep_fill, agg_num_plot, agg_num_plot, operation="mean")
    data_agg_iso["topo_ref"] = spat_agg_2d(topo_raw, agg_num_iso, agg_num_iso,
                                           operation="mean")
    topo_raw += dep_fill
    data_agg_iso["topo_mod"] = spat_agg_2d(topo_raw, agg_num_iso, agg_num_iso,
                                           operation="mean")
    del topo_fill, dep_fill

    # Check iteration break criteria
    iter_num += 1
    if dep_depth_max < 30.0:
        iter_con = False

file_txt.close()

###############################################################################
# Save topography to NetCDF
###############################################################################

if not np.all(np.isnan(topo_raw) == mask_ocean_raw):
    raise ValueError("inconsistency between topography and ocean mask")

# Load unmodified topography
ds = xr.open_dataset(path_in_out + file_netcdf)
topo_unmod = ds["Elevation"].values
topo_env = ds["Elevation_env"].values
ds.close()

# Append to NetCDF file
ncfile = Dataset(filename=path_in_out + file_netcdf, mode="a")
# -----------------------------------------------------------------------------
nc_data = ncfile.createVariable(varname="Elevation_final", datatype="i2",
                                dimensions=("lat", "lon"),
                                fill_value=topo_fill_val)
nc_data.units = "m"
nc_data.long_name = "Final topography (envelope, isostatic adjustment " \
                    + " and depression filling)"
nc_data[:] = np.where(mask_ocean_raw, topo_fill_val, topo_raw).astype(np.int16)
# -----------------------------------------------------------------------------
nc_data = ncfile.createVariable(varname="Elev_diff_1", datatype="i2",
                                dimensions=("lat", "lon"),
                                fill_value=topo_fill_val)
nc_data.units = "m"
nc_data.long_name = "Difference between final and unmodified topography"
nc_data[:] = np.where(mask_ocean_raw, topo_fill_val, topo_raw - topo_unmod) \
    .astype(np.int16)
# -----------------------------------------------------------------------------
nc_data = ncfile.createVariable(varname="Elev_diff_2", datatype="i2",
                                dimensions=("lat", "lon"),
                                fill_value=topo_fill_val)
nc_data.units = "m"
nc_data.long_name = "Difference between final and envelope topography"
nc_data[:] = np.where(mask_ocean_raw, topo_fill_val, topo_raw - topo_env) \
    .astype(np.int16)
# -----------------------------------------------------------------------------
ncfile.close()

del topo_unmod, topo_env

# Save final topography in MERIT tiles
if not os.path.exists(path_in_out + "MERIT_tiles/"):
    os.makedirs(path_in_out + "MERIT_tiles/")
mask_assign = np.zeros(topo_raw.shape, dtype=np.uint8)
for i in files_dem:
    ds = xr.open_dataset(path_dem + i)
    ind_raw_lon = np.where(np.in1d(lon_raw, ds["lon"].values))[0][[0, -1]]
    ind_raw_lat = np.where(np.in1d(lat_raw, ds["lat"].values))[0][[0, -1]]
    ind_tile_lon = np.where(np.in1d(ds["lon"].values, lon_raw))[0][[0, -1]]
    ind_tile_lat = np.where(np.in1d(ds["lat"].values, lat_raw))[0][[0, -1]]
    slic_tile = (slice(ind_tile_lat[0], ind_tile_lat[-1] + 1),
                 slice(ind_tile_lon[0], ind_tile_lon[-1] + 1))
    slic_raw = (slice(ind_raw_lat[0], ind_raw_lat[-1] + 1),
                slice(ind_raw_lon[0], ind_raw_lon[-1] + 1))
    ds["Elevation"][slic_tile] = topo_raw[slic_raw]
    mask_assign[slic_raw] += 1
    ds.to_netcdf(path_in_out + "MERIT_tiles/" + i, format="NETCDF4",
                 encoding={"lat": {"_FillValue": None},
                           "lon": {"_FillValue": None}})
    print("MERIT tile " + i + " modified")
if not np.all(mask_assign):
    raise ValueError("error in assigning modified elevation data "
                     + "to MERIT tiles")

# Spatial aggregation of data
data_agg_plot["topo_final"] = spat_agg_2d(topo_raw, agg_num_plot, agg_num_plot,
                                          operation="mean")
del topo_raw

###############################################################################
# Compute changes in mountain volumes
###############################################################################

# Compute surface area of spatially aggregated grid cells
lon_grid, lat_grid = gridcoord(np.deg2rad(data_agg_plot["lon"]),
                               np.deg2rad(data_agg_plot["lat"]))  # [rad]
d_lon = np.diff(np.deg2rad(data_agg_plot["lon"])).mean()  # [rad]
gc_area = rad_earth ** 2 * d_lon * \
          (np.sin(lat_grid[:-1]) - np.sin(lat_grid[1:]))  # [m2]
gc_area = np.repeat(gc_area[:, np.newaxis], len(data_agg_plot["lon"]), axis=1)

# Sum spatially aggregated adjustments
deflect_all = np.concatenate([data_agg_plot["deflection_" + str(i)]
                              [:, :, np.newaxis] for i in range(1, iter_num)],
                             axis=2).sum(axis=2)
dep_fill_all = np.concatenate([data_agg_plot["dep_fill_" + str(i)]
                               [:, :, np.newaxis] for i in range(1, iter_num)],
                              axis=2).sum(axis=2)
mask_mod_agg = (np.abs(deflect_all) > 10.0) | (dep_fill_all > 0.0)

# Test plot
plt.figure()
plt.pcolormesh(data_agg_plot["lon"], data_agg_plot["lat"],
               # np.abs(deflect_all) > 10.0,
               # dep_fill_all > 0.0,
               # mask_mod_agg,
               # data_agg_plot["mask_mod_0"] > 0.5,
               data_agg_plot["mask_env"] > 0.5,
               shading="auto")
plt.colorbar()

# Compute changes in volume
masks = {"gc_mod": mask_mod_agg,
         "gc_env_bound": data_agg_plot["mask_mod_0"] > 0.5,
         "gc_env": data_agg_plot["mask_env"] > 0.5}
file_txt = open(path_in_out + "Mountain_volume_increase.txt", mode="w")
file_txt.write("Mountain volume increases:\n")
for i in list(masks.keys()):
    vol_raw = (data_agg_plot["topo_raw"] * gc_area)[masks[i]].sum()  # [m3]
    vol_env = (data_agg_plot["topo_env"] * gc_area)[masks[i]].sum()  # [m3]
    vol_fin = (data_agg_plot["topo_final"] * gc_area)[masks[i]].sum()  # [m3]
    txt = (" Domain: " + i + " ").center(60, "-")
    file_txt.write(txt + "\n")
    txt = "Raw -> envelope: % .2f" % ((vol_env / vol_raw - 1.0) * 100.0) + " %"
    file_txt.write(txt + "\n")
    txt = "Raw -> final: % .2f" % ((vol_fin / vol_raw - 1.0) * 100.0) + " %"
    file_txt.write(txt + "\n")
file_txt.close()

###############################################################################
# Plot
###############################################################################

dom_ext_plot = [92.0, 107.5, 20.5, 34.0]

# Plot
fig = plt.figure(figsize=(10.5, 16.0))
gs = gridspec.GridSpec(4, 2, left=0.1, bottom=0.1, right=0.9,
                       top=0.9, hspace=0.09, wspace=0.05,
                       height_ratios=[1.0, 1.0, 1.0, 0.05],
                       width_ratios=[1.0, 1.0])
# -----------------------------------------------------------------------------
ax = plt.subplot(gs[0, 0], projection=ccrs.PlateCarree())
plt.pcolormesh(data_agg_plot["lon"], data_agg_plot["lat"],
               data_agg_plot["topo_raw"], cmap=cmap_topo, norm=norm_topo,
               shading="auto")
ax.set_aspect("auto")
ax.set_extent(dom_ext_plot, crs=ccrs.PlateCarree())
plt.title("Present-day topography", fontsize=10,
          fontweight="bold", y=1.002)
plt.scatter(env_cen[1], env_cen[0], s=50, color="red")
ls = ["-", "--"]
count = 0
transformer = Transformer.from_crs(crs_aeqd, crs_wgs84, always_xy=True)
for i in (env_rad, env_rad + env_bound):
    x_circ, y_circ = i * np.cos(phi), i * np.sin(phi)
    lon_circ, lat_circ = transformer.transform(x_circ, y_circ)
    poly = plt.Polygon(list(zip(lon_circ, lat_circ)), facecolor="none",
                       edgecolor="red", linewidth=1.2, linestyle=ls[count])
    ax.add_patch(poly)
    count += 1
plt.text(0.05, 0.95, "(a)", fontsize=10, fontweight="bold",
         horizontalalignment="center", verticalalignment="center",
         transform=ax.transAxes)
# -----------------------------------------------------------------------------
ax = plt.subplot(gs[1, 0], projection=ccrs.PlateCarree())
plt.pcolormesh(data_agg_plot["lon"], data_agg_plot["lat"],
               data_agg_plot["topo_env"], cmap=cmap_topo, norm=norm_topo,
               shading="auto")
ax.set_aspect("auto")
ax.set_extent(dom_ext_plot, crs=ccrs.PlateCarree())
plt.title("Embedded raw envelope topography", fontsize=10,
          fontweight="bold", y=1.002)
plt.text(0.05, 0.95, "(c)", fontsize=10, fontweight="bold",
         horizontalalignment="center", verticalalignment="center",
         transform=ax.transAxes)
# -----------------------------------------------------------------------------
ax = plt.subplot(gs[2, 0], projection=ccrs.PlateCarree())
plt.pcolormesh(data_agg_plot["lon"], data_agg_plot["lat"],
               data_agg_plot["topo_final"], cmap=cmap_topo, norm=norm_topo,
               shading="auto")
ax.set_aspect("auto")
ax.set_extent(dom_ext_plot, crs=ccrs.PlateCarree())
plt.title("Envelope topography", fontsize=10, fontweight="bold", y=1.002)
plt.text(0.05, 0.95, "(e)", fontsize=10, fontweight="bold",
         horizontalalignment="center", verticalalignment="center",
         transform=ax.transAxes)
# -----------------------------------------------------------------------------
ax = plt.subplot(gs[3, 0])
cbar = mpl.colorbar.ColorbarBase(ax, cmap=cmap_topo, norm=norm_topo,
                                 ticks=ticks_topo, orientation="horizontal")
cbar.ax.tick_params(labelsize=10)
plt.xlabel("Elevation [m]", fontsize=10)
# -----------------------------------------------------------------------------
levels_diff = np.arange(-1500.0, 1700.0, 200.0)
cmap_diff = cm.roma
norm_diff = mpl.colors.BoundaryNorm(levels_diff, ncolors=cmap_diff.N,
                                    extend="both")
ax = plt.subplot(gs[0, 1], projection=ccrs.PlateCarree())
plt.pcolormesh(data_agg_plot["lon"], data_agg_plot["lat"],
               (data_agg_plot["topo_env"] - data_agg_plot["topo_raw"]),
               cmap=cmap_diff, norm=norm_diff, shading="auto")
ax.set_aspect("auto")
ax.set_extent(dom_ext_plot, crs=ccrs.PlateCarree())
plt.title("Difference (embedded raw envelope - present-day)",
          fontsize=10,  fontweight="bold", y=1.002)
plt.text(0.05, 0.95, "(b)", fontsize=10, fontweight="bold",
         horizontalalignment="center", verticalalignment="center",
         transform=ax.transAxes)
# -----------------------------------------------------------------------------
ax = plt.subplot(gs[1, 1], projection=ccrs.PlateCarree())
plt.pcolormesh(data_agg_plot["lon"], data_agg_plot["lat"],
               (data_agg_plot["topo_final"] - data_agg_plot["topo_env"]),
               cmap=cmap_diff, norm=norm_diff, shading="auto")
ax.set_aspect("auto")
ax.set_extent(dom_ext_plot, crs=ccrs.PlateCarree())
plt.title("Difference (envelope - embedded raw envelope)", fontsize=10,
          fontweight="bold", y=1.002)
plt.text(0.05, 0.95, "(d)", fontsize=10, fontweight="bold",
         horizontalalignment="center", verticalalignment="center",
         transform=ax.transAxes)
# -----------------------------------------------------------------------------
ax = plt.subplot(gs[2, 1], projection=ccrs.PlateCarree())
plt.pcolormesh(data_agg_plot["lon"], data_agg_plot["lat"],
               (data_agg_plot["topo_final"] - data_agg_plot["topo_raw"]),
               cmap=cmap_diff, norm=norm_diff, shading="auto")
ax.set_aspect("auto")
ax.set_extent(dom_ext_plot, crs=ccrs.PlateCarree())
plt.title("Difference (envelope - present-day)", fontsize=10,
          fontweight="bold", y=1.002)
plt.text(0.05, 0.95, "(f)", fontsize=10, fontweight="bold",
         horizontalalignment="center", verticalalignment="center",
         transform=ax.transAxes)
# -----------------------------------------------------------------------------
ax = plt.subplot(gs[3, 1])
cbar = mpl.colorbar.ColorbarBase(ax, cmap=cmap_diff, norm=norm_diff,
                                 orientation="horizontal")
cbar.ax.tick_params(labelsize=10)
plt.xlabel("Elevation difference [m]", fontsize=10)
# -----------------------------------------------------------------------------
fig.savefig(path_plot + "Envelope_steps.png", dpi=300, bbox_inches="tight")
plt.close(fig)
