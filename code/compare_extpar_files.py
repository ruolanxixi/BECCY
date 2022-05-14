# Description: Compare EXTPAR files (cross-sections and mountain volumes)
#
# Authors: Christian R. Steger, IAC ETH Zurich

# Load modules
import sys
import xarray as xr
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from pyproj import CRS, Transformer
from cmcrameri import cm

mpl.style.use("classic")

# Load required functions
sys.path.append("/Users/csteger/Downloads/BECCY/code/")
from auxiliary import truncate_colormap
from auxiliary import gridcoord

###############################################################################
# Settings
###############################################################################

# EXTPAR files
path = "/Users/csteger/Dropbox/IAC/Data/Model/BECCY/EXTPAR_files/"
files = {"12km_unmod":  "extpar_EAS_ext_12km_merit_unmod_topo.nc",
         "12km_env":    "extpar_EAS_ext_12km_merit_env_topo.nc",
         "4.4km_unmod": "extpar_BECCY_4.4km_merit_unmod_topo.nc",
         "4.4km_env":   "extpar_BECCY_4.4km_merit_env_topo.nc",
         "2.2km_unmod": "extpar_BECCY_2.2km_merit_unmod_topo.nc",
         "2.2km_env":   "extpar_BECCY_2.2km_merit_env_topo.nc"}

# Latitude transect values
rlat_trans = [-2.27, -1.27, -0.07]  # [degree]
rlon_rang = [-18.0, -9.0]
file_ref = "12km_unmod"  # transects selected according to this file

# Region with envelope topography
env_cen = (26.50, 100.80)  # centre of circle (latitude/longitude) [deg]
env_rad = 500.0 * 1000.0  # radius of circle [m]
env_bound = 100.0 * 1000.0  # boundary zone width [m]

###############################################################################
# Load and process data
###############################################################################

# Load data
data = {}
for i in list(files.keys()):
    data[i] = {}
    ds = xr.open_dataset(path + files[i])
    data[i]["rlon"] = ds["rlon"].values  # [degree]
    data[i]["rlat"] = ds["rlat"].values  # [degree]
    data[i]["lon"] = ds["lon"].values  # [degree]
    data[i]["lat"] = ds["lat"].values  # [degree]
    data[i]["topo"] = ds["HSURF"].values  # [m]
    data[i]["rotated_pole"] = ds["rotated_pole"].attrs
    ds.close()

# Ensure that all input data is on same rotated grid
if any([data[i]["rotated_pole"] != data[list(data.keys())[0]]["rotated_pole"]
        for i in data.keys()]):
    raise ValueError("input data does not share the same rotated grid")

# Map instances
crs_wgs84 = CRS.from_epsg(4326)
crs_aeqd = CRS.from_dict({"proj": "aeqd", "lat_0": env_cen[0],
                          "lon_0": env_cen[1], "datum": "WGS84", "units": "m"})

###############################################################################
# Plot cross-sections along rotated latitude
###############################################################################

# Compute selected rotated latitudes
ind_ref = [np.argmin(np.abs(data[file_ref]["rlat"] - i)) for i in rlat_trans]
rlat_sel = data[file_ref]["rlat"][ind_ref]  # [degree]

# Plot settings
cols = {"12km": "blue", "4.4km": "red", "2.2km": "green"}
ls = {"unmod": "-", "env": "--"}

# Colormap
levels_topo = np.arange(0., 5750.0, 250.0)
ticks_topo = np.arange(0., 6000.0, 1000.0)
cmap_topo = truncate_colormap(cm.bukavu, 0.55, 1.0)
norm_topo = mpl.colors.BoundaryNorm(levels_topo, ncolors=cmap_topo.N,
                                    extend="max")

# Overview map plot
key = "2.2km_unmod"
plt.figure()
plt.pcolormesh(data[key]["rlon"], data[key]["rlat"], data[key]["topo"],
               cmap=cmap_topo, norm=norm_topo, shading="auto")
plt.xlim(np.array(rlon_rang) + np.array([-1.5, 1.5]))
plt.ylim([rlat_sel.min() - 2.5, rlat_sel.max() + 2.5])
for i in rlat_sel:
    plt.plot(rlon_rang, [i, i], lw=1.5, color="black")
    plt.scatter(rlon_rang, [i, i], s=30, color="black")
plt.xlabel("Rotated longitude [degree]")
plt.ylabel("Rotated latitude [degree]")

# Plot cross-sections
fig = plt.figure(figsize=(16.5, 10.5))
gs = gridspec.GridSpec(len(rlat_sel), 1, left=0.1, bottom=0.1, right=0.9,
                       top=0.9, hspace=0.11, wspace=0.05)
count = 0
for i in rlat_sel[::-1]:
    ax = plt.subplot(gs[count])
    for j in list(data.keys()):
        ind = np.argmin(np.abs(data[j]["rlat"] - i))
        res, exp = j.split("_")
        if exp == "unmod":
            plt.plot(data[j]["rlon"], data[j]["topo"][ind, :], lw=1.5,
                     color=cols[res], ls=ls[exp], label=res)
        else:
            plt.plot(data[j]["rlon"], data[j]["topo"][ind, :], lw=1.5,
                     color=cols[res], ls=ls[exp])
        if "unmod" in j:
            plt.scatter(data[j]["rlon"], data[j]["topo"][ind, :], s=20,
                        color=cols[res])
    t = plt.text(0.09, 0.90, "Rot. latitude: %.2f" % i + " deg",
                 fontsize=11,
                 horizontalalignment="center", verticalalignment="center",
                 transform=ax.transAxes, fontweight="bold")
    plt.xlim(rlon_rang)
    plt.xlabel("Rotated longitude [degree]")
    plt.ylabel("Elevation [m]")
    if count == 0:
        plt.legend(fontsize=12, frameon=False)
    count += 1

###############################################################################
# Compute mountain volumes of different experiments
###############################################################################

# Map instances
crs_wgs84 = CRS.from_epsg(4326)
crs_aeqd = CRS.from_dict({"proj": "aeqd", "lat_0": env_cen[0],
                          "lon_0": env_cen[1], "datum": "WGS84", "units": "m"})

# Loop through resolutions
for i in ("12km", "4.4km", "2.2km"):

    # Transform coordinates
    transformer = Transformer.from_crs(crs_wgs84, crs_aeqd, always_xy=True)
    x, y = transformer.transform(data[i + "_unmod"]["lon"],
                                 data[i + "_unmod"]["lat"])
    dist = np.sqrt(x ** 2 + y ** 2)  # distance from centre [m]
    mask = (dist <= env_rad)
    # mask = (dist <= (env_rad + env_bound))

    # # Test plot
    # plt.figure()
    # plt.pcolormesh(x, y, mask)
    # plt.colorbar()
    # plt.figure()
    # plt.pcolormesh(x, y, data[i + "_unmod"]["topo"], cmap=cmap_topo,
    #                norm=norm_topo)
    # plt.colorbar()

    # Compute grid cell area
    rad_earth = 6370997.0  # default PROJ sphere radius [m]
    lon_grid, lat_grid = gridcoord(np.deg2rad(data[i + "_unmod"]["rlon"]),
                                   np.deg2rad(data[i + "_unmod"]["rlat"]))
    d_lon = np.diff(np.deg2rad(data[i + "_unmod"]["rlon"])).mean()  # [rad]
    gc_area = rad_earth ** 2 * d_lon * \
        (np.sin(lat_grid[1:]) - np.sin(lat_grid[:-1]))  # [m2]
    gc_area = np.repeat(gc_area[:, np.newaxis],
                        len(data[i + "_unmod"]["rlon"]), axis=1)

    # Compute relative mountain volume change
    vol_unmod = data[i + "_unmod"]["topo"] * gc_area  # [m3]
    vol_env = data[i + "_env"]["topo"] * gc_area  # [m3]
    rel_chang = (vol_env[mask].sum() / vol_unmod[mask].sum() - 1.0) * 100.0
    print("Mountain volume change (" + i + "): %.2f" % rel_chang + " %")
