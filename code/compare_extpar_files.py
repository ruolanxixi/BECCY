# Description: Compare EXTPAR files (cross-sections and mountain volumes)
#
# Authors: Christian R. Steger, IAC ETH Zurich

# Load modules
import xarray as xr
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from pyproj import CRS, Transformer

mpl.style.use("classic")

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

# Plot cross-sections
fig = plt.figure(figsize=(16.5, 10.5))
gs = gridspec.GridSpec(len(rlat_sel), 1, left=0.1, bottom=0.1, right=0.9,
                       top=0.9, hspace=0.11, wspace=0.05)
count = 0
for i in rlat_sel:
    ax = plt.subplot(gs[count])
    for j in list(data.keys()):
        ind = np.argmin(np.abs(data[j]["rlat"] - i))
        res, exp = j.split("_")

        plt.plot(data[j]["rlon"], data[j]["topo"][ind, :], lw=1.5,
                 color=cols[res], ls=ls[exp])
        if "unmod" in j:
            plt.scatter(data[j]["rlon"], data[j]["topo"][ind, :], s=20,
                        color=cols[res])
    plt.xlim(rlon_rang)
    count += 1

###############################################################################
# Compute mountain volumes of different experiments
###############################################################################
