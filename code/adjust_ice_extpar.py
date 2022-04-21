# Description: Adjust grid cells covered with permanent snow/ice according to
#              modified topography (also adapt linked fields)
#
# Authors: Ruolan Xiang, Christian R. Steger, IAC ETH Zurich

# Load modules
import xarray as xr
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from pyproj import CRS, Transformer
from scipy.spatial import cKDTree
# from netCDF4 import Dataset

mpl.style.use("classic")

###############################################################################
# Settings
###############################################################################

# Elevation threshold for permanent snow/ice coverage
elev_thresh = [3810.0, 4640.0, 5490.0]  # percentiles (5.0, 50.0, 95.0) [m]

# Directions to adjust
adj_ice2soil = True
adj_soil2ice = False

# EXTPAR files
# file_ref = "/Users/csteger/Desktop/extpar_BECCY_4.4km_merit_unmod_topo.nc"
# file_mod = "/Users/csteger/Desktop/extpar_BECCY_4.4km_merit_reduced_topo.nc"
file_ref = "/Users/csteger/Desktop/extpar_EAS_ext_12km_merit_unmod_topo.nc"
file_mod = "/Users/csteger/Desktop/extpar_EAS_ext_12km_merit_reduced_topo.nc"

# Search radius for (nearest) neighbour grid cells
# rad_search = 60.0 * 1000.0  # [m] (4.4 km EXTPAR file)
rad_search = 90.0 * 1000.0  # [m] (12 km EXTPAR file)

###############################################################################
# Adjust EXTPAR file
###############################################################################

# Find variable fields that differ between EXTPAR files
print(" Variable fields that differ " .center(79, "-"))
ds = xr.open_dataset(file_ref)
ds_mod = xr.open_dataset(file_mod)
for i in list(ds.variables):
    if not np.all(ds[i].values == ds_mod[i].values):
        print((i + ":").ljust(12) + ds[i].long_name)
ds.close()
ds_mod.close()
print("-" * 79)

# Load data from EXTPAR files
ds = xr.open_dataset(file_ref)
topo_unmod = ds["HSURF"].values
fr_land = ds["FR_LAND"].values
soiltyp = ds["SOILTYP"].values
ice = ds["ICE"].values
lon = ds["lon"].values  # [degree]
lat = ds["lat"].values  # [degree]
ds.close()
ds = xr.open_dataset(file_mod)
topo_mod = ds["HSURF"].values
ds.close()

# Find adjustable grid cells
mask_mod = (np.abs(topo_unmod - topo_mod) > 0.001)
# mask_ice2soil = (mask_mod & (soiltyp == 1.0) & (topo_mod < elev_thresh[0]))
mask_ice2soil = (mask_mod & (ice > 0.0) & (topo_mod < elev_thresh[0]))
print("Number of grid cells (ice2soil): " + str(mask_ice2soil.sum()))
mask_soil2ice = (mask_mod & (soiltyp != 1.0) & (topo_mod > elev_thresh[2]))
print("Number of grid cells (soil2ice): " + str(mask_soil2ice.sum()))

# Statistics for grid cells that have to be adjusted
print(" Statistics of adjusted grid cells with water fraction "
      .center(79, "-"))
mask = np.zeros(topo_unmod.shape, dtype=bool)
if adj_ice2soil:
    mask[mask_ice2soil] = True
if adj_soil2ice:
    mask[mask_soil2ice] = True
print("Minimal FR_LAND value: %.3f" % fr_land[mask].min())
num = (fr_land[mask] != 1.0).sum()
print("Number of grid cells for which FR_LAND != 1.0: " + str(num)
      + " (%.2f" % (num / mask.sum() * 100.0) + " %)")
print("Number of grid cells for which FR_LAND < 0.99: "
      + str((fr_land[mask] < 0.99).sum()))
print("-" * 79)

# -----------------------------------------------------------------------------
# Notes for adjusting values
# -----------------------------------------------------------------------------

# Not adjusted -> depend on elevation (9)
# - HSURF, FIS
# - S_ORO, Z0, T_CL
# - SSO_STDH, SSO_THETA, SSO_GAMMA, SSO_SIGMA

# Not adjusted -> fixed (2)
# - lon, lat

# Not adjusted -> atmospheric column (5)
# - AER_BC12, AER_DUST12, AER_ORG12, AER_SO412, AER_SS12 (size: 12)

# Values adjusted -> depend on surface/soil (23)
# - ICE
# - PLCOV_MN, PLCOV_MX, LAI_MN, LAI_MX
# - EMIS_RAD
# - RSMIN
# - URBAN
# - FOR_D, FOR_E
# - SKC
# - ROOTDP
# - NDVI_MAX
# - SOILTYP
# - LU_CLASS_FRACTION (size: 23)
# - ALB_DIF12, ALNID12, ALUVD12, NDVI, NDVI_MRAT (size: 12)
# - FR_LAND, FR_LAKE, DEPTH_LK

# Miscellaneous
# - LU_CLASS_FRACTION -> add up to 1.0
# - also replace water/lake properties to keep grid data consistent

# -----------------------------------------------------------------------------

# Compute ECEF coordinates and construct tree
crs_ecef = CRS.from_dict({"proj": "geocent", "ellps": "sphere"})
crs_latlon = CRS.from_dict({"proj": "latlong", "ellps": "sphere"})
trans = Transformer.from_crs(crs_latlon, crs_ecef, always_xy=True)
x_ecef, y_ecef, z_ecef = trans.transform(lon, lat, np.zeros_like(lon))
pts_gc = np.vstack((x_ecef.ravel(), y_ecef.ravel(), z_ecef.ravel())) \
    .transpose()
tree = cKDTree(pts_gc)

# -----------------------------------------------------------------------------
# Derive grid cell replacement indices for ice -> soil
# -----------------------------------------------------------------------------
if adj_ice2soil:

    # Indices for replacement (source -> target)
    num = mask_ice2soil.sum()
    ind_ta = np.where(mask_ice2soil)
    ind_su = (np.empty(num, dtype=np.int32), np.empty(num, dtype=np.int32))
    dist_all = np.empty(num, dtype=np.float32)
    frac_diff_all = np.empty(num, dtype=np.float32)
    for i in range(num):

        ind_2d_ta = (ind_ta[0][i], ind_ta[1][i])
        ind_lin = tree.query_ball_point([x_ecef[ind_2d_ta],
                                         y_ecef[ind_2d_ta],
                                         z_ecef[ind_2d_ta]],
                                        r=rad_search)  # linear indices
        dist = np.sqrt((x_ecef.ravel()[ind_lin] - x_ecef[ind_2d_ta]) ** 2
                       + (y_ecef.ravel()[ind_lin] - y_ecef[ind_2d_ta]) ** 2
                       + (z_ecef.ravel()[ind_lin] - z_ecef[ind_2d_ta]) ** 2)
        # chord length [m]

        # print(ind_2d_ta)
        # for ind, j in enumerate(ind_lin):
        #     print(j // lon.shape[1], j % lon.shape[1], dist[ind])

        # Mask with ice/glacier (1) and water (9) grid cells
        mask_soil = (soiltyp.ravel()[ind_lin] == 1.0) \
            | (ice.ravel()[ind_lin] > 0.0) \
            | (soiltyp.ravel()[ind_lin] == 9.0)
        frac_diff_abs = np.abs(fr_land.ravel()[ind_lin] - fr_land[ind_2d_ta])
        frac_diff_abs[mask_soil] = np.nan
        mask_frac = (frac_diff_abs == np.nanmin(frac_diff_abs))
        dist[~mask_frac] = np.nan
        ind_sel = np.nanargmin(dist)

        ind_su[0][i] = ind_lin[ind_sel] // lon.shape[1]
        ind_su[1][i] = ind_lin[ind_sel] % lon.shape[1]
        dist_all[i] = dist[ind_sel]
        frac_diff_all[i] = fr_land.ravel()[ind_lin[ind_sel]] \
            - fr_land[ind_2d_ta]

    print(" Statistics (ice -> soil) ".center(79, "-"))
    print("Mean distance: %.0f" % dist_all.mean() + " m")
    print("95th percentile distance: %.0f" % np.percentile(dist_all, 95.0)
          + " m")
    print("Maximal distance: %.0f" % dist_all.max() + " m")
    print("Differences in FR_LAND (min, max, mean): "
          + "%.3f" % frac_diff_all.min() + ", %.3f" % frac_diff_all.max()
          + ", %.5f" % frac_diff_all.mean())
    print("Range of soil types used for replacing ice: "
          + str(int(soiltyp[ind_su[0], ind_su[1]].min()))
          + " - " + str(int(soiltyp[ind_su[0], ind_su[1]].max())))
    print("-" * 79)

# -----------------------------------------------------------------------------
# Derive grid cell replacement indices for soil -> ice
# -----------------------------------------------------------------------------

# To do..

# -----------------------------------------------------------------------------
# Adjust EXTPAR file
# -----------------------------------------------------------------------------

# Variables that are replaced
var_rep = {"2d": ("ICE", "PLCOV_MN", "PLCOV_MX", "LAI_MN", "LAI_MX",
                  "EMIS_RAD", "RSMIN", "URBAN", "FOR_D", "FOR_E",
                  "SKC", "ROOTDP", "NDVI_MAX", "SOILTYP",
                  "FR_LAND", "FR_LAKE", "DEPTH_LK"),
           "3d": ("LU_CLASS_FRACTION",
                  "ALB_DIF12", "ALNID12", "ALUVD12", "NDVI", "NDVI_MRAT")}

# Load data from EXTPAR files
ds = xr.open_dataset(file_mod)
for i in var_rep["2d"]:
    ds[i].values[ind_ta[0], ind_ta[1]] = ds[i].values[ind_su[0], ind_su[1]]
for i in var_rep["3d"]:
    ds[i].values[:, ind_ta[0], ind_ta[1]] \
        = ds[i].values[:, ind_su[0], ind_su[1]]
ds.to_netcdf(file_mod[:-3] + "_adj.nc", format="NETCDF4",
             encoding={"time": {"_FillValue": None},
                       "mlev": {"_FillValue": None}})
# significant difference in created EXTPAR file:
# - new dimension "string1 = 1 ;"
# - char rotated_pole ; -> char rotated_pole(string1) ;
