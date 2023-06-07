# Description: Make EXTPAR file consistent with modified land-sea mask and
#              topography (from LGM)
#
# Authors: Christian R. Steger, IAC ETH Zurich

# Load modules
import xarray as xr
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from pyproj import CRS, Transformer
from scipy.spatial import cKDTree

mpl.style.use("classic")

# -----------------------------------------------------------------------------
# Settings
# -----------------------------------------------------------------------------

# In-/output folder
# path = "/Users/csteger/Downloads/extpar/4km/"
path = "/Users/csteger/Downloads/extpar/12km/"

# File names
# extpar_in = "extpar_BECCY_4.4km_merit.nc"
# extpar_out = "extpar_BECCY_4.4km_merit_LGM_final.nc"
extpar_in = "extpar_EAS_ext_12km_merit.nc"
extpar_out = "extpar_EAS_ext_12km_merit_LGM_final.nc"

topo_buffer = "topography_buffer.nc"
lu_buffer = "extpar_landuse_buffer.nc"

# -----------------------------------------------------------------------------
# Notes on fields in EXTPAR actually used in COSMO
# -----------------------------------------------------------------------------

# Both in EXTPAR and laf-file:
# - FR_LAND
# - PLCOV_MN, PLCOV_MX (EXTPAR) -> PLCOV(time, rlat, rlon) (laf-file)
# - LAI_MN, LAI_MX (EXTPAR) -> LAI(time, rlat, rlon) (laf-file)
# - Z0: Roughness length
# - S_ORO -> S_ORO_MAX
# - FOR_D, FOR_E
# - ROOTDP: Root depth
# - HSURF
# - SSO_STDH: standard deviation of subgrid scale height
# - SOILTYP
# - AER_BC, AER_DUST, AER_ORG, AER_SO4, AER_SS

# In EXTPAR but not in laf-file:
# - ICE: Ice fraction due to GLOBCOVER20009 Data
# - EMIS_RAD: longwave surface emissivity
# - RSMIN: Minimal stomata resistance
# - URBAN: urban area fraction
# - SKC: Skin conductivity
# - NDVI_MAX: NDVI yearly maximum for climatology 1998-2003
# - FIS: Geopotential (S)
# - T_CL: CRU near surface temperature climatology
# - FR_LAKE: fraction lake
# - DEPTH_LK: Lake depth
# - LU_CLASS_FRAC: Fraction of GLOBCOVER2009 land use classes in target grid
#   element

# Unclear:
# - ALB_DIF, ALNID, ALUVD
#   (old albedo: ALB_DRY, ALB_SAT)
# - NDVI, NDVI_MRAT (?)

# -----------------------------------------------------------------------------
# Other notes
# -----------------------------------------------------------------------------

# EXTPAR soil types
# 1: ice and glacier
# 2: rock, lithosols
# 3: sand
# 4: sandy loam
# 5: loam (default soil type)
# 6: loamy clay
# 7: clay
# 8: histosols (peat)
# 9: water

# -----------------------------------------------------------------------------
# Process EXTPAR file
# -----------------------------------------------------------------------------

# Read data from topography buffer
data_buffer = {}
ds = xr.open_dataset(path + topo_buffer)
for i in ("HSURF", "SSO_STDH", "SSO_THETA", "SSO_GAMMA", "SSO_SIGMA",
          "FR_LAND_TOPO", "Z0_TOPO", "S_ORO"):
    if i in tuple(ds.variables):
        data_buffer[i] = ds[i].values.squeeze()
ds.close()

# Read data from landuse buffer
ds = xr.open_dataset(path + lu_buffer)
for i in ("Z0", ):
    if i in tuple(ds.variables):
        data_buffer[i] = ds[i].values.squeeze()
ds.close()

# Check equality of topography-derived fields
ds = xr.open_dataset(path + extpar_in)
for i in ("HSURF", "SSO_STDH", "SSO_THETA", "SSO_GAMMA", "SSO_SIGMA", "S_ORO"):
    flag = np.all(ds[i].values == data_buffer[i])
    print(i + " identical: " + str(flag))
    if not flag:
        plt.pcolormesh(ds[i].values - data_buffer[i])
        plt.title(i)
ds.close()

# Load EXTPAR data
ds = xr.open_dataset(path + extpar_in)
fr_land = ds["FR_LAND"].values
fr_lake = ds["FR_LAKE"].values
soiltyp = ds["SOILTYP"].values
rlon = ds["rlon"].values
rlat = ds["rlat"].values
lon = ds["lon"].values
lat = ds["lat"].values
ccrs_rot = ccrs.RotatedPole(
    pole_latitude=ds["rotated_pole"].grid_north_pole_latitude,
    pole_longitude=ds["rotated_pole"].grid_north_pole_longitude)
z0_tot = ds["Z0"].values
ds.close()

# Check consistency between land, sea, lake and soil types
print(fr_land[soiltyp == 9].max())
print(fr_land[soiltyp != 9].min())
print(np.unique(soiltyp[fr_lake > 0.5]))

# Recompute total surface roughness
z0_tot_rec = data_buffer["Z0"] + data_buffer["Z0_TOPO"]
mask = (fr_land < 0.5)
z0_tot_rec[mask] = np.maximum(1.E-6, z0_tot_rec[mask])
z0_tot_rec[~mask] = np.maximum(1.E-2, z0_tot_rec[~mask])
print("Z0_tot identical: " + str(np.all(z0_tot == z0_tot_rec)))

# Mask with sea grid cells that are land in LGM
mask_sea2land = (data_buffer["FR_LAND_TOPO"] >= 0.5) & (fr_land < 0.5) \
                & (fr_lake < 0.5)
if not np.all(soiltyp[mask_sea2land] == 9):
    raise ValueError("Inconsistency between land fraction and soil type")
print("Number of grid cells that are changed from sea to land: "
      + str(mask_sea2land.sum()))

# Remove Caspian Sea from mask
mask_caspian_sea = (ds["HSURF"].values != data_buffer["HSURF"])
mask_sea2land[mask_caspian_sea] = False

# Test plot
data_plot = np.zeros(mask_sea2land.shape, dtype=np.int32)
data_plot[soiltyp != 9] = 1
data_plot[mask_sea2land] = 2
cmap = mpl.colors.ListedColormap(["white", "grey", "orangered"])
bounds = [-0.5, 0.5, 1.5, 2.5]
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
plt.figure(figsize=(11, 6))
ax = plt.axes(projection=ccrs_rot)
plt.pcolormesh(rlon, rlat, data_plot, cmap=cmap, norm=norm)
# ax.set_aspect("auto")
gl = ax.gridlines(crs=ccrs.PlateCarree(), linewidth=0.5, color="black",
                  alpha=0.8, linestyle="-", draw_labels=True, dms=True,
                  x_inline=False, y_inline=False)
gl.top_labels = False
gl.right_labels = False
ax.coastlines(resolution="50m", linewidth=0.5)

# Construct k-d tree
crs_ecef = CRS.from_dict({"proj": "geocent", "ellps": "sphere"})
crs_latlon = CRS.from_dict({"proj": "latlong", "ellps": "sphere"})
trans = Transformer.from_crs(crs_latlon, crs_ecef, always_xy=True)
x_ecef, y_ecef, z_ecef = trans.transform(lon, lat, np.zeros_like(lon))
mask_rep = (soiltyp != 9)
pts_gc = np.vstack((x_ecef[mask_rep], y_ecef[mask_rep],
                    z_ecef[mask_rep])).transpose()
tree = cKDTree(pts_gc)

# Fill land (surface) characteristics from nearest neighbour land cell
var_fill = ("ICE", "PLCOV_MN", "PLCOV_MX", "LAI_MN", "LAI_MX",
            "EMIS_RAD", "RSMIN", "URBAN", "FOR_D", "FOR_E", "SKC",
            "ROOTDP", "NDVI_MAX", "T_C", "SOILTYP", "LU_CASS_FRAC",
            "ALB_DIF12", "ALNID12", "ALUVD12",
            "NDVI", "NDVI_MRAT")  # (21)
# variables not considered:
# - FR_LAND, Z0, S_ORO, HSURF, FIS (5)
# - SSO_STDH, SSO_THETA, SSO_GAMMA, SSO_SIGMA (4)
# - lon, lat (2)
# - FR_LAKE (default: 0.0), DEPTH_K (default: -1.0) (2)
# - AER_BC12, AER_DUST12, AER_ORG12, AER_SO412, AER_SS12 (5)
ind_rep_0, ind_rep_1 = np.where(mask_rep)
distances = np.empty(soiltyp.shape, dtype=np.float32)
distances.fill(np.nan)
ds = xr.open_dataset(path + extpar_in)
count = 0
for ind_0, ind_1 in zip(*np.where(mask_sea2land)):

    # Find indices of source grid cell
    point = np.array([x_ecef[ind_0, ind_1],
                      y_ecef[ind_0, ind_1],
                      z_ecef[ind_0, ind_1]])
    dist, ind = tree.query(point, k=1)
    distances[ind_0, ind_1] = dist / 1000.0  # [km]
    ind_sou_0 = ind_rep_0[ind]
    ind_sou_1 = ind_rep_1[ind]

    # Loop through variables
    for i in var_fill:
        if i in tuple(ds.variables):
            ds[i][..., ind_0, ind_1].values \
                = ds[i][..., ind_sou_0, ind_sou_1].values

    print(count)
    count += 1

ds.to_netcdf(path + extpar_out) ############################################### not yet check if correct -> continue...

# Plot nearest-neighbour distance
levels = np.arange(0.0, 450.0, 50.0)
cmap = plt.get_cmap("Spectral_r")
norm = mpl.colors.BoundaryNorm(levels, ncolors=cmap.N, extend="max")
plt.figure(figsize=(11, 6))
ax = plt.axes(projection=ccrs_rot)
plt.pcolormesh(rlon, rlat, distances, cmap=cmap, norm=norm)
# ax.set_aspect("auto")
gl = ax.gridlines(crs=ccrs.PlateCarree(), linewidth=0.5, color="black",
                  alpha=0.8, linestyle="-", draw_labels=True, dms=True,
                  x_inline=False, y_inline=False)
gl.top_labels = False
gl.right_labels = False
ax.coastlines(resolution="50m", linewidth=0.5)
plt.title("Nearest neighbour distance along chord line [km]",
          fontsize=12, fontweight="bold", y=1.01)
plt.colorbar()









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

# Load data from EXTPAR files
ds = xr.open_dataset(file_ref)
topo_unmod = ds["HSURF"].values
fr_land = ds["FR_LAND"].values
soiltyp = ds["SOILTYP"].values
alb_dif12 = ds["ALB_DIF12"].values
alnid12 = ds["ALNID12"].values
aluvd12 = ds["ALUVD12"].values
lon = ds["lon"].values  # [degree]
lat = ds["lat"].values  # [degree]
ds.close()
ds = xr.open_dataset(file_mod)
topo_mod = ds["HSURF"].values
ds.close()

# Check available soil types in data
print(" Occurring soil types in EXTPAR file " .center(79, "-"))
soil_types_tc = {1: "ice and glacier", 2: "rock, lithosols", 3: "sand",
                 4: "sandy loam", 5: "loam (default soil type)",
                 6: "loamy clay", 7: "clay", 8: "histosols (peat)",
                 9: "water", 10: "alkali flat", 11: "shifting sand, dunes",
                 12: "urban, human disturbed", 255: "unknown"}
for i in np.unique(soiltyp.astype(int)):
    print(str(i).ljust(2) + ": " + soil_types_tc[i])

# Check consistency of soiltyp == 1 (ice and glacier) and albedo values
if (np.any(np.diff((alb_dif12 == 0.7), axis=0))
        or np.any(np.diff((alnid12 == 0.7), axis=0))
        or np.any(np.diff((aluvd12 == 0.7), axis=0))):
    raise ValueError("Albedo values of glaciated area temporally inconsistent")
arr_bool = np.concatenate((soiltyp[np.newaxis, :, :] == 1.0,
                           alb_dif12[:1, :, :] == 0.7,
                           alnid12[:1, :, :] == 0.7,
                           aluvd12[:1, :, :] == 0.7), axis=0)
if np.any(np.diff(arr_bool, axis=0)):
    raise ValueError("Fields of 'soiltyp == 1' and albedo are spatially "
                     + "inconsistent")

# Find adjustable grid cells
mask_mod = (np.abs(topo_unmod - topo_mod) > 0.001)
mask_ice2soil = (mask_mod & (soiltyp == 1.0) & (topo_mod < elev_thresh[0]))
mask_soil2ice = (mask_mod & (soiltyp != 1.0) & (topo_mod > elev_thresh[2]))

# Print statistics of potential correction (ice -> soil and soil -> ice)
masks = {"ice -> soil": mask_ice2soil, "soil -> ice": mask_soil2ice}
for i in list(masks.keys()):
    print((" " + i + ": potential correction ").center(79, "-"))
    txt = str(masks[i].sum()) + ", " \
        + "%.3f" % (masks[i].sum() / mask_mod.sum() * 100.0) + " %, " \
        + "%.3f" % (masks[i].sum() / mask_mod.size * 100.0) + " %"
    print("Grid cell (total, percentage of modified and all): " + txt)
    txt = "[%.0f" % topo_unmod[masks[i]].min() + " m, " \
          + "%.0f" % topo_unmod[masks[i]].max() + " m]" + " -> " \
          + "[%.0f" % topo_mod[masks[i]].min() + " m, " \
          + "%.0f" % topo_mod[masks[i]].max() + " m]"
    print("Elevation range (unmod -> mod): " + txt)

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

# Values adjusted -> depend on surface/soil (21)
# - PLCOV_MN, PLCOV_MX, LAI_MN, LAI_MX
# - EMIS_RAD
# - RSMIN
# - URBAN
# - FOR_D, FOR_E
# - SKC
# - ROOTDP
# - NDVI_MAX
# - SOILTYP
# - ALB_DIF12, ALNID12, ALUVD12, NDVI, NDVI_MRAT (size: 12)
# - FR_LAND, FR_LAKE, DEPTH_LK

# Dropped values -> not used by INT2LM (and inconsistent with modification) (2)
# - ICE
# - LU_CLASS_FRACTION (size: 23)

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
        mask_ice_water = (soiltyp.ravel()[ind_lin] == 1.0) \
            | (soiltyp.ravel()[ind_lin] == 9.0)
        frac_diff_abs = np.abs(fr_land.ravel()[ind_lin] - fr_land[ind_2d_ta])
        frac_diff_abs[mask_ice_water] = np.nan
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

if adj_soil2ice:

    raise ValueError("Adjusting grid cells from soil -> ice not implemented")

# -----------------------------------------------------------------------------
# Adjust EXTPAR file
# -----------------------------------------------------------------------------

# Variables that are replaced
var_rep = {"2d": ("PLCOV_MN", "PLCOV_MX", "LAI_MN", "LAI_MX",
                  "EMIS_RAD", "RSMIN", "URBAN", "FOR_D", "FOR_E",
                  "SKC", "ROOTDP", "NDVI_MAX", "SOILTYP",
                  "FR_LAND", "FR_LAKE", "DEPTH_LK"),
           "3d": ("ALB_DIF12", "ALNID12", "ALUVD12", "NDVI", "NDVI_MRAT")}

# Load data from EXTPAR files
ds = xr.open_dataset(file_mod)
ds = ds.drop(["ICE", "LU_CLASS_FRACTION"])
gc_lake_before = (ds["FR_LAKE"].values != 0).sum()
for i in var_rep["2d"]:
    ds[i].values[ind_ta[0], ind_ta[1]] = ds[i].values[ind_su[0], ind_su[1]]
for i in var_rep["3d"]:
    ds[i].values[:, ind_ta[0], ind_ta[1]] \
        = ds[i].values[:, ind_su[0], ind_su[1]]
gc_lake_after = (ds["FR_LAKE"].values != 0).sum()
ds.to_netcdf(file_mod[:-3] + "_adj.nc", format="NETCDF4",
             encoding={"time": {"_FillValue": None},
                       "mlev": {"_FillValue": None}})

print("Number of grid cells with non-zero lake fraction: ")
print("Before modification: " + str(gc_lake_before))
print("After modification: " + str(gc_lake_after))
print("-" * 79)

# Notes:
# - Difference in modified EXTPAR file:
#   - new dimension "string1 = 1 ;"
#   - char rotated_pole ; -> char rotated_pole(string1) ;
#   -> if issue -> solvable by encoding; e.g. "rotated_pole": {"dtype": ""}
