# Description: Make EXTPAR file consistent with modified land-sea mask and
#              topography (from LGM)
#
# Authors: Christian R. Steger, IAC ETH Zurich

# Load modules
import xarray as xr
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import cartopy.crs as ccrs
from pyproj import CRS, Transformer
from scipy.spatial import cKDTree
from scipy.ndimage import label

mpl.style.use("classic")

# -----------------------------------------------------------------------------
# Settings
# -----------------------------------------------------------------------------

# Input files (LGM: last glacial maximum, PD: present day)
extpar_lgm = "/Users/csteger/Downloads/extpar/12km_LGM/" \
             + "extpar_EAS_ext_12km_merit_LGM.nc"
topo_buffer_lgm = "/Users/csteger/Downloads/extpar/12km_LGM/" \
                  + "topography_buffer_LGM.nc"
lu_buffer_lgm = "/Users/csteger/Downloads/extpar/12km_LGM/" \
                + "extpar_landuse_buffer_LGM.nc"
extpar_pd = "/Users/csteger/Downloads/extpar/12km_PD/" \
            + "extpar_EAS_ext_12km_merit_PD.nc"

# Output file
dir_out = "/Users/csteger/Desktop/output/"
extpar_lgm_out = dir_out + "extpar_EAS_ext_12km_merit_LGM_consistent.nc"

# Miscellaneous
keep_caspian_sea = True  # keep Caspian Sea in LGM
avoid_new_lakes = True
# treat new water bodies disconnected from the ocean as land depressions
num_nn = 1  # number of nearest neighbors grid cells to consider (1, 25)

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
# Load and check data
# -----------------------------------------------------------------------------

# Read data from topography buffer
data_buffer = {}
ds = xr.open_dataset(topo_buffer_lgm)
for i in ("HSURF", "SSO_STDH", "SSO_THETA", "SSO_GAMMA", "SSO_SIGMA",
          "FR_LAND_TOPO", "Z0_TOPO", "S_ORO"):
    if i in tuple(ds.variables):
        data_buffer[i] = ds[i].values.squeeze()
ds.close()

# Read data from landuse buffer
ds = xr.open_dataset(lu_buffer_lgm)
for i in ("Z0", ):
    if i in tuple(ds.variables):
        data_buffer[i] = ds[i].values.squeeze()
ds.close()

# Check for differences in topography-derived fields
ds = xr.open_dataset(extpar_lgm)
rlon = ds["rlon"].values
rlat = ds["rlat"].values
ccrs_rot = ccrs.RotatedPole(
    pole_latitude=ds["rotated_pole"].grid_north_pole_latitude,
    pole_longitude=ds["rotated_pole"].grid_north_pole_longitude)
for i in ("HSURF", "SSO_STDH", "SSO_THETA", "SSO_GAMMA", "SSO_SIGMA", "S_ORO"):
    if i in tuple(ds.variables):
        flag = np.all(ds[i].values == data_buffer[i])
        print(i + " identical: " + str(flag))
        if not flag:
            data_plot = (ds[i].values - data_buffer[i])
            data_plot = np.ma.masked_where(data_plot == 0.0, data_plot)
            plt.figure(figsize=(11, 5))
            ax = plt.axes(projection=ccrs_rot)
            plt.pcolormesh(rlon, rlat, data_plot)
            ax.set_aspect("auto")
            gl = ax.gridlines(crs=ccrs.PlateCarree(), linewidth=0.5,
                              color="black", alpha=0.8, linestyle="-",
                              draw_labels=True, dms=True,
                              x_inline=False, y_inline=False)
            gl.top_labels = False
            gl.right_labels = False
            ax.coastlines(resolution="50m", linewidth=0.5)
            plt.title(i + ": EXTPAR minus buffer")
            plt.colorbar()
ds.close()

# Load EXTPAR data for Last Glacial Maximum
ds = xr.open_dataset(extpar_lgm)
fr_land = ds["FR_LAND"].values
soiltyp = ds["SOILTYP"].values  # 1 - 8: land, 9: water
hsurf = ds["HSURF"].values
lon = ds["lon"].values
lat = ds["lat"].values
z0_tot = ds["Z0"].values
ds.close()

# Load lake fraction from present-day extpar file
ds = xr.open_dataset(extpar_pd)
fr_lake = ds["FR_LAKE"].values
depth_lk = ds["DEPTH_LK"].values
ds.close()
# -> 'FR_LAKE' and 'DEPTH_LK' from LGM EXTPAR file are erroneous!

# Check consistency between land, sea, lakes and soil types in EXTPAR file
print(fr_land[soiltyp == 9].max())  # water
print(fr_land[soiltyp != 9].min())  # land
print(np.unique(soiltyp[fr_lake > 0.5]))

# Recompute total surface roughness
z0_tot_rec = data_buffer["Z0"] + data_buffer["Z0_TOPO"]
mask_water = (fr_land < 0.5)
z0_tot_rec[mask_water] = np.maximum(1.0e-6, z0_tot_rec[mask_water])
z0_tot_rec[~mask_water] = np.maximum(1.0e-2, z0_tot_rec[~mask_water])
print("Z0_tot identical: " + str(np.all(z0_tot == z0_tot_rec)))

# -----------------------------------------------------------------------------
# Derive mask for adapting grid cells with modified land-sea-fractions
# -----------------------------------------------------------------------------

# Compute increase in grid cell's land fraction from PD -> LGM
fr_land_inc = ((data_buffer["FR_LAND_TOPO"] - fr_lake) - fr_land)
# -> assume lake fraction to remain constant...
num_gc = (fr_land_inc < 0.0).sum()
print("Number of grid cells with decreasing land fraction: " + str(num_gc))
print("Smallest, 5th- and 9th-smallest value: "
      + ", ".join(["%.8f" % i for i in
                   np.sort(fr_land_inc.ravel())[[0, 4, 8]]]))
# -> ignore grid cells with decreasing fractions...
if keep_caspian_sea:
    mask_caspian_sea = (hsurf != data_buffer["HSURF"])
    fr_land_inc[mask_caspian_sea] = 0.0
fr_land_lgm = (fr_land + fr_land_inc)
fr_land_lgm_bin = (fr_land_lgm >= 0.5).astype(np.int32)
fr_land_bin = (fr_land >= 0.5).astype(np.int32)
mask_fill = np.zeros(fr_land.shape, dtype=bool)
if avoid_new_lakes:
    binary = (~(fr_land_lgm_bin - fr_land_bin).astype(bool)).astype(int)
    thresh_gc_fill = 100
    structure = [[1, 1, 1], [1, 1, 1], [1, 1, 1]]
    label_arr, num_features = label(binary, structure=structure)
    for i in range(1, num_features + 1):
        mask = (label_arr == i)
        if (mask.sum() <= thresh_gc_fill) and np.all(soiltyp[mask] == 9):
            mask_fill[mask] = True
    print("Number of grid cells that represents holes: "
          + str(mask_fill.sum()))
    fr_land_inc[mask_fill] = (1.0 - fr_land[mask_fill])

# Plot land fractions (PD and LGM)
# -----------------------------------------------------------------------------
fig = plt.figure(figsize=(16, 10))
gs = gridspec.GridSpec(2, 3, left=0.1, bottom=0.1, right=0.9,
                       top=0.9, hspace=0.11, wspace=0.1,
                       width_ratios=[1.0, 1.0, 0.05])
# -----------------------------------------------------------------------------
levels = np.arange(0.0, 1.0, 0.05)
cmap = plt.get_cmap("Greens")
norm = mpl.colors.BoundaryNorm(levels, ncolors=cmap.N)
ax = plt.subplot(gs[0, 0], projection=ccrs_rot)
plt.pcolormesh(rlon, rlat, fr_land, cmap=cmap, norm=norm)
ax.set_aspect("auto")
ax.coastlines(resolution="50m", linewidth=0.5)
plt.title("Fractional land coverage (present-day)", fontsize=12,
          fontweight="bold", y=1.01)
# -----------------------------------------------------------------------------
ax = plt.subplot(gs[1, 0], projection=ccrs_rot)
plt.pcolormesh(rlon, rlat, fr_land_lgm, cmap=cmap, norm=norm)
ax.set_aspect("auto")
ax.coastlines(resolution="50m", linewidth=0.5)
plt.title("Fractional land coverage (Last Glacial Maximum)", fontsize=12,
          fontweight="bold", y=1.01)
# -----------------------------------------------------------------------------
levels = np.arange(0.0, 1.0, 0.05)
cmap = plt.get_cmap("Blues")
norm = mpl.colors.BoundaryNorm(levels, ncolors=cmap.N)
ax = plt.subplot(gs[0, 1], projection=ccrs_rot)
plt.pcolormesh(rlon, rlat, fr_land_inc, cmap=cmap, norm=norm)
ax.set_aspect("auto")
ax.coastlines(resolution="50m", linewidth=0.5)
plt.title("Difference in fractional land coverage", fontsize=12,
          fontweight="bold", y=1.01)
# -----------------------------------------------------------------------------
cmap = mpl.colors.ListedColormap(["white", "forestgreen", "orangered"])
bounds = [-0.5, 0.5, 1.5, 2.5]
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
ax = plt.subplot(gs[1, 1], projection=ccrs_rot)
data_plot = fr_land_lgm_bin.copy()
data_plot[mask_fill] = 2
plt.pcolormesh(rlon, rlat, data_plot, cmap=cmap, norm=norm)
ax.set_aspect("auto")
ax.coastlines(resolution="50m", linewidth=0.5)
plt.title("Binary land coverage (Last Glacial Maximum)", fontsize=12,
          fontweight="bold", y=1.01)
# -----------------------------------------------------------------------------
fig.savefig(dir_out + "Land_fractions.png", dpi=300, bbox_inches="tight")
plt.close(fig)

# -----------------------------------------------------------------------------
# Modify/adapt grid cells
# -----------------------------------------------------------------------------

# Variables already modified (-> derived from modified MERIT DEM)
# - S_ORO, HSURF, FIS (3)
# - SSO_STDH, SSO_THETA, SSO_GAMMA, SSO_SIGMA (4)

# Variables kept constant
# - lon, lat (2)
# - AER_BC12, AER_DUST12, AER_ORG12, AER_SO412, AER_SS12 (5)

# Variables filled from nearest neighbour land cell
var_fill_nn = (
    # ------------------------- 2D -------------------------
    "ICE", "PLCOV_MN", "PLCOV_MX", "LAI_MN", "LAI_MX",
    "EMIS_RAD", "RSMIN", "FOR_D", "FOR_E", "SKC",
    "ROOTDP", "NDVI_MAX", "T_CL", "SOILTYP",
    # ------------------------- 3D -------------------------
    "LU_CLASS_FRACTION",
    "ALB_DIF12", "ALNID12", "ALUVD12",
    "NDVI", "NDVI_MRAT"
    # ------------------------------------------------------
    )  # (20)

# Remaining variables
# - FR_LAND: from data_buffer["FR_LAND_TOPO"] (1)
# - Z0: recompute from data_buffer["Z0"], data_buffer["Z0_TOPO"] (1)
# - URBAN: set to 0.0 everywhere (1)
# - FR_LAKE, DEPTH_K:  replace from unmodified EXTPAR file (2)

# Construct k-d tree
crs_ecef = CRS.from_dict({"proj": "geocent", "ellps": "sphere"})
crs_latlon = CRS.from_dict({"proj": "latlong", "ellps": "sphere"})
trans = Transformer.from_crs(crs_latlon, crs_ecef, always_xy=True)
x_ecef, y_ecef, z_ecef = trans.transform(lon, lat, np.zeros_like(lon))
indices_rep = np.where(soiltyp != 9)  # indices of 'replacement' grid cells
pts_land = np.vstack((x_ecef[indices_rep],
                      y_ecef[indices_rep],
                      z_ecef[indices_rep])).transpose()
tree = cKDTree(pts_land)

# Compute nearest neighbours
mask_gc_change = (fr_land_inc > 1.e-5)  # 0.01 per mille
indices_gc_change = np.where(mask_gc_change)
pts_new = np.vstack([x_ecef[indices_gc_change],
                     y_ecef[indices_gc_change],
                     z_ecef[indices_gc_change]]).transpose()
dist, ind_nn = tree.query(pts_new, k=num_nn)
if num_nn == 1:
    dist = dist[:, np.newaxis]
    ind_nn = ind_nn[:, np.newaxis]

# Plot nearest-neighbour distance
levels = np.arange(0.0, 450.0, 50.0)
cmap = plt.get_cmap("Spectral_r")
norm = mpl.colors.BoundaryNorm(levels, ncolors=cmap.N, extend="max")
fig = plt.figure(figsize=(11, 6))
ax = plt.axes(projection=ccrs_rot)
dist_arr_2d = np.empty(mask_gc_change.shape, dtype=np.float32)
dist_arr_2d.fill(np.nan)
dist_arr_2d[mask_gc_change] = dist.mean(axis=1) / 1000.0  # [km]
plt.pcolormesh(rlon, rlat, dist_arr_2d, cmap=cmap, norm=norm)
# ax.set_aspect("auto")
ax.coastlines(resolution="50m", linewidth=0.5)
plt.title("Distance to nearest land grid cell(s) (along chord line) [km]",
          fontsize=12, fontweight="bold", y=1.01)
plt.colorbar()
fig.savefig(dir_out + "Grid_cell_nearest_neighbour_distance.png",
            dpi=300, bbox_inches="tight")
plt.close(fig)

# Fill or recompute variables
ds = xr.open_dataset(extpar_lgm)
var_not_found = set(var_fill_nn) - set(ds.data_vars)
if len(var_not_found) != 0:
    print("Warning: not all variables are present in the EXTPAR file")
    var_fill_nn = set(var_fill_nn).intersection(set(ds.keys()))
for i in range(mask_gc_change.sum()):

    # Indices of grid cell that is modified
    ind_0 = indices_gc_change[0][i]
    ind_1 = indices_gc_change[1][i]

    # Indices of nearest neighbour land grid cell(s)
    ind_nn_0 = indices_rep[0][ind_nn[i, :]]
    ind_nn_1 = indices_rep[1][ind_nn[i, :]]

    # Fill-in land variables in case grid cell turns from water to land
    fr_land_new = (ds["FR_LAND"].values[ind_0, ind_1]
                   + fr_land_inc[ind_0, ind_1]).clip(max=1.0)
    if (soiltyp[ind_0, ind_1] == 9) and (fr_land_new >= 0.5):
        # ---------------------------------------------------------------------
        for j in var_fill_nn:
            if j != "SOILTYP":
                ds[j].values[..., ind_0, ind_1] \
                    = ds[j].values[..., ind_nn_0, ind_nn_1].mean(axis=-1)
            else:
                unique, unique_counts \
                    = np.unique(ds[j].values[ind_nn_0, ind_nn_1],
                                return_counts=True)
                ds[j].values[..., ind_0, ind_1] \
                    = unique[np.argmax(unique_counts)]
        # ---------------------------------------------------------------------
        z0_tot_rec = data_buffer["Z0"][ind_nn_0, ind_nn_1].mean(axis=-1) \
            + data_buffer["Z0_TOPO"][ind_0, ind_1]
        ds["Z0"].values[ind_0, ind_1] = np.maximum(1.0e-2, z0_tot_rec)

    # Update land fraction
    ds["FR_LAND"].values[ind_0, ind_1] = fr_land_new

    if i % 5000 == 0:
        print("Grid cell processed: " + str(i) + "/"
              + str(mask_gc_change.sum()))

# Remaining variables
ds["URBAN"].values[:] = 0.0
ds["FR_LAKE"].values[:] = fr_lake
ds["DEPTH_LK"].values[:] = depth_lk

ds.to_netcdf(extpar_lgm_out,
             encoding={"time": {"_FillValue": None},
                       "mlev": {"_FillValue": None}})

# -----------------------------------------------------------------------------
# Brief 'sanity check' of final EXTPAR file
# -----------------------------------------------------------------------------

ds = xr.open_dataset(extpar_lgm_out)
print(ds["FR_LAND"].values[ds["SOILTYP"].values == 9].max())  # water
print(ds["FR_LAND"].values[ds["SOILTYP"].values != 9].min())  # land
print(np.unique(ds["SOILTYP"].values[ds["FR_LAKE"].values > 0.5]))
print(np.all(np.isfinite(ds["T_CL"].values[ds["SOILTYP"].values != 9])))
mask_land = (ds["SOILTYP"].values > 1) & (ds["SOILTYP"].values < 9)
print(np.all(ds["PLCOV_MX"].values[mask_land] > 0.0))
print(np.all(ds["LAI_MN"].values[mask_land] > 0.0))
print(np.all(ds["LAI_MX"].values[mask_land] > 0.0))
print(np.unique(ds["SOILTYP"].values))
ds.close()
