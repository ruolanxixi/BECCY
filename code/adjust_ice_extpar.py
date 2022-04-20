# Description: Adjust grid cells covered with permanent snow/ice to modified
#              topography (also adapt linked fields!)
#
# Authors: Ruolan Xiang, Christian R. Steger, IAC ETH Zurich

# Load modules
import xarray as xr
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
# from netCDF4 import Dataset

mpl.style.use("classic")

###############################################################################
# Settings
###############################################################################

# Elevation threshold for permanent snow/ice coverage
elev_thresh = [3800.0, 4640.0, 5200.0]  # percentiles (5.0, 50.0, 95.0) [m]

# Directions to adjust
adj_ice2soil = True
adj_soil2ice = False

# EXTPAR files
file_ref = "/Users/csteger/Desktop/extpar_BECCY_4.4km_merit_unmod_topo.nc"
file_mod = "/Users/csteger/Desktop/extpar_BECCY_4.4km_merit_reduced_topo.nc"

# Window size for searching (nearest) neighbour
wind = 10  # 10 -> (10 + 1 + 10) -> 21 x 21 [-]

###############################################################################
# Adjust EXTPAR file
###############################################################################

# Find variable fields that differ between EXTPAR files
print(" Differing variable fields " .center(79, "-"))
ds = xr.open_dataset(file_ref)
ds_mod = xr.open_dataset(file_mod)
for i in list(ds.variables):
    if not np.all(ds[i].values == ds_mod[i].values):
        print((i + ":").ljust(12) + ds[i].long_name)
ds.close()
ds_mod.close()
print("-" * 79)

# Load file with unmodified topography
ds = xr.open_dataset(file_ref)
topo_unmod = ds["HSURF"].values
fr_land = ds["FR_LAND"].values
ds.close()

ds = xr.open_dataset(file_mod)
topo_mod = ds["HSURF"].values
soiltyp = ds["SOILTYP"].values
ds.close()

# Find grid cells that must be adjusted
mask_mod = (np.abs(topo_unmod - topo_mod) > 0.001)
mask_ice2soil = (mask_mod & (soiltyp == 1.0) & (topo_mod < elev_thresh[0]))
print("Number of grid cells (ice2soil): " + str(mask_ice2soil.sum()))
mask_soil2ice = (mask_mod & (soiltyp != 1.0) & (topo_mod > elev_thresh[2]))
print("Number of grid cells (soil2ice): " + str(mask_soil2ice.sum()))

# Statistics for grid cells that have to be adjusted
print(" Statistics of grid cells that have to be adjusted " .center(79, "-"))
mask = np.zeros(topo_unmod.shape, dtype=bool)
if adj_ice2soil:
    mask[mask_ice2soil] = True
if adj_soil2ice:
    mask[mask_soil2ice] = True
print("Minimal FR_LAND value: %.3f" % fr_land[mask].min())
print("Number of grid cells for which FR_LAND < 0.99: "
      + str((fr_land[mask] < 0.99).sum()))
print("-" * 79)

# -----------------------------------------------------------------------------
# Notes for adjusting values
# -----------------------------------------------------------------------------

# Values not adjusted -> depend on elevation (9)
# - HSURF, FIS
# - S_ORO, Z0, T_CL
# - SSO_STDH, SSO_THETA, SSO_GAMMA, SSO_SIGMA

# Values not adjusted -> fixed (2)
# - lon, lat

# Values not adjusted -> atmosphere (5)
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
# - NDIV_MAX
# - SOILTYP
# - LU_CLASS_FRACTION (size: 23)
# - ALB_DIF12, ALNID12, ALUVD12, NDVI, NDVI_MRAT (size: 12)
# - FR_LAND, FR_LAKE, DEPTH_LK

# Miscellaneous
# - LU_CLASS_FRACTION -> add up to 1.0
# - also replace water/lake properties to keep grid data consistent

# -----------------------------------------------------------------------------
# Adjust EXTPAR fields (ice -> soil)
# -----------------------------------------------------------------------------

# Distance array
x = np.arange(topo_unmod.shape[1], dtype=np.float32)
y = np.arange(topo_unmod.shape[0], dtype=np.float32)
x, y = np.meshgrid(x, y)

if adj_ice2soil:

    indices = zip(*np.where(mask_ice2soil))
    for i, j in indices:

        slic = (slice(i - wind, i + wind + 1), slice(j - wind, j + wind + 1))
        frac_diff_abs = np.abs(fr_land[slic] - fr_land[i, j])  # [-]
        mask = (soiltyp[slic] != 1.0) \
            & (frac_land_diff == frac_land_diff.min())
        dist_sqrt = (x[i, j] - x[slic]) ** 2 + (y[i, j] - y[slic]) ** 2
        dist_sqrt[~mask] = np.nan
        ind_0, ind_1 = np.where(dist_sqrt == np.nanmin(dist_sqrt))
        i_rep, j_rep = ind_0[0] + i, ind_1[0] + j


ds = xr.open_dataset(file_ref)
fr_land = ds["FR_LAND"].values
ice = ds["ICE"].values
ds.close()



ice = ds["ICE"].values
soiltyp = ds["SOILTYP"].values

var = "SKC"
print(ds[var].values[soiltyp == 1].min(), ds[var].values[soiltyp == 1].max())

ds.close()



# -----------------------------------------------------------------------------


# Read EXTPAR file
Path = "/project/pr94/rxiang/data/extpar/"
file = "extpar_12km_878x590_topo1_original.nc"
ds = xr.open_dataset(Path + file)
elev_topo = ds["HSURF"].values
ice_topo = ds["ICE"].values
ds.close()

file = "extpar_12km_878x590.nc"
ds = xr.open_dataset(Path + file)
elev_ctrl = ds["HSURF"].values
lat = ds["lat"].values
lon = ds["lon"].values
ds.close()

elev_diff = elev_ctrl - elev_topo
elev_mask = elev_topo
elev_mask = np.ma.masked_where(elev_diff < 0.1, elev_mask)
elev_mask = np.ma.masked_where(elev_mask > 2500, elev_mask)

mask = np.ma.getmask(elev_mask)
ice_topo = ice_topo * mask

file = "extpar_12km_878x590_topo1.nc"
tar = Dataset(Path + file, 'a')
tar['ICE'][:] = ice_topo[:]
tar.close()
