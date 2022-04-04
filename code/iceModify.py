# This script is used for modifying ICE in the topography-reduced EXTPAR file
# The elevation threshold is set to 2500m

# Load modules
import xarray as xr
import numpy as np
from netCDF4 import Dataset

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
