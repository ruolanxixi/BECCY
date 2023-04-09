# Description: Remap MERIT topography to APHRODITE and CRU grid
#
# Author: Christian R. Steger, April 2023

# Load modules
import os
import numpy as np
import xarray as xr

# Path to folders
path_out = "/project/pr133/csteger/Data/Model/BECCY/remap_topo/"
# output directory

###############################################################################
# Load and process precipitation and temperature data
###############################################################################

# Crop MERIT tiles
file_out = path_out + "MERIT_Eastern_Tibet.nc"
if not os.path.isfile(file_out):
    ds = xr.open_mfdataset("/project/pr133/csteger/Data/DEMs/MERIT/Tiles/"
                           + "MERIT_N??-N??_E0??-E???.nc")
    ds = ds.sel(lon=slice(82.5, 115), lat=slice(42.5, 17.5))
    topo = ds["Elevation"].values
    topo[np.isnan(topo)] = 0.0  # set sea pixels to 0.0 m
    ds["Elevation"].values = topo
    ds.to_netcdf(file_out)
else:
    print("File " + file_out.split("/")[-1] + " already generated")

# Remap MERIT DEM to APHRODITE grid
file_target = "/project/pr133/rxiang/data/obs/tmp/APHRO/day/" \
              + "APHRO.2001-2005.025.mon.nc"
print("cdo griddes " + file_target + " > APHRO_grid.txt")
print("cdo remapcon,APHRO_grid.txt MERIT_Eastern_Tibet.nc"
      + " MERIT_Eastern_Tibet_remap_APHRO.nc")
# Remap MERIT DEM to CRU grid
file_target = "/project/pr133/rxiang/data/obs/tmp/cru/mo/" \
              + "cru.2001-2005.05.mon.nc"
print("cdo griddes " + file_target + " > CRU_grid.txt")
print("cdo remapcon,CRU_grid.txt MERIT_Eastern_Tibet.nc"
      + " MERIT_Eastern_Tibet_remap_CRU.nc")
