# Description: Compute precipitation quantities for the Hengduan Mountains for
#              the control experiment as well as the reduced and envelope
#              topography experiment (and also difference).
#
# Author: Christian R. Steger, April 2023

# Load modules
import os
import numpy as np
import xarray as xr

# Path to folders
path_reg_masks = "/project/pr133/csteger/Data/Model/BECCY/region_masks/"

###############################################################################
# Load and process precipitation and temperature data
###############################################################################

# Load region masks
regions = ["ET", "HM", "HMU", "HMC", "HMUS", "HMUN"]
ds = xr.open_dataset(path_reg_masks + "CTRL04_region_masks.nc")
region_masks = {i: ds[i].values.astype(bool) for i in regions}
ds.close()

# Load precipiation
exp = ["ctrl", "topo1", "topo2"]
prec = {}
for i in exp:
    ds = xr.open_dataset("/project/pr133/rxiang/data/cosmo/EAS04_" + i
                         + "/mon/TOT_PREC/2001-2005.TOT_PREC.nc")
    ds = ds.sel(time=(ds["time.month"] >= 5) & (ds["time.month"] <= 9))
    prec[i] = ds["TOT_PREC"].mean(axis=0).values  # mm/day
    ds.close()

# Compute spatially integrated quantity over regions
for i in ["HM", "HMC", "HMU"]:
    print("--------------" + i + "--------------")
    print("ctrl: %.1f" % prec["ctrl"][region_masks[i]].mean() + " mm/day")
    for j in exp[1:]:
        diff_abs = (prec[j][region_masks[i]].mean()
                    - prec["ctrl"][region_masks[i]].mean())
        diff_rel = (prec[j][region_masks[i]].mean()
                    / prec["ctrl"][region_masks[i]].mean() - 1.0) * 100.0
        print(j + ": %.1f" % prec[j][region_masks[i]].mean()
              + " ({0:+.1f}".format(diff_abs) + ")" + " mm/day"
              + " ({0:+.1f}".format(diff_rel) + " %)")

# continue with other precipiation quantities...
