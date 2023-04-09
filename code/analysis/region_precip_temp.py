# Description: Compute Hengduan Mountains regions masks for different products
#
# Author: Christian R. Steger, April 2023

# Load modules
import os
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.gridspec as gridspec

mpl.style.use("classic")

# Change latex fonts
mpl.rcParams["mathtext.fontset"] = "custom"
# custom mathtext font (set default to Bitstream Vera Sans)
mpl.rcParams["mathtext.default"] = "rm"
mpl.rcParams["mathtext.rm"] = "Bitstream Vera Sans"

# Path to folders
path_reg_masks = "/project/pr133/csteger/Data/Model/BECCY/region_masks/"
path_out = "/project/pr133/csteger/Data/Model/BECCY/region_precip_temp/"
# output directory

###############################################################################
# Load and process precipitation and temperature data
###############################################################################

# -----------------------------------------------------------------------------
# Precipitation
# -----------------------------------------------------------------------------

# List of products
prec_prod = {
    # -------------------------------------------------------------------------
    "APHRODITE": {"file": "/project/pr133/rxiang/data/obs/pr/APHRO/"
                          + "APHRO_2001-2005_ymonmean.nc",
                  "var_name": "precip",
                  "units_in": "mm/day"},
    "GPCC": {"file": "/project/pr133/rxiang/data/obs/pr/GPCC/mo/"
                     + "GPCC.2001-2005.025.mon.nc",
             "var_name": "precip",
             "units_in": "mm/month"},
    "IMERG": {"file": "/project/pr133/rxiang/data/obs/pr/IMERG/day_old/"
                      + "IMERG.ydaymean.2001-2005.mon.nc4",
              "var_name": "pr",
              "units_in": "mm/day"},
    "ERA5": {"file": "/project/pr133/rxiang/data/era5/pr/mo/"
                     + "era5.mo.2001-2005.mon.nc",
             "var_name": "tp",
             "units_in": "m/day"},
    "PBCOR": {"file": "/project/pr133/csteger/Data/Observations/"
                      + "PBCOR_V1.0/CHELSA_V12.nc",
              # CHELSA_V12.nc  CHPclim_V1.nc  WorldClim_V2.nc
              "var_name": "corr_P_monthly",
              "units_in": "mm/month"},
    # -------------------------------------------------------------------------
    "CTRL04": {"file": "/project/pr133/rxiang/data/cosmo/EAS04_ctrl/mon/"
                       + "TOT_PREC/2001-2005.TOT_PREC.nc",
               "var_name": "TOT_PREC",
               "units_in": "mm/day"},
    "CTRL11": {"file": "/project/pr133/rxiang/data/cosmo/EAS11_ctrl/mon/"
                       + "TOT_PREC/2001-2005.TOT_PREC.nc",
               "var_name": "TOT_PREC",
               "units_in": "mm/day"},
    # -------------------------------------------------------------------------
    }

# Compute monthly precipitation for products and regions
regions = ["ET", "HM", "HMU", "HMC", "HMUS", "HMUN"]
months_num = np.arange(1, 13)
mask_rainy = (months_num >= 5) & (months_num <= 9)  # rainy season (MJJAS)
mask_dry = (months_num <= 3) | (months_num >= 11)  # dry season (NDJFM)
months_char = ["J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D"]
prec_prod_reg = {}
prec_prod_reg_unc = {}  # uncertainty
for i in prec_prod.keys():

    print("Process product " + i)
    prec_prod_reg[i] = {}
    prec_prod_reg_unc[i] = {}

    # Load precipitation data
    ds = xr.open_dataset(prec_prod[i]["file"])
    prec = ds[prec_prod[i]["var_name"]].values
    if i == "PBCOR":
        corr_fac_an_lb = ds["corr_fac_annual_lower_bound"].values
        corr_fac_an_ub = ds["corr_fac_annual_upper_bound"].values
        precip_orig_an = ds["orig_P_annual"].values  # [mm/year]
        precip_low = precip_orig_an / 362.25 * corr_fac_an_lb  # [mm/day]
        precip_high = precip_orig_an / 362.25 * corr_fac_an_ub  # [mm/day]
    ds.close()

    # Convert units of precipitation data (if necessary)
    if prec_prod[i]["units_in"] == "mm/day":
        pass
    elif prec_prod[i]["units_in"] == "mm/month":
        days_pm = np.array([31.0, 28.25, 31.0, 30.0, 31.0, 30.0,
                            31.0, 31.0, 30.0, 31.0, 30.0, 31.0])
        prec /= days_pm[:, np.newaxis, np.newaxis]
    elif prec_prod[i]["units_in"] == "m/day":
        prec *= 1000.0
    else:
        raise ValueError("Unknown input units")

    # Load region masks
    ds = xr.open_dataset(path_reg_masks + i + "_region_masks.nc")
    region_masks = {j: ds[j].values.astype(bool) for j in regions}
    ds.close()

    # Compute spatially integrated precipitation for regions
    for j in regions:
        prec_seas = np.empty(12 + 3, dtype=np.float32)
        for k in range(12):
            prec_seas[k] = np.nanmean(prec[k, :, :][region_masks[j]])
        prec_seas[12] = np.nanmean(prec[mask_rainy, :, :].mean(axis=0)
                                   [region_masks[j]])  # rainy season
        prec_seas[13] = np.nanmean(prec[mask_dry, :, :].mean(axis=0)
                                   [region_masks[j]])  # dry season
        prec_seas[14] = np.nanmean(prec.mean(axis=0)[region_masks[j]])
        # yearly mean
        # -> use grid cell area as weights for more accurate spatial averaging!
        prec_prod_reg[i][j] = prec_seas

    # Compute spatially integrated precipitation for regions (-> uncertainty)
    if i == "PBCOR":
        for j in regions:
            prec_an = np.empty(2, dtype=np.float32)
            prec_an[0] = np.nanmean(precip_low[region_masks[j]])
            prec_an[1] = np.nanmean(precip_high[region_masks[j]])
            prec_prod_reg_unc[i][j] = prec_an

# -----------------------------------------------------------------------------
# Temperature
# -----------------------------------------------------------------------------

# List of products
temp_prod = {
    # -------------------------------------------------------------------------
    "APHRODITE": {"file": "/project/pr133/rxiang/data/obs/tmp/APHRO/day/"
                          + "APHRO.2001-2005.025.mon.nc",
                  "var_name": "tave",
                  "units_in": "deg_C"},
    "CRU": {"file": "/project/pr133/rxiang/data/obs/tmp/cru/mo/"
                    + "cru.2001-2005.05.mon.nc",
            "var_name": "tmp",
            "units_in": "deg_C"},
    "ERA5": {"file": "/project/pr133/rxiang/data/era5/pr/mo/"
                     + "era5.mo.2001-2005.mon.nc",
             "var_name": "t2m",
             "units_in": "K"},
    # -------------------------------------------------------------------------
    "CTRL04": {"file": "/project/pr133/rxiang/data/cosmo/EAS04_ctrl/mon/"
                       + "T_2M/2001-2005.T_2M.nc",
               "var_name": "T_2M",
               "units_in": "K"},
    "CTRL11": {"file": "/project/pr133/rxiang/data/cosmo/EAS11_ctrl/mon/"
                       + "T_2M/2001-2005.T_2M.nc",
               "var_name": "T_2M",
               "units_in": "K"},
    # -------------------------------------------------------------------------
    }

# Load topographies
topo = {}
# -----------------------------------------------------------------------------
ds = xr.open_dataset("/project/pr133/csteger/Data/Model/BECCY/remap_topo/"
                     + "MERIT_Eastern_Tibet_remap_APHRO.nc")
topo["APHRODITE"] = ds["Elevation"].values  # [m]
ds.close()
# -----------------------------------------------------------------------------
ds = xr.open_dataset("/project/pr133/csteger/Data/Model/BECCY/remap_topo/"
                     + "MERIT_Eastern_Tibet_remap_CRU.nc")
topo["CRU"] = ds["Elevation"].values  # [m]
ds.close()
# -----------------------------------------------------------------------------
ds = xr.open_dataset("/project/pr133/csteger/Data/Observations/ERA5/"
                     + "ERA5_geopotential_global.nc")
ds = ds.sel(longitude=slice(25.0, 215.0), latitude=slice(77.0, -15.0))
topo["ERA5"] = ds["z"].values.squeeze() / 9.80665  # [m]
ds.close()
# -----------------------------------------------------------------------------
ds = xr.open_dataset("/project/pr133/rxiang/data/extpar/"
                     + "extpar_BECCY_4.4km_merit_unmod_topo.nc")
ds = ds.isel(rlon=slice(30, 680), rlat=slice(30, 680))
topo["CTRL04"] = ds["HSURF"].values
ds.close()
# -----------------------------------------------------------------------------
ds = xr.open_dataset("/project/pr133/rxiang/data/extpar/"
                     + "extpar_EAS_ext_12km_merit_unmod_topo.nc")
ds = ds.isel(rlon=slice(32, 1090), rlat=slice(28, 638))
topo["CTRL11"] = ds["HSURF"].values
# print(ds["rlon"][0].values, ds["rlon"][-1].values)
# print(ds["rlat"][0].values, ds["rlat"][-1].values)
ds.close()
# -----------------------------------------------------------------------------

# Compute monthly temperature for products and regions
temp_prod_reg = {}
temp_prod_elev = {}
elev_bin_edges = np.arange(0.0, 7500.0, 500.0)
elev_bin = elev_bin_edges[:-1] + np.diff(elev_bin_edges) / 2.0
num_bin = len(elev_bin)
for i in temp_prod.keys():

    print("Process product " + i)
    temp_prod_reg[i] = {}

    # Load temperature data
    ds = xr.open_dataset(temp_prod[i]["file"])
    temp = ds[temp_prod[i]["var_name"]].values
    ds.close()

    # Convert units of temperature data (if necessary)
    if temp_prod[i]["units_in"] == "deg_C":
        pass
    elif temp_prod[i]["units_in"] == "K":
        temp -= 273.15
    else:
        raise ValueError("Unknown input units")

    # Load region masks
    ds = xr.open_dataset(path_reg_masks + i + "_region_masks.nc")
    region_masks = {j: ds[j].values.astype(bool) for j in regions}
    ds.close()

    # Compute spatially integrated temperature for regions
    for j in regions:
        temp_seas = np.empty(12 + 2, dtype=np.float32)
        for k in range(12):
            temp_seas[k] = np.nanmean(temp[k, :, :][region_masks[j]])
        temp_seas[12] = np.nanmean(temp[mask_rainy, :, :].mean(axis=0)
                                   [region_masks[j]])  # rainy season (MJJAS)
        temp_seas[13] = np.nanmean(temp[mask_dry, :, :].mean(axis=0)
                                   [region_masks[j]])  # dry season
        # -> use grid cell area as weights for more accurate spatial averaging!
        temp_prod_reg[i][j] = temp_seas

    # Compute temperature lapse rates
    temp_grad_rs = np.empty(num_bin, dtype=np.float32)  # rainy season
    temp_grad_rs.fill(np.nan)
    temp_grad_ds = np.empty(num_bin, dtype=np.float32)  # dry season
    temp_grad_ds.fill(np.nan)
    for j in range(num_bin):
        mask = (topo[i] >= elev_bin_edges[j]) \
               & (topo[i] < elev_bin_edges[j + 1]) & region_masks["HM"]
        if mask.sum() >= 15:
            temp_grad_rs[j] = np.nanmean(temp[mask_rainy, :, :].mean(axis=0)
                                         [mask])
            temp_grad_ds[j] = np.nanmean(temp[mask_dry, :, :].mean(axis=0)
                                         [mask])
    temp_prod_elev[i] = {"rainy_season": temp_grad_rs,
                         "dry_season": temp_grad_ds}

###############################################################################
# Plot
###############################################################################

# Settings
cols_prec = {"APHRODITE": "orange", "GPCC": "red", "IMERG": "forestgreen",
             "ERA5": "plum", "PBCOR": "brown",
             "CTRL04": "royalblue", "CTRL11": "lightskyblue"}
cols_temp = {"APHRODITE": "orange", "CRU": "red",
             "ERA5": "plum",  "CTRL04": "royalblue", "CTRL11": "lightskyblue"}
s = 50
s_small = 30

# Plot
fig = plt.figure(figsize=(12.5, 5.0 + 5.0), dpi=300)
gs = gridspec.GridSpec(3, 8, left=0.1, bottom=0.1, right=0.9, top=0.9,
                       hspace=0.05, wspace=0.0,
                       height_ratios=[1.0, 0.17, 1.0],
                       width_ratios=[0.8, 0.12, 0.06, 0.07,
                                     0.18, 0.18, 0.18, 0.18])
# -----------------------------------------------------------------------------
ax = plt.subplot(gs[0, 0])
for i in prec_prod.keys():
    plt.plot(months_num, prec_prod_reg[i]["ET"][:12], color=cols_prec[i],
             zorder=3)
    plt.scatter(months_num, prec_prod_reg[i]["ET"][:12], s=s_small,
                color=cols_prec[i], label=i, zorder=3)
plt.fill_between(x=[4.5, 9.5], y1=-0.5, y2=15.0, color="black", alpha=0.1)
plt.xticks(months_num, months_char)
plt.text(x=6.0, y=9.3, s="Rainy season", fontsize=10)
plt.ylabel("Precipitation [$mm \; d^{-1}$]", labelpad=5)
plt.axis([0.7, 12.3, 0.0, 10.0])
plt.yticks(np.arange(0, 12, 2), np.arange(0, 12, 2))
plt.title("(a) Eastern Tibet", fontsize=12, fontweight="normal", y=1.01,
          loc="left")
plt.legend(loc="upper left", frameon=False, fontsize=10, ncol=1,
           scatterpoints=1)
# -----------------------------------------------------------------------------
ax = plt.subplot(gs[0, 1:3])
x_3 = np.arange(1, 4)
x_ticks_3 = ["Rainy\nseason", "Dry\nseason", "Annual"]
for i in prec_prod.keys():
    plt.scatter(x_3, prec_prod_reg[i]["ET"][12:], s=s, color=cols_prec[i],
                zorder=3)
plt.fill_between(x=[3.0 - 0.2, 3.0 + 0.2],
                 y1=prec_prod_reg_unc["PBCOR"]["ET"][0],
                 y2=prec_prod_reg_unc["PBCOR"]["ET"][1],
                 color=cols_prec["PBCOR"], facecolor="none", alpha=0.75,
                 zorder=2)
plt.axis([0.5, 3.5, 0.0, 10.0])
plt.xticks(x_3, x_ticks_3, rotation=90, fontsize=9)
plt.yticks(np.arange(0, 12, 2), [""] * 6)
# -----------------------------------------------------------------------------
n = 4
for i in ["HM", "HMC", "HMUN", "HMUS"]:
    ax = plt.subplot(gs[0, n])
    for j in prec_prod.keys():
        plt.scatter(x_3, prec_prod_reg[j][i][12:], s=s, color=cols_prec[j],
                    zorder=3)
    plt.fill_between(x=[3.0 - 0.2, 3.0 + 0.2],
                     y1=prec_prod_reg_unc["PBCOR"][i][0],
                     y2=prec_prod_reg_unc["PBCOR"][i][1],
                     color=cols_prec["PBCOR"], facecolor="none", alpha=0.75,
                     zorder=2)
    plt.text(x=1.40, y=18.1, s=i, fontweight="normal")
    plt.axis([0.5, 3.5, 0.0, 19.5])
    plt.xticks(x_3, x_ticks_3, rotation=90, fontsize=9)
    if n == 4:
        plt.yticks(np.arange(0, 20, 2), np.arange(0, 20, 2))
        plt.title("(b) Hengduan Mountains regions", fontsize=12,
                  fontweight="normal", y=1.01, loc="left")
    else:
        plt.yticks(np.arange(0, 20, 2), [""] * 10)
    n += 1
# -----------------------------------------------------------------------------
ax = plt.subplot(gs[2, 0])
for i in temp_prod.keys():
    plt.plot(months_num, temp_prod_reg[i]["ET"][:12], color=cols_temp[i],
             zorder=3)
    plt.scatter(months_num, temp_prod_reg[i]["ET"][:12], s=s_small,
                color=cols_temp[i], label=i, zorder=3)
plt.fill_between(x=[4.5, 9.5], y1=-2.0, y2=21.0, color="black", alpha=0.1)
plt.xticks(months_num, months_char)
plt.text(x=6.0, y=-0.5, s="Rainy season", fontsize=10)
plt.ylabel("2m temperature [$^{\circ} C$]", labelpad=5)
plt.yticks(np.arange(0, 24, 2), np.arange(0, 24, 2))
plt.axis([0.7, 12.3, -2.0, 21.0])
plt.title("(c) Eastern Tibet", fontsize=12, fontweight="normal", y=1.01,
          loc="left")
plt.legend(loc="upper left", frameon=False, fontsize=10, ncol=1,
           scatterpoints=1)
# -----------------------------------------------------------------------------
ax = plt.subplot(gs[2, 1])
x_2 = np.arange(1, 3)
x_ticks_2 = ["Rainy\nseason", "Dry\nseason"]
for i in temp_prod.keys():
    plt.scatter(x_2, temp_prod_reg[i]["ET"][12:], s=s, color=cols_temp[i],
                zorder=3)
plt.xticks(x_2, x_ticks_2, rotation=90, fontsize=9)
plt.yticks(np.arange(0, 24, 2), [""] * 12)
plt.axis([0.5, 2.5, -2.0, 21.0])
# -----------------------------------------------------------------------------
seas = ("Rainy season", "Dry season")
pos_x = ((4, 6), (6, 8))
for i in range(0, 2):
    ax = plt.subplot(gs[2, pos_x[i][0]:pos_x[i][1]])
    for j in temp_prod_elev.keys():
        k = seas[i].replace(" ", "_").lower()
        plt.plot(temp_prod_elev[j][k], elev_bin, color=cols_temp[j])
        plt.scatter(temp_prod_elev[j][k], elev_bin, s=s_small,
                    color=cols_temp[j], label=j)
    plt.text(x=0.5, y=0.94, s=seas[i], fontsize=10,
             horizontalalignment="center", verticalalignment="center",
             transform=ax.transAxes)
    plt.xticks(np.arange(-20, 30, 10))
    ax.set_xticks(np.arange(-25, 35, 10), minor=True)
    if i == 0:
        plt.yticks(np.arange(0, 7000, 500), np.arange(0, 7.0, 0.5))
        plt.xlabel(" " * 40 + "2m temperature [$^{\circ} C$]")
        plt.ylabel("Elevation [km a.s.l.]", labelpad=7.5)
        plt.title("(d) Hengduan Mountains", fontsize=12, fontweight="normal",
                  y=1.01,
                  loc="left")
    else:
        plt.yticks(np.arange(0, 7000, 500), "" * 14)
    if i == 0:
        plt.xlim([-3.0, 29.0])
    else:
        plt.xlim([-17.5, 20.0])
    plt.ylim([0.0, 6000.0])
# -----------------------------------------------------------------------------
# plt.show()
fig.savefig(path_out + "Precipitation_temperature.png", dpi=300,
            bbox_inches="tight")
# fig.savefig(path_out + "Precipitation_temperature.pdf", bbox_inches="tight")
# plt.close(fig)
