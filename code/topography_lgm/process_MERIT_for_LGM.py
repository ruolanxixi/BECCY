# Description: Process MERIT DEM data for LGM (complement with GEBCO data for
#              ocean grid cells and shift elevation)
#
# Sources for data:
# - MERIT: http://hydro.iis.u-tokyo.ac.jp/~yamadai/MERIT_DEM/
#          -> tiles in NetCDF format process for EXTPAR used
# - GEBCO: https://www.gebco.net/data_and_products/gridded_bathymetry_data/
#          -> GEBCO_2023 Grid (ice surface elevation)
# - PMIP4: https://pmip4.lsce.ipsl.fr/doku.php/exp_design:lgm
#          -> Topography and coastlines -> Peltier ICE-6G-C for PMIP4 ->
#          -> ice_sheet download link ->
#          I6_C.VM5a_10min.0.nc, I6_C.VM5a_10min.21.nc
#
# Authors: Christian R. Steger, IAC ETH Zurich

# Load modules
import os
import numpy as np
import xarray as xr
import subprocess
import time
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.gridspec as gridspec
import scipy
from scipy.interpolate import RegularGridInterpolator
from packaging import version
from utilities.miscellaneous import aggregation_2d, aggregation_1d
from utilities.plot import truncate_colormap

mpl.style.use("classic")

# %matplotlib auto
# %matplotlib auto

# -----------------------------------------------------------------------------
# Settings
# -----------------------------------------------------------------------------

# Relevant MERIT tiles
tiles_merit = (
    # ---------------------------------------------------
    "MERIT_N90-N60_W180-W150", "MERIT_N90-N60_W150-W120",
    "MERIT_N90-N60_W120-W090", "MERIT_N90-N60_W090-W060",
    "MERIT_N90-N60_W060-W030", "MERIT_N90-N60_W030-E000",
    "MERIT_N90-N60_E000-E030", "MERIT_N90-N60_E030-E060",
    "MERIT_N90-N60_E060-E090", "MERIT_N90-N60_E090-E120",
    "MERIT_N90-N60_E120-E150", "MERIT_N90-N60_E150-E180",
    # ---------------------------------------------------
    "MERIT_N60-N30_W180-W150", "MERIT_N60-N30_W150-W120",
    "MERIT_N60-N30_W120-W090", "MERIT_N60-N30_W090-W060",
    "MERIT_N60-N30_W060-W030", "MERIT_N60-N30_W030-E000",
    "MERIT_N60-N30_E000-E030", "MERIT_N60-N30_E030-E060",
    "MERIT_N60-N30_E060-E090", "MERIT_N60-N30_E090-E120",
    "MERIT_N60-N30_E120-E150", "MERIT_N60-N30_E150-E180",
    # ---------------------------------------------------
    "MERIT_N30-N00_W180-W150", "MERIT_N30-N00_W150-W120",
    "MERIT_N30-N00_W120-W090", "MERIT_N30-N00_W090-W060",
    "MERIT_N30-N00_W060-W030", "MERIT_N30-N00_W030-E000",
    "MERIT_N30-N00_E000-E030", "MERIT_N30-N00_E030-E060",
    "MERIT_N30-N00_E060-E090", "MERIT_N30-N00_E090-E120",
    "MERIT_N30-N00_E120-E150", "MERIT_N30-N00_E150-E180",
    # ---------------------------------------------------
    "MERIT_N00-S30_W180-W150", "MERIT_N00-S30_W150-W120",
    "MERIT_N00-S30_W120-W090", "MERIT_N00-S30_W090-W060",
    "MERIT_N00-S30_W060-W030", "MERIT_N00-S30_W030-E000",
    "MERIT_N00-S30_E000-E030", "MERIT_N00-S30_E030-E060",
    "MERIT_N00-S30_E060-E090", "MERIT_N00-S30_E090-E120",
    "MERIT_N00-S30_E120-E150", "MERIT_N00-S30_E150-E180",
    # ---------------------------------------------------
)

# Paths
system = "Desktop"  # "Desktop", "Daint"
if system == "Desktop":
    path_merit = "/Users/csteger/Dropbox/IAC/Data/DEMs/MERIT/Tiles/"
    path_gebco = "/Users/csteger/Dropbox/IAC/Data/DEMs/GEBCO/"
    path_pmip4 = "/Users/csteger/Dropbox/IAC/Data/Model/BECCY/PMIP4/"
    dir_work = "/Users/csteger/Desktop/dir_work/"  # working directory
    dir_out = "/Users/csteger/Desktop/output/"  # output directory
elif system == "Daint":
    path_merit = "/store/c2sm/extpar_raw_data/topo/merit/"
    path_gebco = "/project/pr133/csteger/Data/DEMs/GEBCO/"
    path_pmip4 = "/project/pr133/csteger/Data/Model/BECCY/PMIP4/"
    dir_work = "/scratch/snx3000/csteger/Temp/dir_work/"
    dir_out = "/scratch/snx3000/csteger/Temp/output/"
else:
    raise ValueError("unknown system")

# -----------------------------------------------------------------------------
# Process data
# -----------------------------------------------------------------------------

# Check SciPy version
if version.parse(scipy.__version__) < version.parse("1.10.0"):
    raise ImportError("SciPy version 1.10.0 required. Excessive memory "
                      + "requirements of function 'RegularGridInterpolator'"
                        "in older versions")

# Make GEBCO data periodical along longitude and interpolate elevation at poles
# (-> required for spatial interpolation)
print(" Process GEBCO data ".center(79, "-"))
file_out = dir_work + "GEBCO_2023_periodic_longitude_poles.nc"
if not os.path.isfile(file_out):
    t_beg = time.time()
    ds = xr.open_dataset(path_gebco + "GEBCO_2023.nc")
    print(ds["elevation"].shape)
    # -------------------- make periodical along longitude --------------------
    lon_gebco = ds["lon"].values
    ds = xr.concat((ds.isel(lon=(slice(-1, None))),
                    ds,
                    ds.isel(lon=(slice(0, 1)))), dim="lon")
    ds["lon"].values[0] = lon_gebco[-1] - 360.0
    ds["lon"].values[-1] = lon_gebco[0] + 360.0
    print(np.all(np.diff(ds["lon"].values) > 0.0))
    print(np.max(np.diff(ds["lon"].values))
          - np.min(np.diff(ds["lon"].values)))
    # --------------------- interpolate elevation at poles --------------------
    ds = xr.concat((ds.isel(lat=(slice(0, 1))),  # South Pole
                    ds,
                    ds.isel(lat=(slice(-1, None)))),  # North Pole
                   dim="lat")
    ds["elevation"][0, :].values[:] = ds["elevation"][0, :].values.mean()
    ds["lat"].values[0] = -90.0
    ds["elevation"][-1, :].values[:] = ds["elevation"][-1, :].values.mean()
    ds["lat"].values[-1] = 90.0
    print(np.all(np.diff(ds["lat"].values) > 0.0))
    print(ds["elevation"].shape)
    # -------------------------------------------------------------------------
    extent_gebco = (ds["lon"].values[0], ds["lon"].values[-1],
                    ds["lat"].values[0], ds["lat"].values[-1])
    print("Extent of GEBCO data: "
          + ", ".join(["%.5f" % i for i in extent_gebco]))
    encoding = {"lat": {"_FillValue": None},
                "lon": {"_FillValue": None}}
    ds.to_netcdf(file_out, encoding=encoding)
    print("GEBCO data processed (%.1f" % (time.time() - t_beg) + " s)")
    time.sleep(1.0)
    del ds
else:
    ds = xr.open_dataset(file_out)
    extent_gebco = (ds["lon"].values[0], ds["lon"].values[-1],
                    ds["lat"].values[0], ds["lat"].values[-1])
    ds.close()
    print("Extended GEBCO file already computed")

# Compute delta elevation, make it periodical along longitude and interpolate
# elevation at poles (-> required for spatial interpolation)
print(" Process PMIP elevation delta data ".center(79, "-"))
ds = xr.open_dataset(path_pmip4 + "I6_C.VM5a_10min.21.nc")
topo_delta_pmip = ds["Topo"].values
lon_pmip = ds["lon"].values
lat_pmip = ds["lat"].values
ds.close()
ds = xr.open_dataset(path_pmip4 + "I6_C.VM5a_10min.0.nc")
topo_delta_pmip -= ds["Topo"].values
ds.close()
# ---------------- shift to longitudinal range [-180.0, +180] -----------------
ind_beg = np.where(lon_pmip >= 180.0)[0][0]
lon_pmip = np.append(lon_pmip[ind_beg:] - 360.0, lon_pmip[:ind_beg])
topo_delta_pmip = np.concatenate((topo_delta_pmip[:, ind_beg:],
                                  topo_delta_pmip[:, :ind_beg]), axis=1)
print(topo_delta_pmip.shape)
# ---------------------- make periodical along longitude ----------------------
lon_pmip = np.append(lon_pmip, lon_pmip[:1] + 360)
topo_delta_pmip = np.concatenate((topo_delta_pmip,
                                  topo_delta_pmip[:, :1]), axis=1)
print(np.all(np.diff(lon_pmip) > 0.0))
print(np.max(np.diff(lon_pmip)) - np.min(np.diff(lon_pmip)))
# ----------------------- interpolate elevation at poles ----------------------
lat_pmip = np.concatenate((np.array([-90.0], dtype=np.float32),
                           lat_pmip,
                           np.array([90.0], dtype=np.float32)))
len_lon = topo_delta_pmip.shape[1]
topo_delta_pmip = np.concatenate(
    (np.repeat(topo_delta_pmip[0, :].mean(), len_lon).reshape(1, len_lon),
     topo_delta_pmip,
     np.repeat(topo_delta_pmip[-1, :].mean(), len_lon).reshape(1, len_lon)),
    axis=0)
print(np.all(np.diff(lat_pmip) > 0.0))
print(topo_delta_pmip.shape)
# -----------------------------------------------------------------------------
extent_pmip = (lon_pmip[0], lon_pmip[-1], lat_pmip[0], lat_pmip[-1])
print("Extent of PMIP delta elevation data: "
      + ", ".join(["%.5f" % i for i in extent_pmip]))

# # Test plot
# levels = np.arange(-140.0, 160.0, 20.0)
# cmap = plt.get_cmap("Spectral")
# norm = mpl.colors.BoundaryNorm(levels, ncolors=cmap.N, extend="both")
# plt.figure(figsize=(10, 7))
# plt.pcolormesh(lon_pmip, lat_pmip, topo_delta_pmip, cmap=cmap, norm=norm)
# plt.axis([-180.0, 180.0, -90.0, 90.0])
# plt.colorbar()

# Loop through MERIT tiles and process
for i in tiles_merit:

    print((" Process MERIT tile " + i + " ").center(79, "-"))

    # Unzip MERIT data (only required for Desktop data)
    if system == "Desktop":
        t_beg = time.time()
        cmd = "gunzip -c"
        sf = path_merit + i + ".nc.xz"
        tf = dir_work + i + ".nc"
        subprocess.call(cmd + " " + sf + " > " + tf, shell=True)
        print("File unzipped (%.1f" % (time.time() - t_beg) + " s)")
    else:
        tf = path_merit + i + ".nc"

    # Open MERIT data
    t_beg = time.time()
    ds = xr.open_dataset(tf)
    elevation_merit = ds["Elevation"].values  # ~5.2 GB
    lon_merit = ds["lon"].values
    lat_merit = ds["lat"].values
    ds.close()
    if ((lon_merit[0] < extent_gebco[0])
            or (lon_merit[-1] > extent_gebco[1])
            or (lat_merit[-1] < extent_gebco[2])
            or (lat_merit[0] > extent_gebco[3])):
        raise ValueError("MERIT domain exceeds GEBCO domain")
    if ((lon_merit[0] < extent_pmip[0])
            or (lon_merit[-1] > extent_pmip[1])
            or (lat_merit[-1] < extent_pmip[2])
            or (lat_merit[0] > extent_pmip[3])):
        raise ValueError("MERIT domain exceeds PMIP domain")
    print("File opened (%.1f" % (time.time() - t_beg) + " s)")

    # Open required GEBCO subdomain
    ds = xr.open_dataset(dir_work + "GEBCO_2023_periodic_longitude_poles.nc")
    ds = ds.sel(lon=slice(lon_merit[0] - 0.02, lon_merit[-1] + 0.02),
                lat=slice(lat_merit[-1] - 0.02, lat_merit[0] + 0.02))
    elevation_gebco = ds["elevation"].values
    lon_gebco = ds["lon"].values  # increasing [degree]
    lat_gebco = ds["lat"].values  # increasing [degree]
    ds.close()

    # Interpolation GEBCO bathymetry data to MERIT (bilinear)
    t_beg = time.time()
    f_ip = RegularGridInterpolator((lat_gebco, lon_gebco),
                                   elevation_gebco, method="linear",
                                   bounds_error=True)
    mask_ocean = np.isnan(elevation_merit)  # ~1.3 GB
    frac_ocean = mask_ocean.sum() / mask_ocean.size * 100.0  # [%]
    print("Ocean fraction: %.1f" % frac_ocean + " %")
    # -------------------------------------------------------------------------
    # Interpolate in blocks of 6000 x 6000 due to large memory consumption
    # -------------------------------------------------------------------------
    for j in range(0, 36000, 6000):
        for k in range(0, 36000, 6000):
            slic = (slice(j, j + 6000), slice(k, k + 6000))
            ind_0, ind_1 = np.where(mask_ocean[slic])
            lon_ip, lat_ip = np.meshgrid(lon_merit[slic[1]],
                                         lat_merit[slic[0]])
            lon_ip = lon_ip[mask_ocean[slic]]
            lat_ip = lat_ip[mask_ocean[slic]]
            pts_ip = np.vstack((lat_ip, lon_ip)).transpose()
            elevation_merit[slic][mask_ocean[slic]] = f_ip(pts_ip)
    # -------------------------------------------------------------------------
    print("GEBCO data interpolated (%.1f" % (time.time() - t_beg) + " s)")

    # Check minimal/maximal elevation
    elev_min = np.min(elevation_merit)
    elev_max = np.max(elevation_merit)
    print("Minimal/maximal elevation: %.1f" % elev_min + ", %.1f" % elev_max
          + " m")

    # Plot for visual check
    t_beg = time.time()
    elevation_agg = aggregation_2d(elevation_merit, 10, 10, "mean")
    lon_agg = aggregation_1d(lon_merit, 10, "mean")
    lat_agg = aggregation_1d(lat_merit, 10, "mean")
    # -------------------------------------------------------------------------
    fig = plt.figure(figsize=(12, 6))
    gs = gridspec.GridSpec(2, 2, left=0.1, bottom=0.1, right=0.9,
                           top=0.9, hspace=0.06, wspace=0.06,
                           height_ratios=[1.0, 0.05])
    # -------------------------------------------------------------------------
    levels = np.arange(-6000.0, 500.0, 500.0)
    cmap = plt.get_cmap("YlGnBu_r")
    norm = mpl.colors.BoundaryNorm(levels, ncolors=cmap.N, extend="min")
    ax = plt.subplot(gs[0, 0], projection=ccrs.PlateCarree())
    ax.set_facecolor("grey")
    data_plot = np.ma.masked_where(elevation_agg > 0.0, elevation_agg)
    plt.pcolormesh(lon_agg, lat_agg, data_plot, cmap=cmap, norm=norm)
    ax.set_aspect("auto")
    gl = ax.gridlines(crs=ccrs.PlateCarree(), linewidth=0.5, color="black",
                      alpha=0.8, linestyle="-", draw_labels=True)
    gl.bottom_labels = False
    gl.right_labels = False
    ax.set_extent([lon_agg.min(), lon_agg.max(),
                   lat_agg.min(), lat_agg.max()], crs=ccrs.PlateCarree())
    # -------------------------------------------------------------------------
    ax = plt.subplot(gs[1, 0])
    mpl.colorbar.ColorbarBase(ax, cmap=cmap, norm=norm,
                              orientation="horizontal")
    plt.xlabel("Depth [m]")
    # -------------------------------------------------------------------------
    levels = np.arange(0.0, 6500.0, 500.0)
    cmap = truncate_colormap(plt.get_cmap("terrain"), 0.35, 1.0)
    norm = mpl.colors.BoundaryNorm(levels, ncolors=cmap.N, extend="min")
    ax = plt.subplot(gs[0, 1], projection=ccrs.PlateCarree())
    ax.set_facecolor("grey")
    data_plot = np.ma.masked_where(elevation_agg < 0.0, elevation_agg)
    plt.pcolormesh(lon_agg, lat_agg, data_plot, cmap=cmap, norm=norm)
    ax.set_aspect("auto")
    gl = ax.gridlines(crs=ccrs.PlateCarree(), linewidth=0.5, color="black",
                      alpha=0.8, linestyle="-", draw_labels=True)
    gl.bottom_labels = False
    gl.left_labels = False
    gl.right_labels = False
    ax.set_extent([lon_agg.min(), lon_agg.max(),
                   lat_agg.min(), lat_agg.max()], crs=ccrs.PlateCarree())
    # -------------------------------------------------------------------------
    ax = plt.subplot(gs[1, 1])
    mpl.colorbar.ColorbarBase(ax, cmap=cmap, norm=norm,
                              orientation="horizontal")
    plt.xlabel("Elevation [m]")
    # -------------------------------------------------------------------------
    fig.savefig(dir_out + i + "_merged.png", dpi=300,
                bbox_inches="tight")
    plt.close(fig)
    print("Plot created (%.1f" % (time.time() - t_beg) + " s)")

    # Interpolation elevation delta to merged DEM (bilinear)
    t_beg = time.time()
    f_ip = RegularGridInterpolator((lat_pmip.astype(np.float64),
                                    lon_pmip.astype(np.float64)),
                                   topo_delta_pmip.astype(np.float64),
                                   method="linear", bounds_error=True)
    # -------------------------------------------------------------------------
    # Interpolate in blocks of 6000 x 6000 due to large memory consumption
    # -------------------------------------------------------------------------
    for j in range(0, 36000, 6000):
        for k in range(0, 36000, 6000):
            slic = (slice(j, j + 6000), slice(k, k + 6000))
            lon_ip, lat_ip = np.meshgrid(lon_merit[slic[1]],
                                         lat_merit[slic[0]])
            pts_ip = np.vstack((lat_ip.ravel(), lon_ip.ravel())).transpose()
            elevation_merit[slic] += f_ip(pts_ip).reshape(6000, 6000)
    # -------------------------------------------------------------------------
    print("Delta elevation interpolated (%.1f" % (time.time() - t_beg) + " s)")

    # Plot for visual check
    t_beg = time.time()
    elevation_agg = aggregation_2d(elevation_merit, 10, 10, "mean")
    # -------------------------------------------------------------------------
    fig = plt.figure(figsize=(6, 6))
    gs = gridspec.GridSpec(2, 1, left=0.1, bottom=0.1, right=0.9,
                           top=0.9, hspace=0.06, wspace=0.06,
                           height_ratios=[1.0, 0.05])
    # -------------------------------------------------------------------------
    levels = np.arange(0.0, 6500.0, 500.0)
    cmap = truncate_colormap(plt.get_cmap("terrain"), 0.35, 1.0)
    norm = mpl.colors.BoundaryNorm(levels, ncolors=cmap.N, extend="min")
    ax = plt.subplot(gs[0, 0], projection=ccrs.PlateCarree())
    ax.set_facecolor("grey")
    data_plot = np.ma.masked_where(elevation_agg < 0.0, elevation_agg)
    plt.pcolormesh(lon_agg, lat_agg, data_plot, cmap=cmap, norm=norm)
    ax.set_aspect("auto")
    gl = ax.gridlines(crs=ccrs.PlateCarree(), linewidth=0.5, color="black",
                      alpha=0.8, linestyle="-", draw_labels=True)
    gl.bottom_labels = False
    gl.right_labels = False
    ax.set_extent([lon_agg.min(), lon_agg.max(),
                   lat_agg.min(), lat_agg.max()], crs=ccrs.PlateCarree())
    # -------------------------------------------------------------------------
    ax = plt.subplot(gs[1, 0])
    mpl.colorbar.ColorbarBase(ax, cmap=cmap, norm=norm,
                              orientation="horizontal")
    plt.xlabel("Elevation [m]")
    # -------------------------------------------------------------------------
    fig.savefig(dir_out + i + "_LGM.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("Plot created (%.1f" % (time.time() - t_beg) + " s)")

    # Save modified MERIT DEM
    t_beg = time.time()
    ds = xr.open_dataset(tf)
    elevation_merit[(elevation_merit < 0.0) & mask_ocean] = np.nan
    # only set grid cells to ocean that are already ocean in present-day
    # -> keep land areas below sea level (like the Dead Sea)
    ds["Elevation"].values = elevation_merit
    encoding = {"lat": {"_FillValue": None},
                "lon": {"_FillValue": None}}
    ds.to_netcdf(dir_out + i + ".nc", encoding=encoding)
    print("Modified MERIT tile saved (%.1f" % (time.time() - t_beg) + " s)")

    # Remove unzipped raw MERIT file
    if system == "Desktop":
        time.sleep(1.0)
        os.remove(tf)
