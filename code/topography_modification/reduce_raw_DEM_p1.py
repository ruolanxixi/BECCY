# Description: Compute reduced topography part I (generate reduced elevation
#              topography)
#
# Authors: Ruolan Xiang, Christian R. Steger, IAC ETH Zurich

# Load modules
import xarray as xr
from pyproj import CRS, Transformer
import numpy as np
from scipy import interpolate
import time

###############################################################################
# Settings
###############################################################################

# Terrain reduction settings
lat_0, lon_0 = 33.23000, 95.10000  # reference location [degree]
rad_red = 1800.0  # reduction radius [km]
rad_red_ext = 300.0  # extend reduction radius at flat, amplitude maximum [km]
alpha_0, alpha_1 = -135.0, 45.0  # measured anti-clockwise from East [degree]
fac_amp = 0.90  # 0.75  # amplitude of terrain reduction
topo_min = 500.0  # minimal allowed elevation for terrain reduction

# Map projection
crs_wgs84 = CRS.from_epsg(4326)
crs_aeqd = CRS.from_proj4("+proj=aeqd +lat_0=" + str(lat_0) + " +lon_0="
                          + str(lon_0) + " +datum=WGS84 +units=m")

# DEM tiles
tiles_dem = ("MERIT_N60-N30_E060-E090.nc", "MERIT_N60-N30_E090-E120.nc",
             "MERIT_N30-N00_E060-E090.nc", "MERIT_N30-N00_E090-E120.nc")

# Paths
# path_dem = "/Users/kaktus/Documents/ETH/BECCY/myscripts/topo/"
# path_out = "/Users/kaktus/Documents/ETH/BECCY/myscripts/topo/"
path_dem = "/Users/csteger/Dropbox/IAC/Data/DEMs/MERIT/Tiles/"
path_out = "/Users/csteger/Dropbox/IAC/Data/Model/BECCY/MERIT_tiles_red/"

# Miscellaneous
fac_red_out = True  # output reduction factor

###############################################################################
# Process MERIT data
###############################################################################

# Determine required MERIT domain to load
trans_cart2ellps = Transformer.from_crs(crs_aeqd, crs_wgs84, always_xy=True)
fac_enl = 1.05  # "safety factor" to enlarge domain
ang = np.linspace(0.0, 2.0 * np.pi, 100000)
x_dom = (rad_red + rad_red_ext) * np.cos(ang) * 1000.0 * fac_enl
y_dom = (rad_red + rad_red_ext) * np.sin(ang) * 1000.0 * fac_enl
lon_dom, lat_dom = trans_cart2ellps.transform(x_dom, y_dom)
print("Required rectangular DEM domain: " +
      "longitude: %.2f" % lon_dom.min() + " - %.2f" % lon_dom.max() + " deg, "
      + "latitude: %.2f" % lat_dom.min() + " - %.2f" % lat_dom.max() + " deg")

# Check if DEM tiles cover required domain
ds = xr.open_mfdataset([path_dem + i for i in tiles_dem])
if (ds["lon"].values.min() >= lon_dom.min()
        or ds["lon"].values.max() <= lon_dom.max()
        or ds["lat"].values.min() >= lat_dom.min()
        or ds["lat"].values.max() <= lat_dom.max()):
    raise ValueError("provided DEM tiles do not cover required domain")
ds.close()

# Loop through tiles and process
for i in tiles_dem:

    print((" Process tile " + i + " ").center(60, "#"))

    # Load DEM data
    ds = xr.open_dataset(path_dem + i, mask_and_scale=False)
    lon_tile, lat_tile = ds["lon"].values, ds["lat"].values
    ds = ds.sel(lon=slice(lon_dom.min(), lon_dom.max()),
                lat=slice(lat_dom.max(), lat_dom.min()))
    topo_fill_val = ds["Elevation"]._FillValue
    topo = ds["Elevation"].values  # 16-bit integer
    lon, lat = ds["lon"].values, ds["lat"].values
    ds.close()
    print("Size of DEM data: %.2f" % (topo.nbytes / (10.0 ** 9)) + " GB")
    ind_lon_0 = np.where(lon[0] == lon_tile)[0][0]
    ind_lat_0 = np.where(lat[0] == lat_tile)[0][0]
    del lon_tile, lat_tile

    # -------------------------------------------------------------------------
    # Transform geodetic coordinates to map projection
    # -------------------------------------------------------------------------

    ny, nx = topo.shape
    lon_2d, lat_2d = np.meshgrid(lon, lat)
    trans_ellps2cart = Transformer.from_crs(crs_wgs84, crs_aeqd,
                                            always_xy=True)

    # # Accurate transformation (slow)
    # beg_time = time.time()
    # x, y = trans_ellps2cart.transform(lon_2d, lat_2d)
    # x, y, = x.astype(np.float32), y.astype(np.float32)
    # print("Coordinate transformation: %.2f" % (time.time() - beg_time)
    #       + " sec")
    # del lon_2d, lat_2d

    # Approximation of transformation with bilinear interpolation (fast)
    step = 10  # spacing of accurate transformation
    beg_time = time.time()
    mask = np.zeros((ny, nx), dtype=bool)
    mask[0:ny:step, 0:nx:step] = True
    mask[0:ny:step, -1], mask[-1, 0:nx:step] = True, True
    mask[-1, -1] = True
    print("Fraction of transformed coordinates: %.2f"
          % (mask.sum() / float(mask.size) * 100.0) + " %")
    x_trans, y_trans = trans_ellps2cart.transform(lon_2d[mask], lat_2d[mask])
    shp = mask[:, 0].sum(), mask[0, :].sum()
    x_trans, y_trans = x_trans.reshape(shp), y_trans.reshape(shp)
    ind_y = np.arange(0, lon_2d.shape[0], step)
    if ind_y[-1] != (lon_2d.shape[0] - 1):
        ind_y = np.append(ind_y, lon_2d.shape[0] - 1)
    ind_x = np.arange(0, lon_2d.shape[1], step)
    if ind_x[-1] != (lon_2d.shape[1] - 1):
        ind_x = np.append(ind_x, lon_2d.shape[1] - 1)
    f_ip = interpolate.RectBivariateSpline(ind_y, ind_x, x_trans, kx=1, ky=1)
    x = f_ip(np.arange(lon_2d.shape[0]), np.arange(lon_2d.shape[1])) \
        .astype(np.float32)
    f_ip = interpolate.RectBivariateSpline(ind_y, ind_x, y_trans, kx=1, ky=1)
    y = f_ip(np.arange(lon_2d.shape[0]), np.arange(lon_2d.shape[1])) \
        .astype(np.float32)
    print("Coordinate transformation: %.2f" % (time.time() - beg_time)
          + " sec")
    del lon_2d, lat_2d, mask

    # -------------------------------------------------------------------------

    # Compute distance and azimuth angle
    dist = np.sqrt(x ** 2 + y ** 2) / 1000.0  # [km]
    azim = np.arctan2(y, x)  # [rad]
    del x, y

    # Distance factor
    if rad_red_ext == 0.0:
        fac_dist = np.sin(dist / rad_red * np.pi) ** 2
        fac_dist[dist > rad_red] = 0.0
    else:
        fac_dist = np.zeros((ny, nx), dtype=np.float32)
        mask = (dist <= rad_red / 2.0)
        fac_dist[mask] = np.sin(dist[mask] / rad_red * np.pi) ** 2
        mask = (dist > rad_red / 2.0) & (dist <= rad_red / 2.0 + rad_red_ext)
        fac_dist[mask] = 1.0
        mask = (dist > rad_red / 2.0 + rad_red_ext) \
            & (dist <= rad_red + rad_red_ext)
        fac_dist[mask] = np.sin((dist[mask] - rad_red_ext)
                                / rad_red * np.pi) ** 2
    del dist

    # Azimuth factor (azimuth measured anti-clockwise from East)
    fac_azim = np.sin(((azim - np.deg2rad(alpha_1))
                       / np.deg2rad(alpha_0 - alpha_1)) * np.pi)
    mask = (azim < np.deg2rad(alpha_0)) | (azim > np.deg2rad(alpha_1))
    fac_azim[mask] = 0.0
    del azim, mask

    # Compute total reduction factor
    fac_red = (fac_amp * fac_dist * fac_azim)
    del fac_dist, fac_azim

    # Apply reduction factor to terrain
    mask_water = (topo == topo_fill_val)
    topo = topo.astype(np.float32)
    topo -= ((topo - topo_min) * fac_red).clip(min=0.0)

    # Convert topography back to 16-bit integer
    topo = topo.astype(np.int16)
    if np.any(topo[mask_water] != topo_fill_val):
        print("water grid cell(s) modified -> reset")
        topo[mask_water] = topo_fill_val
    del mask_water

    # Save modified topography in MERIT file
    slic = (slice(ind_lat_0, ind_lat_0 + topo.shape[0]),
            slice(ind_lon_0, ind_lon_0 + topo.shape[1]))
    ds = xr.open_dataset(path_dem + i, mask_and_scale=False)
    ds["Elevation"][slic] = topo
    if fac_red_out:
        fac_red_tile = np.zeros(ds["Elevation"].shape, dtype=np.float32)
        fac_red_tile[slic] = fac_red
        ds["fac_red"] = (("lat", "lon"), fac_red_tile)
    ds.to_netcdf(path_out + i, format="NETCDF4",
                 encoding={"lat": {"_FillValue": None},
                           "lon": {"_FillValue": None}})
    del fac_red
    if fac_red_out:
        del fac_red_tile
