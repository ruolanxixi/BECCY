# Description: Test computation of isostatic adjustment
#
# Author: Christian Steger, April 2022

# Load modules
import os
import sys
import gflex
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
from geographiclib.geodesic import Geodesic
from scipy.special import kei
import time

mpl.style.use("classic")

# Load required functions
sys.path.append("/Users/csteger/Downloads/BECCY/code/")
from auxiliary import gridcoord
sys.path.append("/Users/csteger/Downloads/BECCY/code/Isostasy/")
from isostasy_cy import deflection_xy, deflection_lonlat

###############################################################################
# Reference solution from gFlex
###############################################################################

flex = gflex.F2D()

flex.Quiet = False
flex.Method = "FD"
# Solution method: * FD (finite difference),
# SAS (superposition of analytical solutions), SAS_NG (ungridded SAS)
flex.PlateSolutionType = "vWC1994"
# van Wees and Cloetingh (1994); other option is 'G2009': Govers et al. (2009)
flex.Solver = "direct"  # direct or iterative

flex.g = 9.81  # acceleration due to gravity [m s-2]
flex.E = 65E9  # Young's modulus [Pa]
flex.nu = 0.25  # Poisson's ratio [-]
flex.rho_m = 3400.0  # mantle density  [kg m-3]
flex.rho_fill = 0.0  # infill material density  [kg m-3]

flex.Te = 35000.0 * np.ones((50, 50))  # Elastic thickness [m]
flex.qs = np.zeros((50, 50))  # surface load stresses [Pa]
flex.qs[10:40, 10:40] += 1E6
flex.dx = 5000.0  # grid cell size, x-oriented [m]
flex.dy = 5000.0  # grid cell size, y-oriented [m]

# Boundary conditions can be:
# (FD): 0Slope0Shear, 0Moment0Shear, 0Displacement0Slope, Mirror, or Periodic
# For SAS or SAS_NG, NoOutsideLoads is valid, and no entry defaults to this
flex.BC_W = "0Moment0Shear"  # "'0Displacement0Slope' # west boundary condition
flex.BC_E = "0Moment0Shear"  # '0Moment0Shear' # east boundary condition
flex.BC_S = "0Moment0Shear"  # '0Displacement0Slope' # south boundary condition
flex.BC_N = "0Moment0Shear"  # '0Displacement0Slope' # north boundary condition

# latitude/longitude solutions are exact for SAS, approximate otherwise
# latlon = # true/false: flag to enable lat/lon input. Defaults False.
# PlanetaryRadius = # radius of planet [m], for lat/lon solutions

flex.initialize()
flex.run()
flex.finalize()

# If you want to plot the output
flex.plotChoice = "both"
# An output file for deflections could also be defined here
# flex.wOutFile =
flex.output()
# Plots and/or saves output, or does nothing, depending on whether
# flex.plotChoice and/or flex.wOutFile have been set

deflection_ar = flex.w


###############################################################################
# Functions
###############################################################################

# Great circle distance
def great_circ_dist(lat_1, lon_1, lat_2, lon_2):
    rad_earth = 6370997.0  # earth radius (PROJ default value) [m]
    term = np.sin(lat_1) * np.sin(lat_2) + np.cos(lat_1) * np.cos(lat_2)\
        * np.cos(lon_2 - lon_1)
    if np.isscalar(term):
        if term > 1.0:
            term = 1.0
    else:
        term = term.clip(max=1.0)
    return rad_earth * np.arccos(term)


# -----------------------------------------------------------------------------
# Test function
# -----------------------------------------------------------------------------

lat_1, lon_1 = 39.0, 95.2
lat_2, lon_2 = 39.0, 90.2

geodesic = Geodesic.WGS84.Inverse(lat_1, lon_1, lat_2, lon_2)["s12"]
gc_dist = great_circ_dist(np.deg2rad(lat_1), np.deg2rad(lon_1),
                          np.deg2rad(lat_2), np.deg2rad(lon_2))
print(gc_dist / 1000.0, geodesic / 1000.0)
print((gc_dist - geodesic) / 1000.0)

###############################################################################
# Own implementation
###############################################################################

# Settings
rho_m = 3400.0  # mantle density  [kg m-3]
rho_fill = 0.0  # infill material density (density of air: ~1.2) [kg m-3]
g = 9.81  # acceleration due to gravity [m s-2]
Te = 35000.0  # Elastic thickness [m]
E = 65E9  # Young's modulus [Pa]
nu = 0.25  # Poisson's ratio [-]

# Compute parameters
d = E * Te ** 3 / (12.0 * (1.0 - nu ** 2))  # Flexural rigidity [m2 kg s-2]
alpha = (d / ((rho_m - rho_fill) * g)) ** 0.25  # 2D flexural parameter [m]

# -----------------------------------------------------------------------------
# Reproduce above gFlex output
# -----------------------------------------------------------------------------

# Grid
ny, nx = 50, 50
x_grid = np.linspace(0.0, 250000.0, nx + 1)  # [m]
y_grid = np.linspace(0.0, 250000.0, ny + 1)  # [m]
x = x_grid[:-1] + np.diff(x_grid) / 2.0
y = y_grid[:-1] + np.diff(y_grid) / 2.0
dx, dy = np.diff(x).mean(), np.diff(y).mean()

# Point loads
dh = 30.89  # [m] elevation difference
q = np.zeros((ny, nx))  # point loads [N]
q[10:40, 10:40] += (dh * dy * dy) * rho_m * g

# Python
time_beg = time.time()
x_2d, y_2d = np.meshgrid(x, y)
w = np.empty_like(q)
for i in range(ny):
    for j in range(nx):
        r = np.sqrt((x[j] - x_2d) ** 2 + (y[i] - y_2d) ** 2)
        w[i, j] = (q * alpha ** 2 / (2 * np.pi * d)
                   * kei(r / alpha)).sum()
print("Elapsed time (Python): " + "%.2f" % (time.time() - time_beg) + " sec")

# Cython (planar grid; default)
time_beg = time.time()
w_cy = deflection_xy(x, y, q, alpha, d)
print("Elapsed time (Cython; default): " + "%.2f" % (time.time() - time_beg)
      + " sec")
print(np.abs(w_cy - w).max())

# Cython (planar grid; fast)
time_beg = time.time()
w_cy = deflection_xy(x, y, q, alpha, d, "fast")
print("Elapsed time (Cython; fast): " + "%.2f" % (time.time() - time_beg)
      + " sec")
print(np.abs(w_cy - w).max())

# Test plot
plt.figure()
plt.pcolormesh(x / 1000.0, y / 1000.0, w, shading="auto")
plt.colorbar()

load_thick = q / (rho_m * g * dx * dy)  # [m]
print("Maximal load thickness: %.2f" % load_thick.max())

# -----------------------------------------------------------------------------
# Compare solution on planar/spherical grid
# -----------------------------------------------------------------------------

# Grid
# ny, nx = 51, 51
ny, nx = 501, 501  # ca. 65 sec.
# ny, nx = 1001, 1001
x = np.linspace(-250000.0, 250000.0, nx)  # [m]
y = np.linspace(-250000.0, 250000.0, ny)  # [m]
dx, dy = np.diff(x).mean(), np.diff(y).mean()

# Point loads
dh = 1000.0  # [m] elevation difference
q = np.zeros((ny, nx))  # point loads [N]
rho_nsr = 2600.0  # density of near-surface rock [kg m-3]
# -> http://hyperphysics.phy-astr.gsu.edu/hbase/Geophys/earthstruct.html
q[100:400, 100:400] += (dh * dy * dy) * rho_nsr * g
# q[10:40, 10:40] += (dh * dy * dy) * rho_nsr * g

# Cython (planar grid)
time_beg = time.time()
w = deflection_xy(x, y, q, alpha, d, "fast")
print("Elapsed time (Cython): " + "%.2f" % (time.time() - time_beg) + " sec")

# Test plot
plt.figure()
plt.pcolormesh(x / 1000.0, y / 1000.0, w, shading="auto")
plt.title("Planar grid", fontweight="bold", fontsize=12)
plt.xlabel("x [km]")
plt.ylabel("y [km]")
plt.colorbar()

# Spherical grid
rad_earth = 6370997.0  # earth radius (PROJ default value) [m]
spac = 360.0 / (2.0 * np.pi * rad_earth) * dx  # [deg]
lon = np.deg2rad(np.linspace(-(nx // 2) * spac, (nx // 2) * spac, nx))  # [rad]
lat = np.deg2rad(np.linspace(-(ny // 2) * spac, (ny // 2) * spac, ny))  # [rad]

# Compute surface area of grid cells
lon_grid, lat_grid = gridcoord(lon, lat)
d_lon = np.diff(lon).mean()
if np.diff(lat[:2])[0] < 0.0:  # latitude in decreasing order
    gc_area = rad_earth ** 2 * d_lon * \
              (np.sin(lat_grid[:-1]) - np.sin(lat_grid[1:]))
else:
    gc_area = rad_earth ** 2 * d_lon * \
              (np.sin(lat_grid[1:]) - np.sin(lat_grid[:-1]))
gc_area = np.repeat(gc_area[:, np.newaxis], len(lon), axis=1)  # [m2]

# # Python
# time_beg = time.time()
# lon_2d, lat_2d = np.meshgrid(lon, lat)
# w = np.empty_like(q)
# for i in range(ny):
#     for j in range(nx):
#         r = great_circ_dist(lat[i], lon[j], lat_2d, lon_2d)
#         w[i, j] = (q * gc_area * alpha ** 2 / (2 * np.pi * d)
#                    * kei(r / alpha)).sum()
# print("Elapsed time (Python): " + "%.2f" % (time.time() - time_beg) + " sec")

# Cython (spherical grid)
time_beg = time.time()
w_cy = deflection_lonlat(lon, lat, (q * gc_area), alpha, d, "fast")
print("Elapsed time (Cython): " + "%.2f" % (time.time() - time_beg) + " sec")

# Test plot
plt.figure()
plt.pcolormesh(np.rad2deg(lon), np.rad2deg(lat), w, shading="auto")
plt.title("Spherical grid", fontweight="bold", fontsize=12)
plt.xlabel("Longitude [deg]")
plt.ylabel("Latitude [deg]")
plt.colorbar()

###############################################################################
# Approximate kei function by linear interpolation
# (-> tested; only leads to a ca. 3x acceleration in Cython...)
###############################################################################

r = np.linspace(0.0, 1100.0 * 1000.0, 50001)  # [m]
val_kei = kei(r / alpha)
print("Last value: %.8f" % val_kei[-1])
print("Value for 21.0: %.8f" % kei(21.0))

# Test plot
plt.figure()
plt.scatter((r / alpha), val_kei)
plt.plot((r / alpha), val_kei)
plt.xlabel("r / alpha")
plt.ylabel("kei(r / alpha)")

r_betw = r[:-1] + np.diff(r / 2.0)

val_kei_exact = kei(r_betw / alpha)
val_kei_approx = np.interp(r_betw / alpha, (r / alpha), val_kei)
err_abs_max = np.abs(val_kei_approx - val_kei_exact).max()
print("Maximal absolute error: %.8f" % err_abs_max)

# %timeit -r 5 -n 5 kei(r_betw / alpha)
# %timeit -r 5 -n 5 np.interp(r_betw / alpha, (r / alpha), val_kei)

plt.figure()
plt.plot((r_betw / alpha), val_kei_exact, color="black", lw=1.5)
plt.plot((r_betw / alpha), val_kei_approx, color="red", lw=1.5, ls="--")
plt.xlabel("r / alpha")
plt.ylabel("kei(r / alpha)")

# -----------------------------------------------------------------------------
# Code for Cython implementation
# -----------------------------------------------------------------------------

num = 50000
r_alpha_max = 21.0

spac = r_alpha_max / float(num)
r_alpha = np.empty(num + 1, dtype=np.float64)
kei_val = np.empty(num + 1, dtype=np.float64)
for i in range(num + 1):
    r_alpha[i] = float(i) * spac
    kei_val[i] = kei(r_alpha[i])


r_alpha_ip = 1500.0 * 1000.0 / alpha
if r_alpha_ip >= r_alpha_max:
    kei_val_ip = kei_val[num]
else:
    ind = int(np.floor(r_alpha_ip / spac))
    slope = (kei_val[ind + 1] - kei_val[ind]) \
        / (r_alpha[ind + 1] - r_alpha[ind])
    kei_val_ip = slope * (r_alpha_ip - r_alpha[ind]) + kei_val[ind]

print(kei_val_ip)
print(np.interp(r_alpha_ip, r_alpha, kei_val))
print(kei_val_ip - np.interp(r_alpha_ip, r_alpha, kei_val))

plt.figure()
plt.plot(r_alpha, kei_val, color="black", lw=1.5)
plt.xlabel("r / alpha")
plt.ylabel("kei(r / alpha)")
