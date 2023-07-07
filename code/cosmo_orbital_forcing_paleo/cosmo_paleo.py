# Description: Check orbital parameter implementation in paleo-COSMO
#
# Author: Christian R. Steger, IAC ETH Zurich, July 2023

# Load modules
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import datetime as dt
from skyfield.api import load
from skyfield.framelib import itrs

mpl.style.use("classic")

# %matplotlib auto

# Load required functions
sys.path.append("/Users/csteger/Desktop/cosmo_orbital_forcing_paleo/")
from cosmo_paleo import orbit

###############################################################################
# Compare seasonal cycle
###############################################################################

# # present day values as of 1850
# recc = 0.016764  # Eccentricity for paleo orbit
# robld = 23.4593  # Obliquity for paleo orbit
# rlonp = 280.327  # Longitude of perihelion

# LGM (-> https://pmip4.lsce.ipsl.fr/doku.php/exp_design:lgm)
recc = 0.018994
robld = 22.949
rlonp = 294.42  # (114.42 + 180.0)

jj = 2013
itaja = np.arange(1, 366)
zstunde = 12.0
itype_calendar = 0

zdek = np.empty(itaja.size)
zsocof = np.empty(itaja.size)
zsct_save = np.empty(itaja.size)
zeit0 = np.empty(itaja.size)

# Plot: Solar declination angle
plt.figure(figsize=(9, 5))
for ii, i in enumerate(itaja):
    zdek[ii], zsocof[ii], zsct_save[ii], zeit0[ii] \
        = orbit(jj, i, zstunde, False, itype_calendar, recc, robld, rlonp)
plt.plot(itaja, np.rad2deg(zdek), color="red", label="present day")
for ii, i in enumerate(itaja):
    zdek[ii], zsocof[ii], zsct_save[ii], zeit0[ii] \
        = orbit(jj, i, zstunde, True, itype_calendar, recc, robld, rlonp)
plt.plot(itaja, np.rad2deg(zdek), color="blue", label="LGM")
plt.xlabel("Day of the year")
plt.ylabel("Declination angle of the sun [degree]")
plt.axis([0.0, 366.0, -26.0, 26.0])
plt.legend(frameon=False, fontsize=11)

# Plot: Solar 'constant'
plt.figure(figsize=(9, 5))
for ii, i in enumerate(itaja):
    zdek[ii], zsocof[ii], zsct_save[ii], zeit0[ii] \
        = orbit(jj, i, zstunde, False, itype_calendar, recc, robld, rlonp)
plt.plot(itaja, zsct_save, color="red", label="present day")
for ii, i in enumerate(itaja):
    zdek[ii], zsocof[ii], zsct_save[ii], zeit0[ii] \
        = orbit(jj, i, zstunde, True, itype_calendar, recc, robld, rlonp)
plt.plot(itaja, zsct_save, color="blue", label="LGM")
plt.hlines(y=1368.0, xmin=0.0, xmax=366.0, color="black", lw=1.2)
plt.xlabel("Day of the year")
plt.ylabel("Solar 'constant' [W m-2]")
plt.xlim([0.0, 366.0])
plt.legend(frameon=False, fontsize=11)

###############################################################################
# Compare diurnal cycle
###############################################################################

year = 2013
month = 4
day = 15
doy = (dt.datetime(year, month, day) - dt.datetime(year, 1, 1)).days + 1
hour = np.arange(0.0, 24.0, 0.25)

zdek = np.empty(hour.size)
zsocof = np.empty(hour.size)
zsct_save = np.empty(hour.size)
zeit0 = np.empty(hour.size)

# Skyfield
planets = load("de421.bsp")
sun = planets["sun"]
earth = planets["earth"]

# Plot: Subsolar longitude
plt.figure(figsize=(9, 5))
# -----------------------------------------------------------------------------
ts = load.timescale()
data_sf = np.empty((2, hour.size), dtype=np.float32)  # [degree]
for ii, i in enumerate(hour):
    hour_int = int(np.floor(i))
    minute = int((i - hour_int) * 60.0)
    ta = dt.datetime(year, month, day, hour_int, minute,
                     tzinfo=dt.timezone.utc)
    t = ts.from_datetime(ta)
    position = earth.at(t).observe(sun)
    x, y, z = position.frame_xyz(itrs).au  # unit: astronomical units
    # ITRS reference frame: positions are fixed relative to the Earth’s surface
    dist = np.sqrt(x ** 2 + y ** 2 + z ** 2)
    data_sf[0, ii] = np.rad2deg(np.arcsin(z / dist))
    data_sf[1, ii] = np.rad2deg(np.arctan2(y, x))
plt.plot(hour, data_sf[1, :], color="lightgrey", lw=6.0,
         label="present day (Skyfield)")
# -----------------------------------------------------------------------------
for ii, i in enumerate(hour):
    zdek[ii], zsocof[ii], zsct_save[ii], zeit0[ii] \
        = orbit(year, doy, i, False, itype_calendar, recc, robld, rlonp)
plt.plot(hour, np.rad2deg(zeit0) * (-1.0), color="red", label="present day")
# -----------------------------------------------------------------------------
for ii, i in enumerate(hour):
    zdek[ii], zsocof[ii], zsct_save[ii], zeit0[ii] \
        = orbit(year, doy, i, True, itype_calendar, recc, robld, rlonp)
plt.plot(hour, np.rad2deg(zeit0) * (-1.0), color="blue", ls="--",
         label="LGM")
# -----------------------------------------------------------------------------
plt.xlabel("Hour (UTC)")
plt.ylabel("Subsolar longitude [degree]")
plt.legend(frameon=False, fontsize=11)

# Notes:
# - Approximation: the same 'equation of time' is used for paleo-climate as
#   for the present day
# - Solar declination angle is constant over a day

