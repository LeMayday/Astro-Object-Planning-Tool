# inspired by https://docs.astropy.org/en/stable/generated/examples/coordinates/plot_obs-planning.html

import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u
from astropy.coordinates import AltAz, EarthLocation, SkyCoord, get_body
from astropy.time import Time

def night_half_duration(local_midnight, location):
    # estimate position of sun to determine local solar noon / midnight
    # number of points used to estimate midnight
    n = 1000
    # first looks for noon over the 13 hours before midnight
    range_to_search = local_midnight + np.linspace(-13, 0, n) * u.hour
    local_frame = AltAz(obstime=range_to_search, location=location)
    sun = get_body('sun', range_to_search)
    sunset_to_midnight = (n - np.argmin(np.abs(sun.transform_to(local_frame).alt.degree))) / n * 13
    return sunset_to_midnight

celestial_object = input("Enter the name of the object: ")
celestial_object_coords = SkyCoord.from_name(celestial_object)

local_coords = input("Enter your latitude and longitude in degrees (separated by a space): ")
lat, long = local_coords.split()
location = EarthLocation(lat=float(lat) * u.deg, lon=float(long) * u.deg, height = 0 * u.m)

date = input("Enter the date of your observation (mm/dd/yyyy): ")
month, day, year = date.split("/")
date = Time(year + "-" + month + "-" + day, format='iso', scale='utc')
local_midnight = Time(date.to_value('jd', 'float') + 1 - float(long) / 360, format='jd')

night_half_length = night_half_duration(local_midnight, location)
delta_midnight = np.linspace(-night_half_length, night_half_length, 1000) * u.hour
local_frame = AltAz(obstime=local_midnight + delta_midnight, location=location)

celestial_object_alt_az = celestial_object_coords.transform_to(local_frame)

fig = plt.figure(figsize = (12, 9))

plt.plot(delta_midnight, celestial_object_alt_az.alt.degree, linewidth=2)
plt.title(celestial_object + " on the night of " + month + "/" + day + "/" + year + " at " + str(lat) + ", " + str(long))
plt.xlabel('Local Solar Time')
plt.ylabel('Altitude (degrees)')
plt.xticks([-night_half_length, 0, night_half_length], ["Sunset", "Midnight", "Sunrise"])
plt.ylim([0, 90])
plt.yticks(np.arange(0, 91, 15))

plt.savefig(celestial_object.replace(" ", "_") + "_altitude_" + month + "_" + day + "_" + year)
