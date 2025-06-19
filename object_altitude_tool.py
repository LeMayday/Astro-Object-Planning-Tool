# Author: Mayday
# Date: 6/29/2024
# inspired by https://docs.astropy.org/en/stable/generated/examples/coordinates/plot_obs-planning.html
# requires matplotlib, numpy, astropy, and jplephem modules

from typing import Tuple
import argparse
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 16})
import numpy as np
import astropy.units as u
from astropy.coordinates import AltAz, EarthLocation, SkyCoord, get_body
from astropy.time import Time
from plot_utils import colored_line, colored_line_between_pts, polar_subplots, make_proxy

class Observing_Metrics:
    '''
    Class object that calculates and stores observing metrics for the night
    time_above_horizon
    time_above_15_deg
    time_astro_dark (also above horizon)
    time_optimal (time above 15 deg and during astro dark)
    '''
    def __init__(self, object_alt_az: SkyCoord, local_midnight: Time, location: EarthLocation):
        # # num data points
        # n = np.size(object_alt_az)
        # # fraction of altitudes that match criterion * length of night in hours
        # self.time_above_horizon = np.count_nonzero(object_alt_az.alt.degree > 0) / n * night_half_length * 2
        # self.time_above_15_deg = np.count_nonzero(object_alt_az.alt.degree > 15) / n * night_half_length * 2
        # # fraction of night before astro dark * n to get start index
        # idx1 = np.floor((night_half_length - astro_dark_half_length) / (night_half_length * 2) * n).astype(int)
        # # fraction of night until end of astro dark * n to get end index
        # idx2 = np.ceil((night_half_length + astro_dark_half_length) / (night_half_length * 2) * n).astype(int)
        # co_alt_az_astro_dark = object_alt_az.alt.degree[idx1 : idx2]
        # self.time_astro_dark = np.count_nonzero(co_alt_az_astro_dark > 0) / n * night_half_length * 2
        # # TODO: add time_optimal
        pass

def get_night_start_end(local_midnight: Time, location: EarthLocation) -> Tuple[float, float]:
    '''
    Determines time from sunset to local midnight and astronomical dark to local midnight
    Inputs: Time object with JD format, EarthLocation based on lat, long
    Outputs: time from sunset to midnight and time from astro dark to midnight in hours
    '''
    # number of points used to estimate sunset
    n = 1000
    # given midnight, we know sunset will occur sometime within the 12 hours prior
    # for simplicity, sunset to midnight and midnight to sunrise is assumed to be the same
    range_to_search = local_midnight + np.linspace(-12, 12, n) * u.hour
    sun_alt_az = get_object_alt_az('sun', range_to_search, location)
    night = np.nonzero(sun_alt_az.alt.degree <= 0)
    night_start = night[0] / n * 24 - 12
    night_end = night[-1] / n * 24 - 12
    # astronomical dark occurs when the sun is 18 degrees below the horizon
    astro_dark = np.nonzero(sun_alt_az.alt.degree <= -18)
    astro_dark_start = astro_dark[0] / n * 24 - 12
    astro_dark_end = astro_dark[-1] / n * 24 - 12
    return Tuple(night_start, night_end), Tuple(astro_dark_start, astro_dark_end)

def get_local_midnight(date: Time, long: u.quantity.Quantity = 0*u.deg) -> Time:
    # +1 for night of date, - long / 360 to convert to local midnight
    local_midnight = Time(date.to_value('jd', 'float') + 1 - long.value / 360, format='jd')
    return local_midnight

def get_viewer_location(lat: u.quantity.Quantity, long: u.quantity.Quantity = 0*u.deg) -> EarthLocation:
    location = EarthLocation(lat=lat, lon=long, height = 0 * u.m)
    return location

def get_object_alt_az(object_name: str, time_range: np.ndarray, location: EarthLocation) -> SkyCoord:
    # local frame so coords can be converted to alt, az
    local_frame = AltAz(obstime=time_range, location=location)
    if object_name == 'sun' or object_name == 'moon':
        # uses jpl ephemerides 
        object_alt_az = get_body(object_name, time_range).transform_to(local_frame)
    else:
        # uses Simbad to resolve object names and retrieve coordinates
        object_alt_az = SkyCoord.from_name(name=object_name, frame=local_frame)
    return object_alt_az

def plot_object(object_name: str, date: Time, midnight_deltaT: np.ndarray, location: EarthLocation):

    local_midnight = get_local_midnight(date)

    night_hrs_vec = local_midnight + midnight_deltaT * u.hour
    object_alt_az = get_object_alt_az(object_name, night_hrs_vec, location)
    metrics = Observing_Metrics(object_alt_az, local_midnight, location)
    # determine Moon position
    moon_alt_az = get_object_alt_az('moon', night_hrs_vec, location)
    sun_alt_az = get_object_alt_az('sun', night_hrs_vec, location)
    
    year = str(date.ymdhms[0])
    month = str(date.ymdhms[1])
    day = str(date.ymdhms[2])

    astro_dark_deltaT = midnight_deltaT[sun_alt_az.alt.degree <= -18]
    # , subplot_kw={'projection': 'polar'}
    fig, ax = polar_subplots(figsize = (10, 10), subplot_kw={'projection': 'polar'})

    color = sun_alt_az.alt.degree
    color[sun_alt_az.alt.degree > 0] = 0
    color[sun_alt_az.alt.degree < -18] = -18

    obj_above_horizon = object_alt_az.alt.radian >= 0
    moon_above_horizon = moon_alt_az.alt.radian >= 0

    lines = colored_line_between_pts(object_alt_az.az.radian[obj_above_horizon], np.cos(object_alt_az.alt.radian[obj_above_horizon]), color[obj_above_horizon], ax, linewidth=5, cmap="cividis")
    p = ax.plot(moon_alt_az.az.radian[moon_above_horizon], np.cos(moon_alt_az.alt.radian[moon_above_horizon]), linewidth=3, color='k', label='Moon', alpha=0.8)
    ax.set_rlim([0, 1])
    ax.set_rticks(np.cos(np.radians([15, 30, 45, 60, 75])), labels=['15', '30', '45', '60', '75'])
    ax.set_theta_zero_location('N')
    
    ax.set_title(object_name + " on the night of " + month + "/" + day + "/" + year + " at " + str(location.lat.value) + " Lat")

    handles, labels = ax.get_legend_handles_labels()
    proxy = make_proxy(np.min(color[obj_above_horizon]), lines, linewidth=5)
    handles.insert(0, proxy)
    labels.insert(0, object_name)

    ax.legend(handles=handles, labels=labels, loc="upper left")

    fig.tight_layout()
    fig.savefig(object_name.replace(" ", "_") + "_altitude_" + month + "_" + day + "_" + year, dpi=300)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--name', required=True, help="Name of object")
    parser.add_argument('-l', '--latitude', required=True, type=float, help="Latitude as decimal degrees")
    # parser.add_argument('-m', '--meridian', required=True, type=float, help="Longitude (meridian) as decimal degrees")
    parser.add_argument('-d', '--date', required=True, help="Date of the observation (mm/dd/yyyy)")
    args = parser.parse_args()

    month, day, year = args.date.split("/")
    date = Time(year + "-" + month + "-" + day, format='iso', scale='utc')
    location = get_viewer_location(args.latitude * u.deg)
    local_midnight = get_local_midnight(date)

    # need to go a bit past 24 hrs since local_midnight isn't centered on the solar day
    n = 1000
    range_to_search = np.linspace(-12.5, 12.5, n)
    sun_alt_az = get_object_alt_az('sun', local_midnight + range_to_search * u.hour, location)

    midnight_deltaT = range_to_search[sun_alt_az.alt.degree <= 0]

    if (midnight_deltaT.size == 0):
        print("Sun does not set.")
        exit()

    astro_dark_deltaT = range_to_search[sun_alt_az.alt.degree <= 18]

    if (astro_dark_deltaT.size == 0):
        print("No astronomical dark.")
        exit()

    plot_object(args.name, date, midnight_deltaT, location)

if __name__ == "__main__":
    main()
