# Author: Mayday
# Date: 6/29/2024
# inspired by https://docs.astropy.org/en/stable/generated/examples/coordinates/plot_obs-planning.html
# requires matplotlib, numpy, astropy, and jplephem modules

from typing import Tuple
import argparse
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 16})
import numpy as np
from numpy.typing import NDArray
import astropy.units as u
from astropy.coordinates import AltAz, EarthLocation, SkyCoord, get_body
from astropy.time import Time
import plot_utils

class Observing_Metrics:
    '''
    Class object that calculates and stores observing metrics for the night
    time_above_horizon
    time_above_15_deg
    time_astro_dark (also above horizon)
    time_optimal (time above 15 deg and during astro dark)
    '''
    @staticmethod
    def __compute_deltaT_subject_to(deltaT: np.ndarray, condition: NDArray[np.bool_]):
        deltaT = deltaT[condition]
        if deltaT.size == 0:
            return 0
        return deltaT[-1] - deltaT[0]

    @classmethod
    def __init__(self, midnight_deltaT: np.ndarray):
        self.__midnight_deltaT = midnight_deltaT

    @classmethod
    def __time_to_astro_dark(self, sun_alt_az: SkyCoord) -> float:
        during_astro_dark = sun_alt_az.alt.degree <= -18
        astro_dark_deltaT = self.__midnight_deltaT[during_astro_dark]
        return astro_dark_deltaT[0] - self.__midnight_deltaT[0]

    @classmethod
    def __time_above_horizon_and_astro_dark(self, object_alt_az: SkyCoord, sun_alt_az: SkyCoord) -> float:
        above_horizon = object_alt_az.alt.degree >= 0
        during_astro_dark = sun_alt_az.alt.degree <= -18
        return self.__compute_deltaT_subject_to(self.__midnight_deltaT, above_horizon and during_astro_dark)
    
    @classmethod
    def __time_above_15deg_and_astro_dark(self, object_alt_az: SkyCoord, sun_alt_az: SkyCoord) -> float:
        above_15deg = object_alt_az.alt.degree >= 15
        during_astro_dark = sun_alt_az.alt.degree <= -18
        return self.__compute_deltaT_subject_to(self.__midnight_deltaT, above_15deg and during_astro_dark)

    @classmethod
    def compute_metrics(self, object_alt_az: SkyCoord, sun_alt_az: SkyCoord):
        self.time_to_astro_dark = self.__time_to_astro_dark(sun_alt_az)
        self.time_visible = self.__time_above_horizon_and_astro_dark(object_alt_az, sun_alt_az)
        self.time_optimal = self.__time_above_15deg_and_astro_dark(object_alt_az, sun_alt_az)


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

def plot_object(object_name: str, night_hrs_vec: np.ndarray, date: Time, location: EarthLocation, projection: str):

    object_alt_az = get_object_alt_az(object_name, night_hrs_vec, location)
    # determine Moon position
    moon_alt_az = get_object_alt_az('moon', night_hrs_vec, location)
    sun_alt_az = get_object_alt_az('sun', night_hrs_vec, location)

    fig, ax = plot_utils.polar_subplots(figsize = (10, 10), subplot_kw={'projection': 'polar'})

    color = sun_alt_az.alt.degree
    color[sun_alt_az.alt.degree > 0] = 0
    color[sun_alt_az.alt.degree < -18] = -18

    obj_above_horizon = object_alt_az.alt.radian >= 0
    moon_above_horizon = moon_alt_az.alt.radian >= 0

    lines = plot_utils.colored_line_between_pts(object_alt_az.az.radian[obj_above_horizon], plot_utils.project_onto_polar(object_alt_az.alt.radian[obj_above_horizon], projection),
                                                color[obj_above_horizon], ax, linewidth=5, cmap="cividis")
    p = ax.plot(moon_alt_az.az.radian[moon_above_horizon], plot_utils.project_onto_polar(moon_alt_az.alt.radian[moon_above_horizon], projection), 
                linewidth=3, color='k', label='Moon', alpha=0.8)
    ax.set_rlim([0, 1])
    ax.set_rticks(plot_utils.project_onto_polar(np.radians([0, 15, 30, 45, 60, 75, 90]), projection), labels=['', '15', '30', '45', '60', '75', ''])
    ax.set_theta_zero_location('N')
    
    handles, labels = ax.get_legend_handles_labels()
    proxy = plot_utils.make_proxy(np.mean(color[obj_above_horizon]), lines, linewidth=5)
    handles.insert(0, proxy)
    labels.insert(0, object_name)
    ax.legend(handles=handles, labels=labels, loc="upper left")

    year = str(date.ymdhms[0])
    month = str(date.ymdhms[1])
    day = str(date.ymdhms[2])
    ax.set_title(object_name + " on the night of " + month + "/" + day + "/" + year + " at " + str(location.lat.value) + " Lat")
    fig.tight_layout()
    fig.savefig(object_name.replace(" ", "_") + "_altitude_" + month + "_" + day + "_" + year, dpi=300)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--name', required=True, type= str, help="Name of object")
    parser.add_argument('-l', '--latitude', required=True, type=float, help="Latitude as decimal degrees")
    # parser.add_argument('-m', '--meridian', required=True, type=float, help="Longitude (meridian) as decimal degrees")
    parser.add_argument('-d', '--date', required=True, type=str, help="Date of the observation (mm/dd/yyyy)")
    parser.add_argument('-p', '--projection', required=False, type=str, help="Projection onto polar", choices=plot_utils.projections, default='linear')
    args = parser.parse_args()

    month, day, year = args.date.split("/")
    date = Time(year + "-" + month + "-" + day, format='iso', scale='utc')
    location = get_viewer_location(args.latitude * u.deg)
    local_midnight = get_local_midnight(date)

    # need to go a bit past 24 hrs since local_midnight isn't centered on the solar day
    n_search = 1000
    range_to_search = np.linspace(-12.5, 12.5, n_search)
    sun_alt_az = get_object_alt_az('sun', local_midnight + range_to_search * u.hour, location)

    midnight_deltaT = range_to_search[sun_alt_az.alt.degree <= 0]

    if (midnight_deltaT.size == 0):
        print("Sun does not set.")
        exit()

    astro_dark_deltaT = range_to_search[sun_alt_az.alt.degree <= 18]

    if (astro_dark_deltaT.size == 0):
        print("No astronomical dark.")
        exit()

    n = 5000
    midnight_deltaT = np.linspace(midnight_deltaT[0], midnight_deltaT[-1], n)
    local_midnight = get_local_midnight(date)
    night_hrs_vec = local_midnight + midnight_deltaT * u.hour
    Observing_Metrics(midnight_deltaT)
    plot_object(args.name, date, night_hrs_vec, location, args.projection)

if __name__ == "__main__":
    main()
