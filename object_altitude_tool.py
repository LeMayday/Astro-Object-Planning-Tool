# Date: 6/29/2024
# Edited: 6/20/2025
# inspired by https://docs.astropy.org/en/stable/generated/examples/coordinates/plot_obs-planning.html
# requires matplotlib, numpy, astropy, and jplephem modules

import argparse
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 16})
import numpy as np
import astropy.units as u
from astropy.coordinates import EarthLocation
from astropy.time import Time
from plot_utils import colored_line_between_pts, polar_subplots, make_proxy, project_onto_polar, PROJECTIONS
from observing_utils import Observing_Metrics, get_local_midnight, get_viewer_location, get_object_alt_az

def plot_object(object_name: str, night_hrs_vec: Time, date: Time, location: EarthLocation, projection: str):

    object_alt_az = get_object_alt_az(object_name, night_hrs_vec, location)
    # determine Moon position
    moon_alt_az = get_object_alt_az('moon', night_hrs_vec, location)
    sun_alt_az = get_object_alt_az('sun', night_hrs_vec, location)
    Observing_Metrics.compute_metrics(object_alt_az, sun_alt_az)

    fig, ax = polar_subplots(figsize = (10, 10), subplot_kw={'projection': 'polar'})

    color = sun_alt_az.alt.degree
    color[sun_alt_az.alt.degree > 0] = 0
    color[sun_alt_az.alt.degree < -18] = -18

    obj_above_horizon = object_alt_az.alt.radian >= 0
    moon_above_horizon = moon_alt_az.alt.radian >= 0

    lines = colored_line_between_pts(object_alt_az.az.radian[obj_above_horizon], project_onto_polar(object_alt_az.alt.radian[obj_above_horizon], projection),
                                     color[obj_above_horizon][:-1], ax, linewidth=5, cmap="cividis")
    p = ax.plot(moon_alt_az.az.radian[moon_above_horizon], project_onto_polar(moon_alt_az.alt.radian[moon_above_horizon], projection), 
                linewidth=3, color='k', label='Moon', alpha=0.8)
    ax.set_rlim([0, 1])
    ax.set_rticks(project_onto_polar(np.radians([0, 15, 30, 45, 60, 75, 90]), projection), labels=['', '15', '30', '45', '60', '75', ''])
    ax.set_theta_zero_location('N')

    text_str = f"Time Visible (Dark): {Observing_Metrics.time_visible:.2f} hrs\nTime Above 15 deg (Dark): {Observing_Metrics.time_optimal:.2f} hrs\nTime to Astro Dark: {Observing_Metrics.time_to_astro_dark:.2f} hrs"
    ax.text(np.radians(180), 1.1, text_str, multialignment='right', horizontalalignment='center', verticalalignment='top')
    
    handles, labels = ax.get_legend_handles_labels()
    proxy = make_proxy(np.mean(color[obj_above_horizon]), lines, linewidth=5)
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
    parser.add_argument('-p', '--projection', required=False, type=str, help="Projection onto polar", choices=PROJECTIONS, default='linear')
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
    plot_object(args.name, night_hrs_vec, date, location, args.projection)

if __name__ == "__main__":
    main()
