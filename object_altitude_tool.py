# Date: 6/29/2024
# Edited: 6/20/2025
# inspired by https://docs.astropy.org/en/stable/generated/examples/coordinates/plot_obs-planning.html
# requires matplotlib, numpy, astropy, and jplephem modules

import argparse
from typing import Tuple
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 16})
from matplotlib.figure import Figure
from matplotlib.projections.polar import PolarAxes
import numpy as np
import astropy.units as u
from astropy.coordinates import EarthLocation
from astropy.time import Time
from plot_utils import colored_line_between_pts, polar_subplots, make_proxy, project_onto_polar, PROJECTIONS
from observing_utils import Observing_Metrics, get_local_midnight, get_viewer_location, get_object_alt_az, moon_rising


def plot_object(object_name: str, night_hrs_vec: Time, location: EarthLocation, projection: str) -> Tuple[Figure, PolarAxes]:
    '''
    Plot object and Moon during given night hours
    '''
    object_alt_az = get_object_alt_az(object_name, night_hrs_vec, location)
    moon_alt_az = get_object_alt_az('moon', night_hrs_vec, location)
    sun_alt_az = get_object_alt_az('sun', night_hrs_vec, location)
    Observing_Metrics.compute_metrics(object_alt_az, sun_alt_az)

    fig, ax = polar_subplots(figsize = (10, 10), subplot_kw={'projection': 'polar'})
    # maps colors based on sun altitude (lighest when sun is above horizon, darkest during astro dark)
    color = sun_alt_az.alt.degree
    color[sun_alt_az.alt.degree > 0] = 0
    color[sun_alt_az.alt.degree < -18] = -18
    # define boolean arrays to make indexing more intuitive
    obj_above_horizon = object_alt_az.alt.radian >= 0
    moon_above_horizon = moon_alt_az.alt.radian >= 0
    # plots object alt vs az with sun alt used for color
    lines = colored_line_between_pts(object_alt_az.az.radian[obj_above_horizon], project_onto_polar(object_alt_az.alt.radian[obj_above_horizon], projection),
                                     color[obj_above_horizon][:-1], ax, linewidth=5, cmap="cividis")
    # plots moon alt vs az
    ax.plot(moon_alt_az.az.radian[moon_above_horizon], project_onto_polar(moon_alt_az.alt.radian[moon_above_horizon], projection), 
                linewidth=3, color='k', label='Moon', alpha=0.8)
    ax.set_rlim([0, 1])
    # draws lines of constant altitude
    ax.set_rticks(project_onto_polar(np.radians([0, 15, 30, 45, 60, 75, 90]), projection), 
                  labels=['', r'15$^\degree$', r'30$^\degree$', r'45$^\degree$', r'60$^\degree$', r'75$^\degree$', ''])
    ax.set_theta_zero_location('N')
    
    # adds ticks showing time since sunset to axis
    for i in range(1, int((night_hrs_vec[-1] - night_hrs_vec[0]).value * 24) + 1):
        tick_time = night_hrs_vec[0] + i * u.hour
        tick_ojb_alt_az = get_object_alt_az(object_name, tick_time, location)
        if tick_ojb_alt_az.alt.radian >= 0:
            ax.text(tick_ojb_alt_az.az.radian, project_onto_polar(tick_ojb_alt_az.alt.radian, projection), f"+{i}", 
                    horizontalalignment='left', verticalalignment='bottom')
        tick_moon_alt_az = get_object_alt_az('moon', tick_time, location)
        if tick_moon_alt_az.alt.radian >= 0:
            ax.text(tick_moon_alt_az.az.radian, project_onto_polar(tick_moon_alt_az.alt.radian, projection), f"+{i}", 
                    horizontalalignment='left', verticalalignment='bottom')
    ax.text(np.radians(180), 1.1, "**Note: Ticks along object paths show hours after sunset.", multialignment='right', horizontalalignment='center', verticalalignment='top')

    # report observing metrics
    text_str = f"""Time Visible (during Astro Dark): {Observing_Metrics.time_visible:.1f} hrs
                Time Above 15 deg (during Astro Dark): {Observing_Metrics.time_optimal:.1f} hrs
                Time from Sunset to Astro Dark: {Observing_Metrics.time_to_astro_dark:.1f} hrs"""
    ax.text(np.radians(180), 1.1, text_str, multialignment='right', horizontalalignment='center', verticalalignment='top')
    # configure legend
    handles, labels = ax.get_legend_handles_labels()
    # LineCollection cannot be added to the legend except through proxy artist
    proxy = make_proxy(np.mean(color[obj_above_horizon]), lines, linewidth=5)
    handles.insert(0, proxy)
    labels.insert(0, object_name)
    # changes legend location depending on whether the Moon is primarily rising or setting
    legend_loc = "upper left" if moon_rising(moon_alt_az) else "upper right"
    ax.legend(handles=handles[1:], labels=labels[1:], loc=legend_loc, edgecolor='None', facecolor='#e6e6e6')
    fig.set_facecolor('#e6e6e6')
    return fig, ax

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--name', required=True, type= str, help="Name of object")
    parser.add_argument('-l', '--latitude', required=True, type=float, help="Latitude as decimal degrees")
    # parser.add_argument('-m', '--meridian', required=True, type=float, help="Longitude (meridian) as decimal degrees")
    parser.add_argument('-d', '--date', required=True, type=str, help="Date of the observation (mm/dd/yyyy)")
    parser.add_argument('-p', '--projection', required=False, type=str, help="Projection onto polar", choices=PROJECTIONS, default='linear')
    args = parser.parse_args()
    # define info about the observing session from arguments
    month, day, year = args.date.split("/")
    date = Time(year + "-" + month + "-" + day, format='iso', scale='utc')
    location = get_viewer_location(args.latitude * u.deg)
    local_midnight = get_local_midnight(date)

    # need to go a bit past 24 hrs since local_midnight isn't centered on the solar day
    n_search = 1000
    range_to_search = np.linspace(-12.5, 12.5, n_search)
    sun_alt_az = get_object_alt_az('sun', local_midnight + range_to_search * u.hour, location)
    # time range when sun is below horizon in hrs after local midnight
    midnight_deltaT = range_to_search[sun_alt_az.alt.degree <= 0]
    # exit if the sun never sets
    if (midnight_deltaT.size == 0):
        print("Sun does not set.")
        exit()
    # exit if the sun never goes below 18 degrees
    # this isn't really a deal breaker, but it probably won't be dark enough to get good images
    astro_dark_deltaT = range_to_search[sun_alt_az.alt.degree <= 18]
    if (astro_dark_deltaT.size == 0):
        print("No astronomical dark.")
        exit()

    # redefine midnight_deltaT to have more granularity
    n = 5000
    midnight_deltaT = np.linspace(midnight_deltaT[0], midnight_deltaT[-1], n)
    night_hrs_vec = local_midnight + midnight_deltaT * u.hour
    Observing_Metrics(midnight_deltaT)
    # plot object, set title, and save file
    fig, ax = plot_object(args.name, night_hrs_vec, location, args.projection)
    ax.set_title(args.name + " on the night of " + month + "/" + day + "/" + year + " at " + str(args.latitude) + " Lat")
    fig.tight_layout()
    fig.savefig(args.name.replace(" ", "-") + "_altaz_" + month + "-" + day + "-" + year, dpi=300)

if __name__ == "__main__":
    main()
