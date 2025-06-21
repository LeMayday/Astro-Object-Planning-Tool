# Date: 6/29/2024
# Helper functions for object_altitude_tool.py

import numpy as np
from numpy.typing import NDArray
from astropy.coordinates import AltAz, EarthLocation, SkyCoord, get_body
from astropy.time import Time
import astropy.units as u


class Observing_Metrics:
    '''
    Class object that calculates and stores observing metrics for the night
    time_to_astro_dark (time between sunset and astro dark)
    time_visible (time above horizon and during astro dark)
    time_optimal (time above 15 deg and during astro dark)
    '''
    @staticmethod
    def __compute_deltaT_subject_to(deltaT: np.ndarray, condition: NDArray[np.bool_]) -> float:
        '''
        Computes time difference between beginning and end of array subject to constraint.
        Size of deltaT should be the same as size of condition
        '''
        deltaT = deltaT[condition]
        if deltaT.size == 0:
            return 0
        return deltaT[-1] - deltaT[0]

    @classmethod
    def __init__(self, midnight_deltaT: np.ndarray):
        self.__midnight_deltaT = midnight_deltaT

    @classmethod
    def __time_to_astro_dark(self, sun_alt_az: SkyCoord) -> float:
        '''
        Computes time between sunset and beginning of astro dark
        '''
        during_astro_dark = sun_alt_az.alt.degree <= -18
        astro_dark_deltaT = self.__midnight_deltaT[during_astro_dark]
        return astro_dark_deltaT[0] - self.__midnight_deltaT[0]

    @classmethod
    def __time_above_horizon_and_astro_dark(self, object_alt_az: SkyCoord, sun_alt_az: SkyCoord) -> float:
        '''
        Compute the amount of time the object is above the horizon during astro dark
        '''
        above_horizon = object_alt_az.alt.degree >= 0
        during_astro_dark = sun_alt_az.alt.degree <= -18
        return self.__compute_deltaT_subject_to(self.__midnight_deltaT, np.logical_and(above_horizon, during_astro_dark))
    
    @classmethod
    def __time_above_15deg_and_astro_dark(self, object_alt_az: SkyCoord, sun_alt_az: SkyCoord) -> float:
        '''
        Computes the amount of time the object is above 15 degrees during astro dark
        '''
        above_15deg = object_alt_az.alt.degree >= 15
        during_astro_dark = sun_alt_az.alt.degree <= -18
        return self.__compute_deltaT_subject_to(self.__midnight_deltaT, np.logical_and(above_15deg, during_astro_dark))

    @classmethod
    def compute_metrics(self, object_alt_az: SkyCoord, sun_alt_az: SkyCoord):
        '''
        Computes the desired quality metrics for the observing session
        '''
        self.time_to_astro_dark = self.__time_to_astro_dark(sun_alt_az)
        self.time_visible = self.__time_above_horizon_and_astro_dark(object_alt_az, sun_alt_az)
        self.time_optimal = self.__time_above_15deg_and_astro_dark(object_alt_az, sun_alt_az)


def get_local_midnight(date: Time, long: u.quantity.Quantity = 0*u.deg) -> Time:
    '''
    Gets the local midnight based on date of observation and longitude
    '''
    # +1 for night of date, - long / 360 to convert to local midnight
    local_midnight = Time(date.to_value('jd', 'float') + 1 - long.value / 360, format='jd')
    return local_midnight

def get_viewer_location(lat: u.quantity.Quantity, long: u.quantity.Quantity = 0*u.deg) -> EarthLocation:
    '''
    Gets an object representing the location of the observer
    '''
    location = EarthLocation(lat=lat, lon=long, height = 0 * u.m)
    return location

def get_object_alt_az(object_name: str, time_range: np.ndarray, location: EarthLocation) -> SkyCoord:
    '''
    Gets an object storing altitude and azimuth of the object as viewed from the specified location during the specified time
    Requires an internet connection
    '''
    # local frame so coords can be converted to alt, az
    local_frame = AltAz(obstime=time_range, location=location)
    if object_name == 'sun' or object_name == 'moon':
        # uses jpl ephemerides 
        object_alt_az = get_body(object_name, time_range).transform_to(local_frame)
    else:
        # uses Simbad to resolve object names and retrieve coordinates
        object_alt_az = SkyCoord.from_name(name=object_name, frame=local_frame)
    return object_alt_az

def moon_rising(moon_alt_az: SkyCoord) -> bool:
    '''
    Determines whether the Moon is primarily rising or setting given the alt_az coordinates on the night of observation
    '''
    moon_above_horizon = moon_alt_az.alt.radian >= 0
    moon_az_above_horizon = moon_alt_az.az.degree[moon_above_horizon]
    avg_moon_az = (moon_az_above_horizon[-1] + moon_az_above_horizon[0])/2
    return avg_moon_az <= 180
