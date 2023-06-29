#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 29 14:42:52 2019

@author: yogesh
"""

from skyfield.api import load, Topos
from astropy.time import Time
from skyfield.timelib import Time as SkyfieldTime
import numpy as np
# Constants
c = 299792.458  # Speed of light in km/s

def heliocentric_velocity_correction(galaxy_velocity, observer_latitude, observer_longitude, observation_time):
    """
    Corrects the velocity of a distant galaxy from the Earth's frame of rest to the Sun's frame of rest,
    accounting for the Doppler shift due to the Earth's rotation, Earth's orbital velocity, and the galaxy's velocity in the Earth's frame,
    considering the precise time of observation.
    
    Args:
        galaxy_velocity (float): Velocity of the distant galaxy in km/s in the Earth's frame of rest.
        observer_latitude (float): Latitude of the observer on Earth in degrees.
        observer_longitude (float): Longitude of the observer on Earth in degrees.
        observation_time (float): Time of observation in Julian Date (JD) or any other suitable time format.
    
    Returns:
        float: Corrected velocity of the distant galaxy in km/s in the Sun's frame of rest.
    """
    # Load the ephemeris data for the Earth and Sun
    filepath='/sugon/MaGroup/chandola/fastpipeline/process/'
    planets = load(filepath+'de421.bsp')
    earth = planets['earth']
    sun = planets['sun']
    ts=load.timescale()
    tjd=ts.ut1_jd(observation_time.jd)
    #t=Time(observation_time,format='mjd',scale='utc')
    #print(tjd)
   # observation_datetime=t.utc.datetime
    
   # t_utc=SkyfieldTime(datetime=observation_datetime,scale='utc')
    
    # Get the Earth's position at the observation time
    observer_position = earth + Topos(latitude_degrees=observer_latitude, longitude_degrees=observer_longitude)
    observer_ephem = observer_position.at(tjd)
    observer_xyz = observer_ephem.position.km
    
    # Get the Sun's position at the observation time
    sun_ephem = sun.at(tjd)
    sun_xyz = sun_ephem.position.km
    
    # Calculate the velocity due to the Earth's rotation
    observer_velocity = observer_ephem.velocity.km_per_s
    rotation_velocity = np.dot(observer_velocity, observer_xyz) / np.linalg.norm(observer_xyz)
    
    # Calculate the velocity due to the Earth's orbital motion
    orbital_velocity = np.dot(observer_velocity, sun_xyz) / np.linalg.norm(observer_xyz)
    print('Julian day:',tjd,'Doppler shift in velocity:',rotation_velocity+orbital_velocity)
    # Calculate the velocity of the distant galaxy in the Sun's frame of rest
    heliocentric_velocity = galaxy_velocity + rotation_velocity + orbital_velocity
    
    return heliocentric_velocity