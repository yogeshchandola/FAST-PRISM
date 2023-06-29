#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  2 18:45:12 2023

@author: chandola
"""

from astropy.time import Time
from astropy import units as u
#from astroplan import download_IERS_A
from astropy.utils import iers
from astropy.coordinates import SkyCoord, SpectralCoord,EarthLocation, FK5 

c=2.99792458*100000

def obsloc(x,y,z,t):
       location=EarthLocation.from_geocentric(x,y,z,unit=u.m)
       lon,lat,h=location.to_geodetic()   
       time = Time(t, format='mjd')
       iers.Conf.iers_auto_url.set('ftp://cddis.gsfc.nasa.gov/pub/products/iers/finals2000A.all')
            #iers.Conf.iers_auto_url.set('https://datacenter.iers.org/eop.php')
       time.delta_ut1_utc
       return lon, lat, time

def freqtovel(freq,mode,restfreq):
    if(mode=='radio'):
        velocity=((restfreq-freq)/restfreq)*c
    elif(mode=='optical'):
        velocity=((restfreq-freq)/freq)*c
    return velocity

