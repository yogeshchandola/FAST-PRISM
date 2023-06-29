#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 29 14:38:36 2019

@author: yogesh
"""
import numpy as np
#import scipy as sc
#import math as mh
#import ephem as em
from astropy.coordinates import SkyCoord, SpectralCoord,EarthLocation, FK5 #Angle
import astropy
#from astropy.coordinates import EarthLocation
from astropy.time import Time
from astropy import units as u
#from astroplan import download_IERS_A
from astropy.utils import iers
import process.modules.feedposition as fd
import math
from astropy.utils.iers import conf
  
conf.auto_max_age = None
######################### Co-ordinate calculation ##########################

beams0 = {'beam_name':['M01','M02','M03','M04','M05','M06','M07','M08','M09','M10','M11','M12','M13','M14',
                      'M15','M16','M17','M18','M19'],
         'beam_offset_ra':[0.,5.74,2.88,-2.86,-5.74,-2.88,2.86,
                          11.5,8.61,5.75,0.018,-5.71,-8.6,-11.5,
                          -8.63,-5.77,-0.0181,5.73,8.61],
         'beam_offset_dec':[0.,0.00811,-4.97,-4.98,-0.0127,4.97,4.98,
                           0.0116,-4.96,-9.93,-9.94,-9.96,-4.99,-0.03,
                           4.95,9.93,9.94,9.95,4.98]}

###################### position of different beams relative to central beam (Peng Jiang et al. 2020)

#beams0 = {'beam_name':['M01','M02','M03','M04','M05','M06','M07','M08','M09','M10','M11','M12','M13','M14',
#                      'M15','M16','M17','M18','M19'],
#         'beam_offset_ra':[-0.04,5.76,2.86,-2.89,-5.78,-2.92,2.85,
#                          11.55,8.65,5.78,0.018,-5.78,-8.67,-11.61,
#                          -8.68,-5.83,-0.03,5.76,8.63],
#         'beam_offset_dec':[-0.02,-0.01,-4.98,-5.01,-0.02,4.98,4.97,
#                           -0.02,-5.01,-10.02,-10.04,-10.07,-5.05,-0.02,
#                           4.99,10.02,10.00,10.00,4.98]}



#pi = 3.1415926535



class coordinate:
        def __init__(self,b,beamsize,rang,x,y,z,t): # constructor
            self.b=b
            self.beamsize=beamsize
            self.rang=rang
            self.x=x
            self.y=y
            self.z=z   
            self.t=t
    
            
        def cor(self): 
            # position of different beams relative to the central beam after rotation
            bra = np.array(beams0['beam_offset_ra'])
            bdec = np.array(beams0['beam_offset_dec'])
            bang = (np.arctan(bdec/bra)/np.pi) * 180.0
            ind2 = np.where(bra<0)[0]
            ind4 = np.where((bra>0) & (bdec<0))[0]
            bang[ind2] = bang[ind2] + 180.
            bang[ind4] = bang[ind4] + 360.
            bang[0]=0.
            bang_rot=bang+self.rang
            beamdist = np.sqrt(bra**2 + bdec**2)
            beams = {'beam_name':beams0['beam_name'], \
                          'beam_offset_ra':list(beamdist*np.cos(bang_rot/180.*np.pi)),\
                          'beam_offset_dec':list(beamdist*np.sin(bang_rot/180.*np.pi))}
            
    
            # method to estimate ra,dec  of central beam            
            location=EarthLocation.from_geocentric(self.x,self.y,self.z,unit=u.m)
            lon,lat,h=location.to_geodetic()   
            time = Time(self.t, format='mjd', location =location)
            f=fd.feedposition(self.t)
            el,az=f.altaz()
            iers.Conf.iers_auto_url.set('ftp://cddis.gsfc.nasa.gov/pub/products/iers/finals2000A.all')
            #iers.Conf.iers_auto_url.set('https://datacenter.iers.org/eop.php')
            time.delta_ut1_utc
            coor=SkyCoord(alt=el,az=az, unit=(u.deg, u.deg),frame='altaz',obstime=time,location=location) 
            coord=coor.transform_to(frame='fk5') 
            ra=coord.ra.hourangle
            dec=coord.dec.deg 
            dec=dec+beams['beam_offset_dec'][self.b]/60 
            ra=ra+(beams['beam_offset_ra'][self.b]/(15*60)/np.cos(dec*np.pi/180))
          #  r=2*self.beamsize
          #  if(self.b == 0):
          #      dR=0.0
         #       dD=0.0
         #   if(self.b > 0 & self.b <=6): 
         #       #dR=r*math.cos(np.pi+((self.b-1)*(np.pi/3))+self.rang*np.pi/180)
         #       #dD=r*math.sin(np.pi+((self.b-1)*(np.pi/3))+self.rang*np.pi/180)
         #       dR=r*math.cos(((1-self.b)*(np.pi/3))-self.rang*np.pi/180)
         #       dD=r*math.sin(((1-self.b)*(np.pi/3))-self.rang*np.pi/180)
         #   if(self.b > 6):
         #       if(self.b%2!=0):
         #        #   dR=2*r*math.cos(np.pi+((self.b-7)*(np.pi/6))+self.rang*np.pi/180)
         #       #   dD=2*r*math.sin(np.pi+((self.b-7)*(np.pi/6))+self.rang*np.pi/180) 
         #           dR=2*r*math.cos(((7-self.b)*(np.pi/6))-self.rang*np.pi/180)
         #           dD=2*r*math.sin(((7-self.b)*(np.pi/6))-self.rang*np.pi/180) 
         #       if(self.b%2==0):
         #        #   dR=2*r*math.cos(np.pi/6)*math.cos(np.pi+((self.b-7)*(np.pi/6))+self.rang*np.pi/180)
         #        #   dD=2*r*math.cos(np.pi/6)*math.sin(np.pi+((self.b-7)*(np.pi/6))+self.rang*np.pi/180)
         #           dR=2*r*math.cos(np.pi/6)*math.cos(((7-self.b)*(np.pi/6))-self.rang*np.pi/180)
         #           dD=2*r*math.cos(np.pi/6)*math.sin(((7-self.b)*(np.pi/6))-self.rang*np.pi/180)
            #print(dR,dD)
           # ra=ra+(dR/(15*60)/np.cos(dec*np.pi/180))
           # dec=dec+(dD/60)
            #coord1=SkyCoord(ra=ra,dec=dec,unit=(u.hourangle, u.deg),obstime=time)    
            #j2000c=FK5(equinox=Time(2000.0, format='jyear'))
            #coord2=coord1.transform_to(j2000c) 
            #ra=coord2.ra.hourangle
            #dec=coord2.dec.deg
           # sour=SpectralCoord(self.freq*u.MHz, observer=location.get_itrs(obstime=time),target=[ra,dec])
           # newsour=sour.with_observer_stationary_relative_to('hcrs')
           # print(newsour)
           # dopcor=self.freq[0]-new
            return el, az, ra, dec
        
        
        def dopplercorrection(self):
            location=EarthLocation.from_geocentric(self.x,self.y,self.z,unit=u.m)
            time = Time(self.t, format='mjd', location =location)
            #f=fd.feedposition(self.t)
            #el,az,rag=f.altaz()
            iers.Conf.iers_auto_url.set('ftp://cddis.gsfc.nasa.gov/pub/products/iers/finals2000A.all')
            #iers.Conf.iers_auto_url.set('https://datacenter.iers.org/eop.php')
            time.delta_ut1_utc
            
            return newsour
      
        def cordrift(self): # method to estimate ra dec for driftscan
            location=EarthLocation.from_geocentric(self.x,self.y,self.z,unit=u.m)
            lon,lat,h=location.to_geodetic()
            time = Time(self.t, format='mjd', location =location)
            #f=fd.feedposition(self.t)
            #el,az,rag=f.altaz()
            iers.Conf.iers_auto_url.set('ftp://cddis.gsfc.nasa.gov/pub/products/iers/finals2000A.all')
            #iers.Conf.iers_auto_url.set('https://datacenter.iers.org/eop.php')
            time.delta_ut1_utc
            lst_apparent = time.sidereal_time('apparent') 
            dec='+31:14:12.9'
            coor=SkyCoord(ra=lst_apparent,dec=dec, unit=(u.hourangle, u.deg),frame='fk5',obstime=time,location=location)
            j2000c=FK5(equinox=Time(2000.0, format='jyear'))
            coord=coor.transform_to(j2000c)    
            ra=coord.ra.hourangle
            dec=coord.dec.deg  
        
            r=2*self.beamsize
            if(self.b == 0):
                dR=0.0
                dD=0.0
            if(self.b > 0 & self.b <=6): 
                dR=r*math.cos(np.pi+((self.b-1)*(np.pi/3))+self.rang*np.pi/180)
                dD=r*math.sin(np.pi+((self.b-1)*(np.pi/3))+self.rang*np.pi/180)
            if(self.b > 6):
                if(self.b%2!=0):
                    dR=2*r*math.cos(np.pi+((self.b-7)*(np.pi/6))+self.rang*np.pi/180)
                    dD=2*r*math.sin(np.pi+((self.b-7)*(np.pi/6))+self.rang*np.pi/180)                 
                if(self.b%2==0):
                    dR=2*r*math.cos(np.pi/6)*math.cos(np.pi+((self.b-7)*(np.pi/6))+self.rang*np.pi/180)
                    dD=2*r*math.cos(np.pi/6)*math.sin(np.pi+((self.b-7)*(np.pi/6))+self.rang*np.pi/180)
                              
            #print(dR,dD)
            ra=ra+(dR/(15*60))
            dec=dec+(dD/60)
            return ra, dec            
            #h=np.sqrt(np.power(self.x,2)+np.power(self.y,2))
            #dec=mh.atan(self.z/h)*(180/np.pi)
           # coord= SkyCoord(0,lat,unit=(u.hourangle,u.deg),frame='fk5')
           # coord=coord.to_string('hmsdms')
            #dec=coord.dec.deg  #to_string('dms')
            #return dec
            #time = Time('2016-09-15T15:32:39.000Z', format='isot',scale='utc', location = location)
            #time.ut1              
            #iers.IERS_A_URL = 'http://toshi.nofs.navy.mil/ser7/finals2000A.all'
            #download_IERS_A()
            #IERS-A default file name, URL, and ReadMe with content description
            #iers.IERS_A_FILE = 'finals2000A.all'
            #iers.IERS_A_URL = 'ftp://cddis.gsfc.nasa.gov/pub/products/iers/finals2000A.all'
            #iers.IERS_A_URL_MIRROR = 'https://datacenter.iers.org/data/9/finals2000A.all'
            #iers.IERS_A_README ='data/ReadMe.finals2000A'
            #iers.Conf.iers_auto_url.set('https://cddis.gsfc.nasa.gov/pub/products/iers/finals2000A.all')
            #time.delta_ut1_utc
            #lst_apparent = time.sidereal_time('apparent')   

        

