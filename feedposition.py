#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  2 11:22:14 2020

@author: yogesh
"""

import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
from astropy.coordinates import EarthLocation, SkyCoord, AltAz
from astropy import coordinates as coord
from astropy.time import Time
from astropy import units as u

#from astropy.utils import iers
#eloc= EarthLocation.of_site('Apache Point Observatory')  
#<EarthLocation (-1463969.30185172, -5166673.34223433,  3434985.71204565) m>

#lon,lat,h=location.to_geodetic()

convert=np.pi/180

#filepath='/sugon/MaGroup/chandola/PT2020_0079/KY/'
#feed=pd.read_excel(filepath+"R1d1_2020_09_14_13_30_15_000.xlsx",engine='openpyxl')
#feed=pd.read_excel(filepath+"R1d2_2020_09_11_13_42_02_000.xlsx",engine='openpyxl')
#feed=pd.read_excel(filepath+"R1d2_2020_09_18_13_16_27_000.xlsx",engine='openpyxl')
#feed=pd.read_excel(filepath+"R1d2_2020_10_02_12_19_28_000.xlsx",engine='openpyxl')
#feed=pd.read_excel(filepath+"R1d3_2020_09_15_13_26_19_000.xlsx",engine='openpyxl')
#feed=pd.read_excel(filepath+"R1d4_2020_09_16_13_26_22_000.xlsx",engine='openpyxl')
#feed=pd.read_excel(filepath+"R2d1_2020_09_14_14_28_01_000.xlsx",engine='openpyxl')
#feed2=pd.read_excel(filepath+"R2d1_2020_09_27_13_36_54_000.xlsx",engine='openpyxl')
#feed3=pd.read_excel(filepath+"R2d1_2020_10_06_13_01_31_000.xlsx",engine='openpyxl')
#feed=pd.read_excel(filepath+"R2d2_2020_09_25_13_44_46_000.xlsx",engine='openpyxl')
#feed=pd.read_excel(filepath+"R2d2_2020_09_29_13_29_02_000.xlsx",engine='openpyxl')
#feed=pd.read_excel(filepath+"R2d2_2020_10_07_12_57_35_000.xlsx",engine='openpyxl')
#feed=pd.read_excel(filepath+"R2d3_2020_09_24_13_48_42_000.xlsx",engine='openpyxl')
#feed=pd.read_excel(filepath+"R2d3_2020_10_02_13_17_15_000.xlsx",engine='openpyxl')
#feed=pd.read_excel(filepath+"R2d3_2020_10_09_12_49_43_000.xlsx",engine='openpyxl')

#filepath='/sugon/MaGroup/chandola/FAST_DATA/KY/PT2021_0034/'
#feed=pd.read_excel(filepath+"J022246.94-093848.7_2022_01_29_18_35_00_000.xlsx")
#feed=pd.read_excel(filepath+"J081437.98+172208.3_2021_10_02_06_20_00_000.xlsx")
#feed=pd.read_excel(filepath+"J081755.21+312827.4_2021_10_05_06_44_00_000.xlsx")
#feed=pd.read_excel(filepath+"J083216.03+183212.1_2021_10_20_06_20_00_000.xlsx")
#feed=pd.read_excel(filepath+"J083548.14+151717.0_2021_10_20_06_55_00_000.xlsx")
#feed=pd.read_excel(filepath+"J090410.36+024744.8_2021_10_23_06_53_00_000.xlsx")
#feed=pd.read_excel(filepath+"J090410.36+024744.8_2022_01_29_23_23_00_000.xlsx")
#feed=pd.read_excel(filepath+"J092527.55+072641.6_2021_10_22_06_34_00_000.xlsx")
#feed=pd.read_excel(filepath+"J092924.92+193421.0_2021_10_18_06_36_00_000.xlsx")
#feed=pd.read_excel(filepath+"J093242.81-003948.8-2_2022_01_30_01_25_00_000.xlsx")
#feed=pd.read_excel(filepath+"J094310.82+295203.6_2021_10_23_07_07_00_000.xlsx")
#feed=pd.read_excel(filepath+"J095058.69+375758.8_2021_11_10_05_53_00_000.xlsx")
#feed=pd.read_excel(filepath+"J102453.63+233234.0_2021_10_19_06_50_00_000.xlsx")
#feed=pd.read_excel(filepath+"J110701.20+182548.8_2021_11_11_06_35_00_000.xlsx")
#feed=pd.read_excel(filepath+"J114538.51+442021.9_2021_11_18_06_40_00_000.xlsx")
#feed=pd.read_excel(filepath+"J115712.38-032107.7_2022_01_30_04_46_18_000.xlsx")
#feed=pd.read_excel(filepath+"J121755.30-033723.3_2022_02_03_04_21_00_000.xlsx")
#feed=pd.read_excel(filepath+"J122113.25-024859.5_2021_12_02_06_53_00_000.xlsx")
#feed=pd.read_excel(filepath+"J122228.47+171437.3_2022_01_30_02_10_00_000.xlsx")
#feed=pd.read_excel(filepath+"J124419.96+405136.8_2021_11_22_06_53_30_000.xlsx")
#feed=pd.read_excel(filepath+"J124707.32+490017.9_2021_11_23_07_00_00_000.xlsx")
#feed=pd.read_excel(filepath+"J132522.00+035848.9_2021_12_19_07_45_00_000.xlsx")
#feed=pd.read_excel(filepath+"J132859.25+173842.3_2021_12_17_07_13_00_000.xlsx")
#feed=pd.read_excel(filepath+"J133242.53+134253.8_2022_02_02_05_31_00_000.xlsx")
#feed=pd.read_excel(filepath+"J135223.46-015648.4_2021_12_17_07_47_00_000.xlsx")
#feed=pd.read_excel(filepath+"J141327.22+550529.2_2021_12_20_06_55_00_000.xlsx")
#feed=pd.read_excel(filepath+"J143806.13+190954.9_2022_01_31_04_13_00_000.xlsx")
#feed=pd.read_excel(filepath+"J144920.71+422101.2_2021_12_20_07_31_00_000.xlsx")
#feed=pd.read_excel(filepath+"J145844.79+372021.5_2022_01_03_07_20_00_000.xlsx")
#feed=pd.read_excel(filepath+"J152142.58+181438.2_2022_01_03_07_56_00_000.xlsx")
#feed=pd.read_excel(filepath+"J153016.25+375831.2_2021_12_20_07_50_00_000.xlsx")
#feed=pd.read_excel(filepath+"J153229.40+015133.7_2022_01_04_07_47_00_000.xlsx")
#feed=pd.read_excel(filepath+"J153836.11+552541.4_2022_01_07_07_36_00_000.xlsx")
#feed=pd.read_excel(filepath+"J154345.80+110935.9_2022_01_04_07_14_00_000.xlsx")
#feed=pd.read_excel(filepath+"J155903.43+230828.7_2022_02_01_06_12_30_000.xlsx")
#feed=pd.read_excel(filepath+"J155927.67+533054.4_2022_01_05_07_35_00_000.xlsx")
#feed=pd.read_excel(filepath+"J162033.43+173955.5_2022_01_05_07_53_30_000.xlsx")
#feed=pd.read_excel(filepath+"J213333.31-071249.2_2021_09_20_22_53_00_000.xlsx")
#feed=pd.read_excel(filepath+"J230551.18-104052.2_2021_09_27_23_55_00_000.xlsx")
#feed=pd.read_excel(filepath+"J235400.91-003449.5_2021_09_28_00_28_00_000.xlsx")
#feed=pd.read_excel(filepath+"J233515.92-011216.8_2022_06_01_06_02_00_000.xlsx")

filepath='/sugon/MaGroup/chandola/FAST_DATA/KY/PT2022_0192/'
#feed=pd.read_excel(filepath+"BB10_2022_08_24_19_32_00_000.xlsx")
#feed=pd.read_excel(filepath+"BB106_2023_01_11_07_00_00_000.xlsx")
#feed=pd.read_excel(filepath+"BB12_2022_08_24_19_50_30_000.xlsx")
#feed=pd.read_excel(filepath+"BB13_2023_01_04_07_24_30_000.xlsx")
#feed=pd.read_excel(filepath+"BB27_2022_08_22_04_38_11_000.xlsx")
#feed=pd.read_excel(filepath+"BB28_2022_08_22_05_02_30_000.xlsx")
#feed=pd.read_excel(filepath+"BB31_2022_08_21_06_38_00_000.xlsx")
#feed=pd.read_excel(filepath+"BB3_2022_10_07_05_45_00_000.xlsx")
#feed=pd.read_excel(filepath+"BB34_2022_08_24_06_17_00_000.xlsx")
#feed=pd.read_excel(filepath+"BB39_2022_10_07_06_18_00_000.xlsx")
#feed=pd.read_excel(filepath+"BB41_2022_10_08_06_35_00_000.xlsx")
#feed=pd.read_excel(filepath+"BB4_2022_10_08_07_08_00_000.xlsx")
#feed=pd.read_excel(filepath+"BB47_2022_10_26_07_08_00_000.xlsx")
#feed=pd.read_excel(filepath+"BB5_2022_10_08_07_55_00_000.xlsx")
#feed=pd.read_excel(filepath+"BB57_2022_10_23_07_15_00_000.xlsx")
#feed=pd.read_excel(filepath+"BB6_2022_12_18_06_57_00_000.xlsx")
#feed=pd.read_excel(filepath+"BB66_2022_10_23_07_52_00_000.xlsx")
#feed=pd.read_excel(filepath+"BB67_2022_10_23_07_33_30_000.xlsx")
feed=pd.read_excel(filepath+"BB68_2022_11_07_06_57_00_000.xlsx")
#feed=pd.read_excel(filepath+"BB69_2022_11_04_06_41_00_000.xlsx")
#feed=pd.read_excel(filepath+"BB70_2022_11_17_07_26_30_000.xlsx")
#feed=pd.read_excel(filepath+"BB7_2022_12_18_07_16_00_000.xlsx")
#feed=pd.read_excel(filepath+"BB88_2022_12_13_06_50_00_000.xlsx")
#feed=pd.read_excel(filepath+"BB9_2022_12_24_07_19_30_000.xlsx")
#feed=pd.read_excel(filepath+"BB2_2023_04_12_19_30_00_000.xlsx")
#feed=pd.read_excel(filepath+"BB90_2023_04_12_00_47_17_000.xlsx")
#feed=pd.read_excel(filepath+"BB122_2023_05_18_06_00_00_000.xlsx")





#filepath='/sugon/MaGroup/chandola/FAST_DATA/KY/PT2021_0137/'
#feed=pd.read_excel(filepath+"0713+438_2021_10_14_06_13_30_000.xlsx")
#feed=pd.read_excel(filepath+"1120+143_2021_11_11_07_09_00_000.xlsx")
#feed=pd.read_excel(filepath+"1219+484_2021_11_13_07_59_00_000.xlsx")
#feed=pd.read_excel(filepath+"254+116_2021_11_29_06_50_00_000.xlsx")
#feed=pd.read_excel(filepath+"2139+143_2021_10_01_00_26_30_000.xlsx")
#feed=pd.read_excel(filepath+"J082018.44+562134.1_2021_10_14_05_57_00_000.xlsx")
#feed=pd.read_excel(filepath+"J102310.56+475145.6_2021_11_13_07_42_00_000.xlsx")
#feed=pd.read_excel(filepath+"J105817.90+195150.9_2021_11_11_06_54_00_000.xlsx")
#feed=pd.read_excel(filepath+"J123932.75+044305.2_2021_11_29_07_05_00_000.xlsx")
#feed=pd.read_excel(filepath+"J222646.53+005211.3_2021_10_01_00_08_00_000.xlsx")



#filepath='/sugon/MaGroup/chandola/FAST_DATA/PT2020_0082/KY/'
#feed=pd.read_excel(filepath+"BB3_2021_02_14_22_50_00_000.xlsx")
#feed=pd.read_excel(filepath+"BB4_2021_02_15_23_10_00_000.xlsx")
#feed=pd.read_excel(filepath+"BB5_2020_09_19_10_40_00_000.xlsx")
#feed=pd.read_excel(filepath+"BB6_2020_09_06_13_55_00_000.xlsx")
#feed=pd.read_excel(filepath+"BB7_2021_02_25_03_15_29_000.xlsx")
#feed=pd.read_excel(filepath+"BB9_2021_02_25_04_10_00_000.xlsx")
#feed=pd.read_excel(filepath+"BB10_2020_09_05_15_24_00_000.xlsx")
#feed=pd.read_excel(filepath+"BB12_2020_09_19_15_30_23_000.xlsx")
#feed=pd.read_excel(filepath+"BB13_2020_09_07_16_07_00_000.xlsx")
#feed=pd.read_excel(filepath+"BB27_2020_09_19_01_10_00_000.xlsx")
#feed=pd.read_excel(filepath+"BB31_2021_02_13_17_29_00_000.xlsx")
#feed=pd.read_excel(filepath+"BB34_2020_09_07_07_30_38_000.xlsx")
#feed=pd.read_excel(filepath+"BB39_2021_02_16_23_01_41_000.xlsx")
#feed=pd.read_excel(filepath+"BB41_2021_02_15_21_10_00_000.xlsx")
#feed=pd.read_excel(filepath+"BB57_2020_09_19_12_40_00_000.xlsx")
#feed=pd.read_excel(filepath+"BB66_2020_09_19_11_10_00_000.xlsx")
#feed=pd.read_excel(filepath+"BB67_2020_09_05_14_54_00_000.xlsx")
#feed=pd.read_excel(filepath+"BB70_2020_09_05_14_30_29_000.xlsx")
#feed=pd.read_excel(filepath+"BB88_2020_09_05_15_44_00_000.xlsx")
#feed=pd.read_excel(filepath+"BB106_2020_09_05_16_45_19_000.xlsx")
#feed1=pd.read_excel(filepath+"BB54_2021_02_15_00_45_00_000.xlsx")
#feed2=pd.read_excel(filepath+"BB54-1_2021_02_15_01_01_15_000.xlsx")

#filepath='/sugon/MaGroup/chandola/'
#feed1=pd.read_excel(filepath+"GCPair_A_1_2019_09_11_15_13_15_000.xlsx",engine='openpyxl')
#feed2=pd.read_excel(filepath+"GCPair_A_2-1_2019_09_11_15_43_55_000.xlsx",engine='openpyxl')
#feed3=pd.read_excel(filepath+"GCPair_A_2_2019_09_11_16_21_27_000.xlsx",engine='openpyxl')
#feed4=pd.read_excel(filepath+"GCPair_A_3-1_2019_09_11_16_05_00_000.xlsx",engine='openpyxl')
#feed1=pd.read_excel(filepath+"GCPair_B_1-1_2019_09_15_23_32_39_000.xlsx",engine='openpyxl')
#feed2=pd.read_excel(filepath+"GCPair_B_1_2019_09_15_23_00_00_000.xlsx",engine='openpyxl')
#feed3=pd.read_excel(filepath+"GCPair_B_2-2_2019_09_15_23_17_00_000.xlsx",engine='openpyxl')
#feed4=pd.read_excel(filepath+"GCPair_B_3-1_2019_09_15_23_25_00_000.xlsx",engine='openpyxl')


#x1=feed1['SDP_PhaPos_X']
#y1=feed1['SDP_PhaPos_Y']
#z1=feed1['SDP_PhaPos_Z']
#t1=feed1['SysTime']

#x2=feed2['SDP_PhaPos_X']
#y2=feed2['SDP_PhaPos_Y']
#z2=feed2['SDP_PhaPos_Z']
#t2=feed2['SysTime']

#x3=feed3['SDP_PhaPos_X']
#y3=feed3['SDP_PhaPos_Y']
#z3=feed3['SDP_PhaPos_Z']
#t3=feed3['SysTime']

#x4=feed4['SDP_PhaPos_X']
#y4=feed4['SDP_PhaPos_Y']
#z4=feed4['SDP_PhaPos_Z']
#t4=feed4['SysTime']

#x=np.concatenate((x1,x2))
#y=np.concatenate((y1,y2))
#z=np.concatenate((z1,z2))
#t=np.concatenate((t1,t2))



#x=np.concatenate((x1,x2,x3,x4))
#y=np.concatenate((y1,y2,y3,y4))
#z=np.concatenate((z1,z2,z3,z4))
#t=np.concatenate((t1,t2,t3,t4))


x=feed['SDP_PhaPos_X']
y=feed['SDP_PhaPos_Y']
z=feed['SDP_PhaPos_Z']
t=feed['SysTime']

#rang=feed1['SDP_AngleM']
MJD=[]
i=0
for el in t:
    time = Time(t[i],format='iso',scale='local')
    MJD.append(time.mjd-8/24)
    i=i+1 
print("Feed position file uploaded")    

#print(np.shape(MJD))
        
class feedposition:
    def __init__(self,T):
        self.T=T
    def altaz(self): 
        A=np.array(MJD)-self.T
        #i=np.where(np.min(A)) 
        l=np.where(abs(A)==np.min(abs(A)))
        #print(l)
        i=l[0][0]
        #print(i,self.T,MJD[i])
        R=np.sqrt(x[i]**2+y[i]**2+z[i]**2)
        z0=-z[i]
        y0=-x[i]
        x0=-y[i]
        az0 = 0.0
        if (np.abs(x0) < 1.0e-8):
           if (y0>0):
               az0 = 90.0
           else:
               az0 = 270.0
        else:
           tempaz0 = np.arctan(y0/x0)/convert
        if (x0>0 and y0>0):
           az0 = tempaz0
        if ((x0>0 and y0<0) or (x0>0 and y0==0)):
           az0 = tempaz0+360.0
        if (x0<0 and y0>0):
           az0 = tempaz0+180.0
        if ((x0<0 and y0<0) or (x0<0 and y0==0)):
           az0 = tempaz0+180.0
        el0 = np.arcsin(z0/R)/convert
        #print(az0,el0)    
      #  rag=rang[i]/convert
        return el0, az0

    #print(az0,el0) 
     #az.append(az0)
     #alt.append(el0)
     #time = Time(t[i],format='iso',scale='local',location=pos)
     #MJD.append(time.mjd-8/24)
     #time2=Time(time.mjd-8/24,format='mjd',scale='utc',location=pos)
    
     #time = Time('2016-09-15T15:32:39.000Z', format='isot',scale='utc', location = location)
     #time.ut1              
     #iers.IERS_A_URL = 'http://toshi.nofs.navy.mil/ser7/finals2000A.all'
     #download_IERS_A()
     #IERS-A default file name, URL, and ReadMe with content description
     #iers.IERS_A_FILE = 'finals2000A.all'
     #iers.IERS_A_URL = 'ftp://cddis.gsfc.nasa.gov/pub/products/iers/finals2000A.all'
     #iers.IERS_A_URL_MIRROR = 'https://datacenter.iers.org/data/9/finals2000A.all'
     #iers.IERS_A_README ='data/ReadMe.finals2000A'
     #iers.Conf.iers_auto_url.set('ftp://cddis.gsfc.nasa.gov/pub/products/iers/finals2000A.all')
     #time.delta_ut1_utc=8.0
     #lst_apparent = time2.sidereal_time('apparent')
     #coor=SkyCoord(alt=el0,az=az0, unit=(u.deg, u.deg),frame='altaz',obstime=time2,location=pos)
     #r=coor.transform_to(frame='fk5')
     #coord=SkyCoord(lst_apparent,0, unit=(u.hourangle, u.deg),frame='fk5')
     #coord2=coord.to_string('hmsdms')
     #r1=r.ra.hourangle
     #d=r.dec.deg
     #print(r1,d)
     #ra.append(r1)
     #dec.append(d)
     #i=i+1

#plt.plot(t[0:10],alt[0:10],'o')
#plt.plot(t[0:10],az[0:10],'o')
#plt.plot(t,alt,'o')
#plt.plot(t,az,'o')
#plt.plot(MJD,ra,'o')
#plt.plot(MJD,dec,'o')
#plt.plot(ra[500:1600],dec[500:1600],'o')
#plt.show()
