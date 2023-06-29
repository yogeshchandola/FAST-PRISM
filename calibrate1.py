#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 29 14:41:25 2019

@author: yogesh
"""
import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

class cal:
    
    def __init__(self,beamname,obsdate,non,noff,nchan,F,ondata,offdata,Tcal_path):
            self.beamname=beamname
            self.obsdate=obsdate
            self.non=non
            self.noff=noff
            self.nchan=nchan
            self.F=F
            self.ondata=ondata
            self.offdata=offdata
            self.Tcal_path=Tcal_path
        
    
    def get_Tcal(beamname,obsdate,Tcal_path):
        # extract Tcal for the specific beam 
        beams = {'M01':0,'M02':1,'M03':2,'M04':3,'M05':4,'M06':5,'M07':6,'M08':7,'M09':8,'M10':9,
                 'M11':10,'M12':11,'M13':12,'M14':13,'M15':14,'M16':15,'M17':16,'M18':17,'M19':18}
    
        caltype='psr'
        if (int(obsdate) < 20200601):
        #if obsdate >= 20200601:
       # badbeam = ['M06','M14']
       # if beamname in badbeam:
            file = Tcal_path+ 'lowcal_20200531_'+caltype+'_tny.npz'
            tmp = np.load(file)
            freq = tmp['freq']
            tcals = tmp['tcal']
            tcal = tcals[:,:,beams[beamname]]
        #if (obsdate > 20200101) and (obsdate < 20200601):
        #else:
        elif(int(obsdate) >= 20200601):   
            file = Tcal_path+'lowcal_20201014_'+caltype+'_tny.npz'
            tmp = np.load(file)
            freq = tmp['freq']
            tcals = tmp['tcal']
            tcal = tcals[:,:,beams[beamname]]
        
        return freq, np.transpose(tcal)
         
    def interp(y,x,newx0,smooth=0):
        indx = np.where((newx0 > x[0]) & (newx0 < x[-1]))[0]
        newx = newx0[indx]
        #newy = y.copy()
        yshape = y.shape[0:-1] + newx0.shape
        #print(y.shape,yshape)
        newy = np.zeros(shape=yshape)
        if len(yshape) == 1:
            f = interp1d(x, y)
            newy[indx] = f(newx)
            if smooth != 0:
                newy = gaussian_filter1d(newy, smooth)
        if len(yshape) == 3:
            for i in range(yshape[0]):
                for j in range(yshape[1]):
                    f = interp1d(x, y[i,j,:])
                    tmpy = f(newx)
                    if smooth != 0:
                        tmpy = gaussian_filter1d(tmpy, smooth)
                    newy[i,j,indx] = tmpy
        if len(yshape) == 2:
            for i in range(yshape[0]):
                f = interp1d(x, y[i,:])
                tmpy = f(newx)
                if smooth != 0:
                    tmpy = gaussian_filter1d(tmpy, smooth)
                newy[i,indx] = tmpy
        return newy
    
    def tonoff(self): 
        pon=self.ondata
        poff=self.offdata
        obfreq=self.F
        freq, tcal= cal.get_Tcal(self.beamname,self.obsdate, self.Tcal_path)
        if(self.beamname=='M14'):
            tcal[0]=tcal[1]
        tcaln=cal.interp(tcal, freq, obfreq)   
        tcaln2=np.concatenate([tcaln,tcaln],axis=0)
        tcalt=gaussian_filter1d(tcaln2,201,axis=1).T
        #tcalt=np.average(tcaln2,axis=0)
        #tcalt=10.0*tcalt 
        #pcal=np.median(pon-poff,axis=0)
        pcal= np.nanmedian(pon,axis=0)-np.nanmedian(poff,axis=0)
        gshape=tcalt/gaussian_filter1d(pcal,201,axis=0)
        Ta=np.empty([self.noff,int(self.nchan),4],dtype=float)
        t=0  
        for el in Ta:  
            G=gshape
            Ta_off=poff[t]*G
            if(self.noff==self.non):
                    Ta_on=pon[t]*G-tcalt
                    Ta[t]=(Ta_on+Ta_off)/2
            else:
                Ta[t]=Ta_off
            t=t+1 

        return Ta,tcalt
    
    def Flux(self):
        t=cal(self.non,self.noff,self.nchan,self.ondata,self.offdata) 
        Ta=t.tonoff()
        Gain=0.59*25.2
        Flux=Ta/Gain
        return Flux
