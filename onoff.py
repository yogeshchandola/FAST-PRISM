#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 16 20:04:42 2022

@author: chandola
"""
import time
import concurrent.futures
from astropy.io import fits
from astropy.modeling import models, fitting
from astropy.modeling.models import Polynomial1D
from astropy.modeling.models import Gaussian1D
from scipy.ndimage import gaussian_filter1d
from scipy.optimize import curve_fit
from scipy import fftpack, fft, ifft
import scipy.stats as st
import numpy as np
import os
os.environ['KMP_DUPLICATE_LIB_OK']='True'
import logging
from datetime import datetime
#import modules.coordinate as cor
#import modules.calibrate as cal
#import modules.createcube as cube
#import modules.dopset as dpset
#import modules.baselinesub as blsub
import process.modules.rearrangedata as rarr
#import process.modules.rfiremoval as fl
#import modules.plotspect as plt
from astropy.table import Table
import astropy.io.ascii as ast
import astropy.units as u
#from modules import coordinate
#from coordinate import ra, dec
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, AutoMinorLocator,IndexLocator, FormatStrFormatter
#import process.processF as pr
import process.processonoff as onoff
import spectral_cube
from spectral_cube import SpectralCube


#import matplotlib.cm as cm
#from matplotlib.ticker import MultipleLocator, AutoMinorLocator,IndexLocator, FormatStrFormatter
############################## process the OTF file #######################################
def sinf(x,A,B,C):  
    return A*np.sin(B*x+C)

def poly3D(x,A,B,C,D):
    return A*(x**3)+B*(x**2)+C*x+D

def poly2D(x,A,B,C):
    return A*(x**2)+B*x+C

def poly1D(x,A,B):
    return A*x+B
    
def Tsys(x,P0,P1,P2,n):
    return P0*np.arctan(np.sqrt(1+np.power(x,n))-P1)+P2
##################################################################################
start=time.perf_counter()
###################### Feed position file ##################################################
##################################################################################
filepath='/sugon/MaGroup/chandola/FAST_DATA/PT2022_0192/'
datapath='/sugon/MaGroup/chandola/FAST_DATA/PT2022_0192/'
tcalpath='/sugon/MaGroup/chandola/pipeline3.14/FASTpipeline/cube/process/modules/'
driftname='BB68'
driftdate=['20221107']
s=ast.read('beamsystemvar.txt',"r")
cycleontime=180
slewtime=30
#nstav=input("Number of time stamps to average:")
#nchav=input("Number of channels to average:")
ncalon=2
nstav=1
nchav=1
#sfreq=input("Start frequency:")
#subband=input("Subband size")
sfreq=1348
subband=15
#flag=input("Flagging ON or OFF:")
flag="OFF"
#avpol=input("Average two polarisations; ON or OFF")
avpol="ON"
#if(avpol=='ON'):
#    pol=0
if(avpol=='OFF'):
    print("No average of two polarisations")
    pol=input("Input the polarisation to show: 0 or 1")
contsub='ON'
flagfilepath='/sugon/MaGroup/chandola/FAST_DATA/PT2022_0192/'
flagfile=flagfilepath+driftname+'/flagfile.txt'
flags=ast.read(flagfile)
#maskfreq=[]
maskfreq=[[1352.0,1355.0]]
#badfreq=[[1372.0,13794.0],[1376.0,1378.0],[1393.0,1395.0]]
#badfreq=[[1372.0,1374.0],[1376.0,1382.0],[1393.0,1395.0]]
#badfreq=[[1393.0,1396.0],[1409.0,1413.0]]
##########################################################################################
#string='outdata_'+driftname+'_1_'+driftdate+'_'
string='outdata_'+driftname+'_'+driftdate[0]+'_'
#string='outdata_'+driftname+'_'
prefix=string+str(sfreq)+'_'+str(nstav)+'_'+str(nchav)
sfreq=float(sfreq)
efreq=sfreq+float(subband)
beamsize=2.9*(1420.405752/sfreq)
#maskarray=[]
rang=0.0
cores=1
#var=np.arange(19)
var=[0,7]
#var=[0,7,13]
#var=[0,1,2,3,4,5,6,13,10]
#######################################################################################################

def pprocess(arg): 
    i=arg[0]
    rang=arg[1]
    beamsize=arg[2]
    datapath=arg[3]
    tcalpath=arg[4]
    driftname=arg[5]
    driftdate=arg[6]
    nstav=arg[7]
    nchav=arg[8]
    flag=arg[9]
    
    t=onoff.process(i,rang,beamsize,datapath,tcalpath,driftname,driftdate,cycleontime,slewtime,ncalon,nstav,nchav,flag,sfreq,efreq,flags,maskfreq,contsub)
    
    #R,D,EL,AZ,RAN,Ton,Toff,Ta,Fluxon,Fluxoff,Flux,F,C,TON,TOFF,MJD,freq1,chw=t.process()
   # R,D,EL,AZ,Ta,Flux,mask,F,C,T, MJD, freq1,chw
    R,D,EL,AZ,Flux,F,V,C,T,MJD,freq1,chw=t.process()
    #t=pr.process(i,rang,beamsize,datapath,calpath,tcalpath,driftname,driftdate,nstav,nchav,flag,sfreq,efreq)
    #R,D,dAVFT,Flux,F,C,TON, MJD, freq1,chw=t.process()
    #output=[i,R,D,EL,AZ,Ta,Flux,mask,F,C,T, MJD, freq1,chw]
    output=[i,R,D,EL,AZ,Flux,F,V,C,T,MJD,freq1,chw]
    return output


arglist=[]
for i in var:
    fix=[rang,beamsize,datapath,tcalpath,driftname,driftdate,nstav,nchav,flag]
    fix.insert(0,i)
    print(fix)
    arglist.append(fix)

with concurrent.futures.ProcessPoolExecutor(max_workers=cores) as executor: 
     results=list(executor.map(pprocess,arglist))

####################################################################################################
j=0
for el in results:   
    i=results[j][0]
    R=results[j][1]
    D=results[j][2]
    EL=results[j][3]
    AZ=results[j][4]
    flux=np.array(results[j][5])  # Flux in mJy
  #  Ta=1000*np.array(results[j][6]) #Flux in mJy
  #  mask=results[j][7]
    F=results[j][6]
    V=results[j][7]
    C=results[j][8] 
    T=results[j][9]
    MJD=results[j][10]
    freq1=results[j][11]
    chw=results[j][12]
    Ta=flux
    mask=Ta < -1000000000
    Ta=np.ma.MaskedArray(Ta,mask=mask,fill_value=np.nan)
    #####################################################################################################################
    #u=int(len(TON)/4)
    #Tc=Ta.copy()
    #for k in range(4):
    #         if((k==3)|(k==1)):
    #             Ta_av=np.nanmedian(Tc[(k-1)*u:k*u],axis=0) 
    #             Ta[k*u:(k+1)*u]=Tc[k*u:(k+1)*u] -Ta_av
    #         else:   
    #             Ta_av=np.nanmedian(Tc[(k+1)*u:(k+2)*u],axis=0) 
    #             Ta[k*u:(k+1)*u]=Tc[k*u:(k+1)*u] -Ta_av
    ###### subtact off position data and average polarisations ###########################################
    
    ################################################################################
    if(avpol=='ON'):
        Ta= np.average(Ta[:,:,0:1],axis=2)
            #bvp3=rarr.rearr(ncalon,0,len(F),len(T),1,1,MJD,sfreq,chw,Ta)
            #Ta=bvp3.averagep() 
    elif(avpol=='OFF'):
        Ta=Ta[:,:,int(pol)]
     
    if(j==0): 
        coord=np.empty([len(results),len(R),2],dtype=float)
        elaz=np.empty([len(results),len(R),2],dtype=float)
        tsys=np.empty([len(results),len(R)],dtype=float)
        dataarray=np.empty([len(results),len(T),len(C)],dtype=float)
   #     maskarray=np.empty([len(Files),len(TON),len(C)])
        
    coord[j]=np.vstack([R,D]).T    
    elaz[j]=np.vstack([EL,AZ]).T
    tsys[j]=Tsys(90-np.array(EL), s['P0'][i], s['P1'][i], s['P2'][i], s['n'][i])
    dataarray[j]=Ta
   
    Z=Ta
   #################################################################################################################
    if(i < 9):
        beam='M0'+str(i+1)
    elif(i >=9):
        beam='M'+str(i+1)  
    
    X,Y=np.meshgrid(len(F),len(T))
    fig1, ax = plt.subplots()
    ax.set_title(beam)
    ax.set_xlabel('Channel number')
    ax.set_ylabel('Time stamp')
   # im = ax.imshow(Z,aspect='auto',cmap='rainbow',vmin=np.nanmin(Z),vmax=np.nanmax(Z),origin='lower')
    im = ax.imshow(Z,aspect='auto',cmap='rainbow',vmin=-40.0,vmax=80.0,origin='lower')
    fig1.colorbar(im,ax=ax,label='Flux (mJy)')
    if(avpol=='OFF'):
        plt.savefig(filepath+driftname+'/'+driftdate[0]+'/'+prefix+'M'+str(i+1)+'_'+str(pol)+'_colorplot.png')
    elif(avpol=='ON'):
        plt.savefig(filepath+driftname+'/'+driftdate[0]+'/'+prefix+'M'+str(i+1)+'_avpol'+'_colorplot.png')
    plt.show()
    ###################################################################################
    
    ###################################################################################
    fig2,ax = plt.subplots(figsize=(7, 7),nrows=1)
    ax.set_title(beam)
    ax.set_xlabel('Right Ascension (hrs)')
    ax.set_ylabel('Declination [deg.]')
    weights=np.nanmean(Ta.T,axis=0)
    #plt.scatter(np.array(R),np.array(D),cmap='rainbow',vmin=np.min(weights),vmax=np.max(weights),c=weights)
    plt.scatter(np.array(R),np.array(D),cmap='rainbow',vmin=-40.0,vmax=80.0,c=weights)
    ax.invert_xaxis()
    plt.legend(loc='upper left')
    plt.colorbar(label='Flux average in frequency (mJy)')
    if(avpol=='OFF'):
        plt.savefig(filepath+driftname+'/'+driftdate[0]+'/'+prefix+'M'+str(i+1)+'_'+str(pol)+'_avplot_COR.png')
    elif(avpol=='ON'):
        plt.savefig(filepath+driftname+'/'+driftdate[0]+'/'+prefix+'M'+str(i+1)+'_avpol'+'_avplot_COR.png')
    plt.show()
    ###########################################################################
    fig3, (ax1, ax2)=plt.subplots(figsize=(10, 7), nrows=2)
    ax1.plot(np.array(F),np.nanmean(Ta,axis=0), ls='-',drawstyle='steps',label='Flux average in time')
    ax2.plot(R,np.nanmean(Ta,axis=1),label='Flux average in frequency')
    ax1.legend(loc = 'upper left', frameon = False)
    ax2.legend(loc = 'upper left', frameon = False)
    ax2.invert_xaxis()
    ax1.set_title(beam)
    ax1.set_xlabel('Frequency (MHz)')
    ax2.set_xlabel('Right Ascension (hours)')
    if(avpol=='OFF'):
        plt.savefig(filepath+driftname+'/'+driftdate[0]+'/'+prefix+'M'+str(i+1)+'_'+str(pol)+'_avplots.png')
    elif(avpol=='ON'):
        plt.savefig(filepath+driftname+'/'+driftdate[0]+'/'+prefix+'M'+str(i+1)+'_avpol'+'_avplots.png')
    plt.show()
    

    ########################## Table out #########################################
    Data=Table([R,D], names=['RA','Dec'])
    ast.write(Data,filepath+driftname+'/'+driftdate[0]+'/'+prefix+'M'+str(i+1)+'_cor.dat',overwrite=True)
    Data2=Table([np.array(F),np.array(V),np.array(np.nanmean(Ta,axis=0))],names=['Frequency','Velocity','Fluxdensity'])
    ast.write(Data2,filepath+driftname+'/'+driftdate[0]+'/'+prefix+'M'+str(i+1)+'_spectra.dat',overwrite=True)
##########################################################################################
    j=j+1
    

######################################################################################
################################# select only drifting part ######################
nra=[]
d=0
for el in R:
    #print(d)
    #if((R[d] > R[d-1])&((R[d]-R[d-1])>0.00006)):
    nra.append(d)
    d=d+1
#nra=   
#######################################################################################
fig5,ax = plt.subplots(figsize=(7, 7),nrows=1)
weights=np.nanmean(dataarray.T,axis=0)
#plt.xlim(14.0,13.1)
ax.invert_xaxis()
#weights=np.average(dataarray.T,axis=0)
#plt.xlim(23.33,23.27)
plt.scatter(coord.T[0,nra,:],coord.T[1,nra,:],cmap='rainbow',c=weights[nra,:])
plt.colorbar(label='Flux (mJy)')
plt.xlabel('RA (hrs)')
plt.ylabel('DEC. (deg.)')
if(avpol=='OFF'):
    for pol in range(2):
        plt.savefig(filepath+driftname+'/'+driftdate[0]+'/'+prefix+str(pol)+'_avplot_COR.png')
elif(avpol=='ON'):
    plt.savefig(filepath+driftname+'/'+driftdate[0]+'/'+prefix+'avpol'+'_avplot_COR.png')
plt.show()


#hd.writeto(outdatapath+'coord.out.fits',overwrite=True)

####################################################################################
