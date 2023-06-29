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
#from modules import coordinate
#from coordinate import ra, dec
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, AutoMinorLocator,IndexLocator, FormatStrFormatter
#import process.processF as pr
import process.processMBOTFa as prmbotf
import process.processCRAFTS as prcrafts
import spectral_analysis as spectral
import gc
import functools


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
mode='mbotf'
filepath='/sugon/MaGroup/chandola/FAST_DATA/3057/'
datapath='/sugon/MaGroup/chandola/FAST_DATA/3057/'
tcalpath='/sugon/MaGroup/chandola/pipeline3.14/FASTpipeline/cube/process/modules/'
driftname='GCPair_A'
driftdate=['20190911']

#filepath='/sugon/MaGroup/chandola/PT2020_0079/'
#datapath='/sugon/PT2020_0079/'
#calpath='/sugon/PT2020_0079/CAL/'
#tcalpath='/sugon/MaGroup/chandola/pipeline3.14/FASTpipeline/cube/process/modules/'
#driftname='R2d3'
#driftdate=['20200925','20200929','20201007']
#driftdate=['20200924']
s=ast.read('beamsystemvar.txt',"r")

#nstav=input("Number of time stamps to average:")
#nchav=input("Number of channels to average:")
nstav=1
nchav=1
#sfreq=input("Start frequency:")
#subband=input("Subband size")
sfreq=1365
subband=5
#flag=input("Flagging ON or OFF:")
flag="OFF"
#avpol=input("Average two polarisations; ON or OFF")
avpol="ON"
#if(avpol=='ON'):
#    pol=0
if(avpol=='OFF'):
    print("No average of two polarisations")
    pol=input("Input the polarisation to show: 0 or 1")
contsub='OFF'
flagfilepath='/sugon/MaGroup/chandola/FAST_DATA/3057/GCPair_A/20190911/'
flagfile=flagfilepath+'flagfile.txt'
flags=ast.read(flagfile)
badfreq=[]
#lffreq=[1363.1,1363.5,1365.5,1365.8] # line free frequency for cont. sub.
lffreq=[1366.1,1366.7,1368.1,1368.8]
#lffreq=[1371.1,1372.0,1373.0,1373.4]
#badfreq=np.array([[1312.0,1315.0],[1323.0,1326.0],[1328.0,1330.0],[1339.0,1341.0]])
#badfreq=[[1339.0,1341],[1345.0,1347.0],[1356.0,1358.0]]
#badfreq=[[1372.0,1374.0],[1376.0,1378.0],[1393.0,1395.0]]
#badfreq=[[1372.0,1374.0],[1376.0,1382.0],[1393.0,1395.0]]
#badfreq=[[1393.0,1396.0],[1409.0,1413.0]]
##########################################################################################
#string='outdata_'+driftname+'_M01_'+driftdate[0]+'_'
string='outdata_'+driftname+'_1_'+driftdate[0]+'_'
#string='outdata_'+driftname+'_'
prefix=string+str(sfreq)+'_'+str(nstav)+'_'+str(nchav)
sfreq=float(sfreq)
efreq=sfreq+float(subband)
beamsize=2.9*(1420.405752/sfreq)
#maskarray=[]
rang=23.4
cores=1

#var=np.arange(19)
#var=[0,4,13,10]
var=[11]
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
    if(mode=='mbotf'):
        t=prmbotf.process(i,rang,beamsize,datapath,tcalpath,driftname,driftdate,nstav,nchav,flag,sfreq,efreq,flags,badfreq,contsub,lffreq)
    elif(mode=='crafts'):
        t=prcrafts.process(i,rang,beamsize,datapath,calpath,tcalpath,driftname,driftdate,nstav,nchav,flag,sfreq,efreq,flags,badfreq,contsub,lffreq)
    #R,D,EL,AZ,RAN,Ton,Toff,Ta,Fluxon,Fluxoff,Flux,F,C,TON,TOFF,MJD,freq1,chw=t.process()
    R,D,EL,AZ,Ta,Flux,mask,F,C,T, MJD, freq1,chw=t.process()
    #t=pr.process(i,rang,beamsize,datapath,calpath,tcalpath,driftname,driftdate,nstav,nchav,flag,sfreq,efreq)
    #R,D,dAVFT,Flux,F,C,TON, MJD, freq1,chw=t.process()
    output=[i,R,D,EL,AZ,Ta,Flux,mask,F,C,T, MJD, freq1,chw]
    
    return output


arglist=[]
for i in var:
    fix=[rang,beamsize,datapath,tcalpath,driftname,driftdate,nstav,nchav,flag]
    fix.insert(0,i)
    print(fix)
    arglist.append(fix)

with concurrent.futures.ProcessPoolExecutor(max_workers=cores) as executor: 
     results=list(executor.map(pprocess,arglist))
     gc.collect() 
    # All objects collected
     objects = [i for i in gc.get_objects() 
           if isinstance(i, functools._lru_cache_wrapper)]
    # All objects cleared
     for object in objects:
        object.cache_clear()
     

####################################################################################################
j=0
for el in results:   
    i=results[j][0]
    R=results[j][1]
    D=results[j][2]
    EL=results[j][3]
    AZ=results[j][4]
    Ta0=results[j][5]   # Flux in K
    Ta=1000*np.array(results[j][6]) #Flux in mJy
    mask=results[j][7]
    F=results[j][8]
    C=results[j][9] 
    T=results[j][10]
    MJD=results[j][11]
    freq1=results[j][12]
    chw=results[j][13]
    
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
            bvp3=rarr.rearr(len(T),0,len(F),len(T),1,1,MJD,sfreq,chw,Ta)
            Ta=bvp3.averagep() 
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
    ax2.set_xlabel('Right Ascension (hrs)')
    if(avpol=='OFF'):
        plt.savefig(filepath+driftname+'/'+driftdate[0]+'/'+prefix+'M'+str(i+1)+'_'+str(pol)+'_avplots.png')
    elif(avpol=='ON'):
        plt.savefig(filepath+driftname+'/'+driftdate[0]+'/'+prefix+'M'+str(i+1)+'_avpol'+'_avplots.png')
    plt.show()
    
    
   
    
    ########################## Table out #########################################
    Data=Table([R,D], names=['RA','Dec'])
    ast.write(Data,filepath+driftname+'/'+driftdate[0]+'/'+prefix+'M'+str(i+1)+'_cor.dat',overwrite=True)
##########################################################################################
    gc.collect()
  
    # All objects collected
    objects = [i for i in gc.get_objects() 
           if isinstance(i, functools._lru_cache_wrapper)]
    # All objects cleared
    for object in objects:
        object.cache_clear()
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


####################################################################################
scale=np.empty([len(results),len(R)],dtype=float)
i=0
for el in tsys:
    scale[i]=tsys[i]/tsys[0]
    i=i+1
    
###################################################################################
fig6,ax = plt.subplots(figsize=(7, 7),nrows=1)
ax.invert_xaxis()
weights=tsys.T
#plt.xlim(23.33,23.27)
plt.scatter(coord.T[0,nra,:],coord.T[1,nra,:],cmap='rainbow',c=weights[nra,:])
plt.colorbar(label='Tsys (K)')

plt.xlabel('RA (hrs)')
plt.ylabel('DEC. (deg.)')
plt.savefig(filepath+driftname+'/'+driftdate[0]+'/'+prefix+'avpol'+'_Tsys_COR.png')
plt.show()

####################################################################################
fig7,ax = plt.subplots(figsize=(7, 7),nrows=1)
ax.invert_xaxis()
weights=scale.T
#plt.xlim(23.33,23.27)
plt.scatter(coord.T[0,nra,:],coord.T[1,nra,:],cmap='rainbow',c=weights[nra,:])
plt.colorbar(label='Scale')
plt.xlabel('RA (hrs)')
plt.ylabel('DEC. (deg.)')
plt.savefig(filepath+driftname+'/'+driftdate[0]+'/'+prefix+'avpol'+'_Scale_COR.png')
plt.show()
####################################################################################
fig8,ax = plt.subplots(figsize=(7, 7),nrows=1)
ax.invert_xaxis()
dav=np.nanmean(dataarray.T,axis=0)
weights=dav/scale.T
#plt.xlim(23.33,23.27)
plt.scatter(coord.T[0,nra,:],coord.T[1,nra,:],cmap='rainbow',c=weights[nra,:])
plt.colorbar(label='Flux (mJy) corrected')

plt.xlabel('RA (hrs)')
plt.ylabel('DEC. (deg.)')
plt.savefig(filepath+driftname+'/'+driftdate[0]+'/'+prefix+'avpol'+'_Ta_corrected_COR.png')
plt.show()
#####################################################################################
fig9,ax = plt.subplots(figsize=(7, 7),nrows=1)
ax.invert_xaxis()
weights=elaz.T[0,:,:]
#plt.xlim(23.33,23.27)
plt.scatter(coord.T[0,nra,:],coord.T[1,nra,:],cmap='rainbow',c=weights[nra,:])
plt.colorbar(label='Elev.(deg.)')
plt.xlabel('RA (hrs)')
plt.ylabel('DEC. (deg.)')
plt.savefig(filepath+driftname+'/'+driftdate[0]+'/'+prefix+'avpol'+'_ELEV_COR.png')
plt.show()
#####################################################################################
fig10,ax = plt.subplots(figsize=(7, 7),nrows=1)
ax.invert_xaxis()
weights=90-elaz.T[0,:,:]
#plt.xlim(23.33,23.27)
plt.scatter(coord.T[0,nra,:],coord.T[1,nra,:],cmap='rainbow',c=weights[nra,:])
plt.colorbar(label='ZEN.(deg.)')
plt.xlabel('RA (hrs)')
plt.ylabel('DEC. (deg.)')
plt.savefig(filepath+driftname+'/'+driftdate[0]+'/'+prefix+'avpol'+'_ZEN_COR.png')
plt.show()
#####################################################################################
fig11,ax = plt.subplots(figsize=(7, 7),nrows=1)
ax.invert_xaxis()
weights=elaz.T[1,:,:]
#plt.xlim(23.33,23.27)
plt.scatter(coord.T[0,nra,:],coord.T[1,nra,:],cmap='rainbow',c=weights[nra,:])
plt.colorbar(label='Az.(deg.)')
plt.xlabel('RA (hrs)')
plt.ylabel('DEC. (deg.)')
plt.savefig(filepath+driftname+'/'+driftdate[0]+'/'+prefix+'avpol'+'_AZ_COR.png')
plt.show()


fig12,ax = plt.subplots(figsize=(7, 7),nrows=1)
tarr=np.arange(len(elaz[0,:,0]))+1
a1=len(elaz[0,:,0])/3
plt.scatter(tarr,90-elaz.T[0,:,0])
plt.axvline(x=a1,color='b',label='Scan 1 end')
plt.axvline(x=a1*2,color='r',label='Scan 2 end')
plt.xlabel('Time stamps')
plt.ylabel('Zenith angle (deg.)')
plt.legend()
plt.savefig(filepath+driftname+'/'+driftdate[0]+'/'+prefix+'ZEN_vs_time_.png')
plt.show()


###################################################################################
#if(flag=='ON'):
#    finalarray=np.ma.MaskedArray(dataarray,maskarray,fillvalue=np.nan)
#elif(flag=='OFF'):
    
finalarray=dataarray
#from astropy.wcs import WCS
#from astropy.utils.misc import NumpyRNGContext
import cygrid

if(mode=='mbotf'):
######################## Rearrange the cube ###################################    
    ramin= np.nanmin(coord[:,:,0])#min RA
    ramax= np.nanmax(coord[:,:,0])#max RA  
    decmin= np.nanmin(coord[:,:,1])# min dec
    decmax= np.nanmax(coord[:,:,1])# max dec

#hd.writeto(outdatapath+'coord.out.fits',overwrite=True)
##################### regridding ##############################################
    


    beamsize_fwhm=beamsize/60
    pixsize = beamsize_fwhm/3.0
    dnaxis1 = int(np.round((ramax-ramin)*15 / pixsize))
    dnaxis2 = int(np.round((decmax-decmin) / pixsize))
    dnaxis3= len(C)
#dnaxis3=len(F)

    hd1=fits.PrimaryHDU()
    hd1.header['NAXIS']= 3
    hd1.header['NAXIS1']= dnaxis1
    hd1.header['NAXIS2']= dnaxis2
    hd1.header['NAXIS3']= dnaxis3
    hd1.header['CTYPE1']= 'RA---TAN'
    hd1.header['CTYPE2']= 'DEC--TAN'
    hd1.header['CTYPE3']= 'freq'
    hd1.header['CUNIT1']= 'deg'
    hd1.header['CUNIT2']= 'deg'
    hd1.header['CUNIT3']= 'Hz'
    hd1.header['CDELT1']= -pixsize
    hd1.header['CDELT2']= pixsize
    hd1.header['CDELT3']= (F[1]-F[0])*1.0E+06
  #  hd1.header['CDELT3']= chw*1.0E+06
    hd1.header['CRPIX1']= 1
    hd1.header['CRPIX2']= 1
    hd1.header['CRPIX3']= 1
    hd1.header['CRVAL1']= ramax*15
    hd1.header['CRVAL2']= decmin
    hd1.header['CRVAL3']= F[0]*1.0E+06
#hd1.header['CRVAL3']=F[0]*1.0E+06
    
    gridder = cygrid.WcsGrid(hd1.header) # 

    kernelsize_fwhm = beamsize_fwhm/2 # degrees
    kernelsize_sigma = kernelsize_fwhm/np.log(8*np.sqrt(2))
    sphere_radius = 3*kernelsize_sigma
    PA=0.0
    gridder.set_kernel('gauss2d',(kernelsize_sigma, kernelsize_sigma, PA),sphere_radius,kernelsize_sigma/2.0)
#gridder.set_kernel('gauss2d',(0.62,0.62,0.0),1.85,0.31)
    ra=[]
    dec=[]
    signal=[]
    i=0
    for el in coord[0,nra,0]:
#    print(i)
        j=0
        for el in coord[:,0,0]:
            ra.append(15*coord[j,i,0])
            dec.append(coord[j,i,1])
            sig=np.ma.masked_values(finalarray[j,i,:], value=np.nan)
            signal.append(sig)
#        print(j)
            j=j+1
        i=i+1    

    ra1=np.array(ra,dtype=float) 
    dec1=np.array(dec,dtype=float)
    signal1=np.array(signal,dtype=float)       
        
#signal2=np.array([signal,signal])
    gridder.grid(ra1, dec1, signal1)
    head=gridder.get_header()
    unweightdata=gridder.get_unweighted_datacube()
    weights=gridder.get_weights()
    gridded_map = gridder.get_datacube()


#######################################cube##############################################
    hd2=fits.PrimaryHDU()
    hd2.data=gridded_map
    OBJECT  = driftname          #/ Name of object                                 
    hd2.header['telescop'] = 'FAST'
    hd2.header['INSTRUME']= '19beam'          # / Instrument/detector used          	
    hd2.header['OBSERVER']= 'Yogesh'           #/ Astronomical observer                          
    hd2.header['DATE-OBS']	= '2019-09-11'    #     / Date (yyyy-mm-dd) of observations              
    hd2.header['BUNIT   ']	= 'mJy/beam'       #    / Physical units of image  
    hd2.header['BMAJ']=beamsize_fwhm
    hd2.header['BMIN']=beamsize_fwhm
    hd2.header['BPA']=0.0                      
    hd2.header['EPOCH   ']	=            2.00E+03 #/ Celestial coordiate equinox                    
    hd2.header['ALTRPIX ']	=         1.000000E+00 #/ Frequency alternate reference pixel                 
#hd2.header['DATAMAX ']	=  np.max(gridded_map) #/ Maximum physical value in image                
#hd2.header['DATAMIN ']	=  np.min(gridded_map) 
    hd2.header['CTYPE1']= 'RA---TAN'
    hd2.header['CTYPE2']= 'DEC--TAN'
    hd2.header['CTYPE3']= 'FREQ'
    hd2.header['CUNIT1']= 'deg'
    hd2.header['CUNIT2']= 'deg'
    hd2.header['CUNIT3']= 'Hz'
    hd2.header['CDELT1']= -pixsize
    hd2.header['CDELT2']= pixsize
    hd2.header['CDELT3']= (F[1]-F[0])*1.0E+06
 #   hd2.header['CDELT3']= chw*1.0E+06
    hd2.header['CRPIX1']= 1
    hd2.header['CRPIX2']= 1
    hd2.header['CRPIX3']= 1
    hd2.header['CROTA1']=      0.000000E+00 #/ Axis coordinate rotation                       
    hd2.header['CROTA2']=      0.000000E+00                                                  
    hd2.header['CROTA3']=      0.000000E+00   
    hd2.header['CRVAL1']= ramax*15
    hd2.header['CRVAL2']= decmin
    hd2.header['CRVAL3']= F[0]*1.0E+06
    if(contsub=='ON'):
        if(avpol=='ON'):
            hd2.writeto(filepath+driftname+'/'+driftdate[0]+'/'+prefix+'avpol'+'_0_contsub_cube.fits',overwrite=True)
        elif(avpol=='OFF'):
            hd2.writeto(filepath+driftname+'/'+driftdate[0]+'/'+prefix+str(pol)+'_0_contsub_cube.fits',overwrite=True)
    elif(contsub=='OFF'):
        if(avpol=='ON'):
            hd2.writeto(filepath+driftname+'/'+driftdate[0]+'/'+prefix+'avpol'+'_0_cube.fits',overwrite=True)
        elif(avpol=='OFF'):
            hd2.writeto(filepath+driftname+'/'+driftdate[0]+'/'+prefix+str(pol)+'_0_cube.fits',overwrite=True)

########################################map###############################################
    target_wcs = gridder.get_wcs()

    fig = plt.figure(figsize=(7, 7))
    ax = fig.add_subplot(111, projection=target_wcs.celestial)
#im = ax.imshow(    np.nanmean(gridded_map,axis=0), cmap='rainbow', aspect='auto',vmin=23,vmax=30, origin='lower', interpolation='nearest')
    im = ax.imshow(
        np.nanmean(gridded_map,axis=0), cmap='rainbow', aspect='auto',
        origin='lower', interpolation='nearest'
        )
#ax.invert_xaxis()
    fig.colorbar(im,ax=ax,label='Flux (mJy)')
    ra1, dec1 = ax.coords
    ra1.set_axislabel('R.A. [h:m:s]')
    dec1.set_axislabel('Dec [deg]')
    if(contsub=='ON'):
        if(avpol=='OFF'):
            plt.savefig(filepath+driftname+'/'+driftdate[0]+'/'+prefix+'_contsub_colormap.png')
        elif(avpol=='ON'):
            plt.savefig(filepath+driftname+'/'+driftdate[0]+'/'+prefix+'avpol'+'_contsub_colormap.png')
    elif(contsub=='OFF'):
        if(avpol=='OFF'):
            plt.savefig(filepath+driftname+'/'+driftdate[0]+'/'+prefix+'_colormap.png')
        elif(avpol=='ON'):
            plt.savefig(filepath+driftname+'/'+driftdate[0]+'/'+prefix+'avpol'+'_colormap.png')
    plt.show()


    hd3=fits.PrimaryHDU()
    hd3.data=np.nanmean(gridded_map,axis=0)
    OBJECT  = driftname           #/ Name of object                                 
    hd3.header['telescop'] = 'FAST'
    hd3.header['INSTRUME']= '19beam'          # / Instrument/detector used          	
    hd3.header['OBSERVER']= 'Yogesh'           #/ Astronomical observer                          
    hd3.header['DATE-OBS']	= '2019-09-11'    #     / Date (yyyy-mm-dd) of observations              
    hd3.header['BUNIT   ']	= 'mJy/beam'       #    / Physical units of image  
    hd3.header['BMAJ']=beamsize_fwhm
    hd3.header['BMIN']=beamsize_fwhm
    hd3.header['BPA']=0.0                      
    hd3.header['EPOCH   ']	=            2.00E+03 #/ Celestial coordiate equinox                    
    hd3.header['ALTRPIX ']	=         1.000000E+00 #/ Frequency alternate reference pixel                 
#hd2.header['DATAMAX ']	=  np.max(gridded_map) #/ Maximum physical value in image                
#hd2.header['DATAMIN ']	=  np.min(gridded_map) 
    hd3.header['CTYPE1']= 'RA---TAN'
    hd3.header['CTYPE2']= 'DEC--TAN'
#hd2.header['CTYPE3']= 'FREQ'
    hd3.header['CUNIT1']= 'deg'
    hd3.header['CUNIT2']= 'deg'
#hd2.header['CUNIT3']= 'Hz'
    hd3.header['CDELT1']= -1*pixsize
    hd3.header['CDELT2']= pixsize
#hd2.header['CDELT3']= (F[1]-F[0])*1.0E+06
    hd3.header['CRPIX1']= 1
    hd3.header['CRPIX2']= 1
#hd2.header['CRPIX3']= 1
    hd3.header['CROTA1']=      0.000000E+00 #/ Axis coordinate rotation                       
    hd3.header['CROTA2']=      0.000000E+00                                                  
#hd2.header['CROTA3']=      0.000000E+00   
    hd3.header['CRVAL1']= ramax*15
    hd3.header['CRVAL2']= decmin
#hd2.header['CRVAL3']= F[0]*1.0E+06
    if(contsub=='ON'):
        if(avpol=='ON'):
            hd3.writeto(filepath+driftname+'/'+driftdate[0]+'/'+prefix+'avpol'+'_contsub_collapsed.fits',overwrite=True)
            spectral.spectralanalysis(filepath, driftname, driftdate, prefix,0)
        elif(avpol=='OFF'):
            hd3.writeto(filepath+driftname+'/'+driftdate[0]+'/'+prefix+str(pol)+'_contsub_collapsed.fits',overwrite=True)
    elif(contsub=='OFF'):
        if(avpol=='ON'):
            hd3.writeto(filepath+driftname+'/'+driftdate[0]+'/'+prefix+'avpol'+'_collapsed.fits',overwrite=True)
        elif(avpol=='OFF'):
            hd3.writeto(filepath+driftname+'/'+driftdate[0]+'/'+prefix+str(pol)+'_collapsed.fits',overwrite=True)


##########################################################
####################


####################################################################################################
#################################################################################



if(mode=='crafts'):
    ramin= np.min(coord[:,:,0])#min RA
    ramax= np.max(coord[:,:,0])#max RA  
    
    radiff=ramax-ramin
    if(radiff < 0.15):
          N=1
    elif(radiff > 0.15):
          N=int(radiff/0.15)+1
    for ele in range(N):
            ramin= np.min(coord[:,:,0])#min RA
            ramax= np.max(coord[:,:,0])#max RA 
            if(ramin+(ele+1)*0.15 < ramax):
                index=np.where((coord[:,:,0]> ramin+ele*0.15)& (coord[:,:,0]< (ele+1)*0.15+ramin ))
            elif(ramin+(ele+1)*0.15 > ramax):
                index=np.where((coord[:,:,0]> ramin+ele*0.15)& (coord[:,:,0]< ramax ))
            print(index)
            fig,ax = plt.subplots(figsize=(7, 7),nrows=1)
            ax.invert_xaxis()   
            plt.plot(coord[:,index[1],0],coord[:,index[1],1],'o')
            plt.xlabel("R.A. (hrs)")
            plt.ylabel("Dec. (deg.)")
            plt.show()
            
           
######################## Rearrange the cube ###################################    
            ramin= np.min(coord[:,index[1],0])#min RA
            ramax= np.max(coord[:,index[1],0])#max RA  
            decmin= np.min(coord[:,index[1],1])# min dec
            decmax= np.max(coord[:,index[1],1])# max dec

#hd.writeto(outdatapath+'coord.out.fits',overwrite=True)
##################### regridding ##############################################
            beamsize_fwhm=beamsize/60.0
            pixsize = beamsize_fwhm/3.0
            dnaxis1 = int(np.round((ramax-ramin)*15 / pixsize))
            dnaxis2 = int(np.round((decmax-decmin) / pixsize))
            dnaxis3= len(C)
            hd1=fits.PrimaryHDU()
            hd1.header['NAXIS']= 3
            hd1.header['NAXIS1']= dnaxis1
            hd1.header['NAXIS2']= dnaxis2
            hd1.header['NAXIS3']= dnaxis3
            hd1.header['CTYPE1']= 'RA---TAN'
            hd1.header['CTYPE2']= 'DEC--TAN'
            hd1.header['CTYPE3']= 'FREQ'
            hd1.header['CUNIT1']= 'deg'
            hd1.header['CUNIT2']= 'deg'
            hd1.header['CUNIT3']= 'Hz'
            hd1.header['CDELT1']= -1*pixsize
            hd1.header['CDELT2']= pixsize
            hd1.header['CDELT3']= (F[1]-F[0])*1.0E+06
            #hd1.header['CDELT3']= chw*1.0E+06
            hd1.header['CRPIX1']= 1
            hd1.header['CRPIX2']= 1
            hd1.header['CRPIX3']= 1
            hd1.header['CRVAL1']= ramax*15
            hd1.header['CRVAL2']= decmin
            hd1.header['CRVAL3']= F[0]*1.0E+06
#hd1.header['CRVAL3']=F[0]*1.0E+06
            gridder = cygrid.WcsGrid(hd1.header) #
            kernelsize_fwhm = beamsize_fwhm/2 # degrees
            kernelsize_sigma = kernelsize_fwhm/np.log(8*np.sqrt(2))
            sphere_radius = 3*kernelsize_sigma
            PA=0.0
           # print(kernelsize_sigma)
            kernel_args=('gauss2d',(kernelsize_sigma, kernelsize_sigma,PA),sphere_radius,kernelsize_sigma/2.0)
            gridder.set_kernel(*kernel_args)
           # gridder.set_kernel('gauss1d',(kernelsize_sigma,),sphere_radius,kernelsize_sigma/2.0)
#gridder.set_kernel('gauss2d',(0.62,0.62,0.0),1.85,0.31)
            ra=[]
            dec=[]
            signal=[]
            i=0
            for el in index[1]:
                j=0
                for el in coord[:,0,0]:
                    ra.append(15*coord[j,index[1][i],0])
                    dec.append(coord[j,index[1][i],1])
                    sig=np.ma.masked_values(finalarray[j,index[1][i],:], value=np.nan)
                    signal.append(sig)
                    j=j+1
                i=i+1    

            ra1=np.array(ra,dtype=float) 
            dec1=np.array(dec,dtype=float)
            signal1=np.array(signal,dtype=float)       
        
#signal2=np.array([signal,signal])
            gridder.grid(ra1, dec1, signal1)
            head=gridder.get_header()
            unweightdata=gridder.get_unweighted_datacube()
            weights=gridder.get_weights()
            gridded_map = gridder.get_datacube()


#######################################cube##############################################
            hd2=fits.PrimaryHDU()
            hd2.data=gridded_map
            OBJECT  =driftname           #/ Name of object                                 
            hd2.header['telescop'] = 'FAST'
            hd2.header['INSTRUME']= '19beam'          # / Instrument/detector used          	
            hd2.header['OBSERVER']= 'Yogesh'           #/ Astronomical observer                          
            hd2.header['DATE-OBS']	= '2020-09-14' #driftdate    #     / Date (yyyy-mm-dd) of observations              
            hd2.header['BUNIT   ']	= 'mJy/beam'       #    / Physical units of image  
            hd2.header['BMAJ']=beamsize/60
            hd2.header['BMIN']=beamsize/60
            hd2.header['BPA']=0.0                      
            hd2.header['EPOCH   ']	=            2.00E+03 #/ Celestial coordiate equinox                    
            hd2.header['ALTRPIX ']	=         1.000000E+00 #/ Frequency alternate reference pixel                 
#hd2.header['DATAMAX ']	=  np.max(gridded_map) #/ Maximum physical value in image                
#hd2.header['DATAMIN ']	=  np.min(gridded_map) 
            hd2.header['CTYPE1']= 'RA---TAN'
            hd2.header['CTYPE2']= 'DEC--TAN'
            hd2.header['CTYPE3']= 'FREQ'
            hd2.header['CUNIT1']= 'deg'
            hd2.header['CUNIT2']= 'deg'
            hd2.header['CUNIT3']= 'Hz'
            hd2.header['CDELT1']= -pixsize
            hd2.header['CDELT2']= pixsize
            hd2.header['CDELT3']= (F[1]-F[0])*1.0E+06
           # hd2.header['CDELT3']= chw*1.0E+06
            hd2.header['CRPIX1']= 1
            hd2.header['CRPIX2']= 1
            hd2.header['CRPIX3']= 1
            hd2.header['CROTA1']=      0.000000E+00 #/ Axis coordinate rotation                       
            hd2.header['CROTA2']=      0.000000E+00                                                  
            hd2.header['CROTA3']=      0.000000E+00   
            hd2.header['CRVAL1']= ramax*15
            hd2.header['CRVAL2']= decmin
            hd2.header['CRVAL3']= F[0]*1.0E+06
            if(contsub=='ON'):
                if(avpol=='ON'):
                    hd2.writeto(filepath+driftname+'/'+driftdate[0]+'/'+prefix+'avpol_'+str(ele)+'_contsub_cube.fits',overwrite=True)
                elif(avpol=='OFF'):
                    hd2.writeto(filepath+driftname+'/'+driftdate[0]+'/'+prefix+str(pol)+'_'+str(ele)+'_contsub_cube.fits',overwrite=True)
            elif(contsub=='OFF'):
                if(avpol=='ON'):
                    hd2.writeto(filepath+driftname+'/'+driftdate[0]+'/'+prefix+'avpol_'+str(ele)+'_cube.fits',overwrite=True)
                elif(avpol=='OFF'):
                    hd2.writeto(filepath+driftname+'/'+driftdate[0]+'/'+prefix+str(pol)+'_'+str(ele)+'_cube.fits',overwrite=True)

########################################map####################################str###########
            target_wcs = gridder.get_wcs()

            fig = plt.figure(figsize=(7, 7))
            ax = fig.add_subplot(111, projection=target_wcs.celestial)
            #im = ax.imshow(
            #np.nanmean(gridded_map,axis=0), cmap='rainbow', aspect='auto',vmin=np.min(np.nanmean(gridded_map,axis=0)),vmax=np.max(np.nanmean(gridded_map,axis=0)),
            #origin='lower', interpolation='nearest')
            im = ax.imshow(
            np.nanmean(gridded_map,axis=0), cmap='rainbow', aspect='auto',origin='lower', interpolation='nearest')
            fig.colorbar(im,ax=ax,label='Flux (mJy)')
            ra1, dec1 = ax.coords
            ra1.set_axislabel('R.A. [h:m:s]')
            dec1.set_axislabel('Dec [deg]')
            if(contsub=='ON'):
                if(avpol=='OFF'):
                    plt.savefig(filepath+driftname+'/'+driftdate[0]+'/'+prefix+'_'+str(ele)+'_contsub_colormap.png')
                elif(avpol=='ON'):
                    plt.savefig(filepath+driftname+'/'+driftdate[0]+'/'+prefix+'avpol_'+str(ele)+'_contsub_colormap.png')
            elif(contsub=='OFF'):
                if(avpol=='OFF'):
                    plt.savefig(filepath+driftname+'/'+driftdate[0]+'/'+prefix+'_'+str(ele)+'_colormap.png')
                elif(avpol=='ON'):
                    plt.savefig(filepath+driftname+'/'+driftdate[0]+'/'+prefix+'avpol_'+str(ele)+'_colormap.png')
            plt.show()

            hd3=fits.PrimaryHDU()
            hd3.data=np.nanmean(gridded_map,axis=0)
            OBJECT  = driftname          #/ Name of object                                 
            hd3.header['telescop'] = 'FAST'
            hd3.header['INSTRUME']= '19beam'          # / Instrument/detector used          	
            hd3.header['OBSERVER']= 'Yogesh'           #/ Astronomical observer                          
            hd3.header['DATE-OBS']  = '2020-09-14' #     / Date (yyyy-mm-dd) of observations              
            hd3.header['BUNIT   ']	= 'mJy/beam'       #    / Physical units of image  
            hd3.header['BMAJ']=beamsize/60
            hd3.header['BMIN']=beamsize/60
            hd3.header['BPA']=0.0                      
            hd3.header['EPOCH   ']	=            2.000E+03 #/ Celestial coordiate equinox                    
            hd3.header['ALTRPIX ']	=         1.000000E+00 #/ Frequency alternate reference pixel                 
#hd2.header['DATAMAX ']	=  np.max(gridded_map) #/ Maximum physical value in image                
#hd2.header['DATAMIN ']	=  np.min(gridded_map) 
            hd3.header['CTYPE1']= 'RA---TAN'
            hd3.header['CTYPE2']= 'DEC--TAN'
#hd2.header['CTYPE3']= 'FREQ'
            hd3.header['CUNIT1']= 'deg'
            hd3.header['CUNIT2']= 'deg'
#hd2.header['CUNIT3']= 'Hz'
            hd3.header['CDELT1']= -pixsize
            hd3.header['CDELT2']= pixsize
#hd2.header['CDELT3']= (F[1]-F[0])*1.0E+06
            hd3.header['CRPIX1']= 1
            hd3.header['CRPIX2']= 1
#hd2.header['CRPIX3']= 1
            hd3.header['CROTA1']=      0.000000E+00 #/ Axis coordinate rotation                       
            hd3.header['CROTA2']=      0.000000E+00                                                  
#hd2.header['CROTA3']=      0.000000E+00   
            hd3.header['CRVAL1']= ramax*15
            hd3.header['CRVAL2']= decmin
#hd2.header['CRVAL3']= F[0]*1.0E+06
            if(contsub=='ON'):
                if(avpol=='ON'):
                    hd3.writeto(filepath+driftname+'/'+driftdate[0]+'/'+prefix+'avpol_'+str(ele)+'_contsub_collapsed.fits',overwrite=True)
                    spectral.spectralanalysis(filepath, driftname, driftdate, prefix,ele)
                elif(avpol=='OFF'):
                    hd3.writeto(filepath+driftname+'/'+driftdate[0]+'/'+prefix+str(pol)+'_'+str(ele)+'_contsub_collapsed.fits',overwrite=True)
            elif(contsub=='OFF'):
                if(avpol=='ON'):
                    hd3.writeto(filepath+driftname+'/'+driftdate[0]+'/'+prefix+'avpol_'+str(ele)+'_collapsed.fits',overwrite=True)
                elif(avpol=='OFF'):
                    hd3.writeto(filepath+driftname+'/'+driftdate[0]+'/'+prefix+str(pol)+'_'+str(ele)+'_collapsed.fits',overwrite=True)
            
gc.collect() 
    # All objects collected
objects = [i for i in gc.get_objects() 
           if isinstance(i, functools._lru_cache_wrapper)]
    # All objects cleared
for object in objects:
        object.cache_clear()               

##########################################################
####################


 #if (flag=='ON'):
 #    x=np.ma.MaskedArray(Ta[:,SC[0]],mask=mask,fill_value=np.nan)
 #    maskarray.append(mask)
 #    Z=x
 #elif(flag=='OFF'):
 #    Z=Ta[:,SC[0]]

 #################################################################################
# if (flag=='ON'):
#     fig4, ax4 = plt.subplots()
#     im = ax4.imshow(mask,aspect='auto',cmap='rainbow',origin='lower')
 #im = ax.imshow(Z,aspect='auto',cmap='seismic',origin='lower')
 #im = ax.imshow(Z,aspect='auto',cmap='Spectral',origin='lower')
 #im = ax.imshow(Z,aspect='auto',cmap='RdYlBu',origin='lower')
 #im = ax.imshow(Z,aspect='auto',cmap='bwr',origin='lower')
 #im = ax.imshow(Z,aspect='auto',cmap='gist_gray',origin='lower')
 #plt.legend(loc='upper right')
 #    fig3.colorbar(im,ax=ax4)
 #    if(avpol=='OFF'):
 #        plt.savefig(filepath+driftname+'/'+driftdate+'/'+prefix+'M'+str(i+1)+'_'+str(pol)+'_mask_colorplot.png')
 #    elif(avpol=='ON'):
 #        plt.savefig(filepath+driftname+'/'+driftdate+'/'+prefix+'M'+str(i+1)+'_avpol'+'_mask_colorplot.png')
 #    plt.show()
 
 ############################### fitting linear equation ##############################################################
 #l_init = Linear1D(slope=0.32,intercept=5.6) #initial parameters
 #fit_l=fitting.LinearLSQFitter()
 #fit_l = fitting.LevMarLSQFitter()
 #fit_l = fitting.SLSQPLSQFitter()
 #L = fit_l(l_init, logSM3,metal3)
 #param_cov=fit_l.fit_info['param_cov']
 #print L.slope,L.intercept #, param_cov
 #plt.plot(sm,L(sm),'-b')
 #print(Flux.shape,Flux1.shape)
