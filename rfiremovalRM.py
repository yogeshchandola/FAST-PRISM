#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 27 16:06:04 2022

@author: chandola
"""
import scipy as sci
import numpy as np
import numpy.ma as ma
from scipy.ndimage import gaussian_filter1d
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
 
 #    if(i==3):
  #        Ta=Ta[:,:,1]
  #        Toff=Toff[:,:,1]
  #    if(i==4):
  #        Ta=Ta[:,:,0]
  #        Toff=Toff[:,:,0]
  #    elif(i!=3 and i!=4):
  #        bvp1=rarr.rearr(len(F),len(TOFF),1,1,MJD,freq1,chw,Ton)
  #        Ton=bvp1.averagep()
  #        bvp2=rarr.rearr(len(F),len(TOFF),1,1,MJD,freq1,chw,Toff)
  #        Toff=bvp2.averagep()

def flagbadbeam(beam,data):
        print("Flagging bad beam", beam,"\n")
        mask=data < -1000000000000000000000000000000000000000
        mask[:,:,:]=True
        x=np.ma.MaskedArray(data,mask=mask,fill_value=np.nan)
        return mask, x 

def flagbadpol(beam,data,pol):
        print("Flagging bad polarisation",pol,"for beam", beam,"\n")
        mask=data < -1000000000000000000000000000000000000000
        mask[:,:,:]=False
        mask[:,:,pol]=True
        x=np.ma.MaskedArray(data,mask=mask,fill_value=np.nan)
        return mask, x 


def flagbadchan1(beam,data,nchan,freq,badfreq):
        mask= data < -1000000000000000000000000000000000000000
        mask[:,:,:]=False
        i=0
        for el in badfreq:
            a=np.min(np.where(np.array(freq) > badfreq[i][0]))
            b=np.max(np.where(np.array(freq) < badfreq[i][1]))
            mask[:,a:b,:]=True
            i=i+1
        x=np.ma.MaskedArray(data,mask=mask,fill_value=np.nan)
        return mask, x

def flagbadchan2(beam,data,nchan,freq):
        mask= data < -1000000000000000000000000000000000000000
        mask[:,:,:]=False
        #print(mask.shape)
        #print(mask)
        Dav=np.median(data,axis=0)
########## smooth the spectrum using a Gaussian kernel ################################
        GausssmoothF=gaussian_filter1d(Dav,100,axis=0)
        a=Dav-GausssmoothF
        #fig,(ax1,ax2) = plt.subplots(figsize=(7, 7),nrows=2)
        #ax1.plot(freq,Dav[:,0:1],label=beam,color='green')
        #ax1.plot(freq,GausssmoothF[:,0:1],color='red',label='gausssmooth')
        #ax2.plot(freq,a[:,0:1],color='blue',label='residual')
        #ax1.legend(loc = 'upper left', frameon = False)
        #ax2.legend(loc = 'upper left', frameon = False)
        #plt.xlabel('Frequency (MHz)')
        #plt.ylabel('Ta (K)')
        #plt.show()
        #c=0    
        #w=300
        #u=int(nchan/3)
        #for el in range(3):
            #for it in range(1):
               # print(mask[0,:,0])
        #Ma=ma.MaskedArray(a,mask[0,:,:],fill_value=np.nan)
        med=np.nanmedian(np.nanmean(a[:,0:1],axis=1))
        mad=np.nanmedian(np.absolute(np.nanmean(a[:,0:1].data,axis=1)-med))
        limit1=med+(3*mad)
        limit2=med-(3*mad)
              #  print("Median:",med,"MAD:",mad,"for beam:",beam,". \n")
              #  print("Limits for beam", beam,'are', limit1, "and",limit2,"\n.")
                #print('med:',med,'mad:',mad,'ulimit:',limit1,'llimit:',limit2)
        t1=np.where(np.nanmean(a[:,0:1],axis=1) > limit1) 
        t2=np.where(np.nanmean(a[:,0:1],axis=1 ) < limit2)
        t=np.concatenate((t1[0],t2[0]),axis=0)
              #  fig2, ax3=plt.subplots(figsize=(7, 7),nrows=1)
              #  ax3.plot(freq,np.nanmean(Ma[:,0:1],axis=1),color='grey',label=str(el)+''+beam)
              #  ax3.axhline(y=limit1)
              #  ax3.axhline(y=limit2)
              #  ax3.legend(loc = 'upper left', frameon = False)
              #  plt.xlabel('Frequency (MHz)')
              #  plt.ylabel('Ta (K)')
                #print(t,'\n')
        if(len(t)>1):
                mask[:,t,:]=True
                #    ax3.axvline(x=freq[(c-1)*u+t])
              #  plt.show()
       #         c=c+1
        
       # w=10
       # k=0
       # for el in range(nchan):
       #     flgc=np.where(mask[0,k:k+w-1,0]==True)
       #     if(len(flgc[0])> 5):
       #         mask[:,k:k+w-1,:]=True
       #     k=k+1
    
        x=np.ma.MaskedArray(data,mask=mask,fill_value=np.nan)
        return mask, x

def flagbaddata(data,badtime,badchan,badpol):
     mask=data <-1000000000
     a=len(badtime)
     b=len(badchan)
     c=len(badpol)
     if(a>1 and b >1):
         i=0
         for el in range(int(len(badtime)/2)):
            j=0
            for el in range(int(len(badchan)/2)):
                k=0
                for el in range(int(len(badpol))):
                     mask[int(badtime[i]):int(badtime[i+1]),int(badchan[j]):int(badchan[j+1]),int(badpol[k])]=True
                     k=k+1
                j=j+2
            i=i+2
     elif(a==1 and b>1):
           j=0
           for el in range(int(len(badchan)/2)):
               k=0
               for el in range(int(len(badpol))):
                   if(int(badtime[0])==-1):
                       mask[:,int(badchan[j]):int(badchan[j+1]),int(badpol[k])]=True
                   else:
                       mask[int(badtime[0]),int(badchan[j]):int(badchan[j+1]),int(badpol[k])]=True
                   k=k+1
               j=j+2
     elif(a>1 and b==1):
        i=0
        for el in range(int(len(badtime)/2)):
               k=0
               for el in range(int(len(badpol))):
                   if(int(badchan[0])==-1):
                       mask[int(badtime[i]):int(badtime[i+1]),:,int(badpol[k])]=True
                   else:
                       mask[int(badtime[i]):int(badtime[i+1]),int(badchan[0]),int(badpol[k])]=True
                   k=k+1
               i=i+2 
     elif(a==1 and b==1):
              k=0
              for el in range(int(len(badpol))):
                  if(int(badchan[0])==-1 & int(badtime[0])==-1):
                      mask[:,:,int(badpol[k])]=True
                  else:
                      mask[int(badtime[0]),int(badchan[0]),int(badpol[k])]=True
                  k=k+1
            
            
     x=np.ma.MaskedArray(data,mask=mask,fill_value=np.nan)
     return mask,x
 
    
 
    
def flagsfile(flags,beam,i,Ta):
    a=0
    badchanArray=[]
    for el in flags['badbeam']: 
                  #print(el[a])
                  #a=np.where(flags['badbeam']==badbeam) 
                  #badchan=[]
                  badbeam=str(flags['badbeam'][a]).split(',') 
                  print(badbeam)
                  if(beam in badbeam):
                      badscan=str(flags['badscan'][a]).split(',')
                      print(badscan)
                      if(str(i) in badscan):
                        #a=np.where(flags['badbeam']==beam)
                        #b=np.where(badscan[a[0]]==i)
                            badtime=flags['badtime'][a]#[b[0]]
                            badtime=str(badtime).split(',')
                            badchan=flags['badchan'][a]#[b[0]]
                            badchan=str(badchan).split(',')
                            badpol=flags['badpol'][a]#[b[0]]
                            badpol=str(badpol).split(',')
                            print('Flagged',badtime,badchan,badpol)
                            mask, mTa=flagbaddata(Ta,badtime,badchan,badpol)
                      elif('-1' in badscan):
                       # a=np.where(flags['badbeam']==beam)     
                            badtime=flags['badtime'][a]#[a[0]]
                            badtime=str(badtime).split(',')
                            badchan=flags['badchan'][a]#[a[0]]
                            badchan=str(badchan).split(',')
                            badpol=flags['badpol'][a]#[a[0]]
                            badpol=str(badpol).split(',')
                            print('Flagged',badtime,badchan,badpol)
                            mask, mTa=flagbaddata(Ta,badtime,badchan,badpol)           
                  elif('-1' in badbeam):
                      #a=np.where(flags['badbeam']==badbeam)
                      badscan=str(flags['badscan'][a]).split(',')
                      if(str(i) in badscan):
                       # a=np.where(flags['badbeam']==-1)
                       # b=np.where(flags['badscan'][a[0]]==i)
                            badtime=flags['badtime'][a]#[b[0]]
                            badtime=str(badtime).split(',')
                            badchan=flags['badchan'][a]#[b[0]]
                            badchan=str(badchan).split(',')
                            badpol=flags['badpol'][a]#[b[0]]
                            badpol=str(badpol).split(',')
                            print('Flagged',badtime,badchan,badpol)
                            mask, mTa=flagbaddata(Ta,badtime,badchan,badpol)
                        
                      elif('-1' in badscan):
                           #a=np.where(flags['badbeam']==-1) 
                           badtime=flags['badtime'][a]#[a[0]]
                           badtime=str(badtime).split(',')
                           badchan=flags['badchan'][a]#[a[0]]
                           badchan=str(badchan).split(',')
                           badpol=flags['badpol'][a]#[a[0]]
                           badpol=str(badpol).split(',')
                           print('Flagged',badtime,badchan,badpol)
                           mask, mTa=flagbaddata(Ta,badtime,badchan,badpol)
                  print(badchan)
                  if(a==0):
                      badchanArray=badchan
                      badtimeArray=badtime
                  else:    
                      badchanArray=np.concatenate([badchanArray,badchan])
                      badtimeArray=np.concatenate([badtimeArray,badtime])
                  a=a+1
    return mask,mTa,badchanArray,badtimeArray