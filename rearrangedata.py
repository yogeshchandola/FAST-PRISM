#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Nov  2 19:59:28 2019

@author: yogesh
"""

import numpy as np

class rearr:
    def __init__(self, ncalon,si,nchan,nstamp,nchav,nstav,MJD,freq1,chw,data): # constructor
            self.ncalon=ncalon
            self.si=si
            self.nchav=nchav
            self.nstav=nstav
            self.nchan=nchan    
            self.nstamp=nstamp
            self.MJD=MJD
            self.freq1=freq1
            self.chw=chw
            self.data=data

    def ONOFFt(self):     
        TCALON=[]
        TCALOFF=[]
        i=0
        for el in range(int(self.nstamp)):
            if((self.si+i)%self.ncalon ==0):
                    TCALON.append(i)
            elif((self.si+i)%self.ncalon!=0):
                    TCALOFF.append(i)
            i=i+1
        return TCALON,TCALOFF    
    
    def ONOFF2(self):
        TCALON=[]
        TCALOFF=[]
        i=0
        for el in range(int(self.nstamp)):
            if(i%2 ==0):
                    TCALON.append(i)
            elif(i%2!=0):
                    TCALOFF.append(i)
            i=i+1
        return TCALON,TCALOFF 

    def dcalonoff(self):
        t=rearr(self.ncalon,self.si,self.nchan,self.nstamp,self.nchav,self.nstav,self.MJD,self.freq1,self.chw,self.data)
        TCALON,TCALOFF= t.ONOFFt() 
        #print(TCALON,TCALOFF) 
        #print(self.data)
        DCALON= self.data[TCALON]
        DCALOFF= self.data[TCALOFF]  
        MJDON=self.MJD[TCALON]
        MJDOFF=self.MJD[TCALOFF]
        return DCALON, DCALOFF, MJDON,MJDOFF        
                       
######################### Average data along time ############################
    def averaget(self):
           #print(self.nstamp,self.nstav)
           nst=int(self.nstamp)
           nstav=int(self.nstav)          
           ft=nst/nstav
           print(ft)
           AVT=np.empty([int(ft),int(self.nchan),4],dtype=float)
           T=[]
           i=0
           a=0
           b=nstav
           for el in range(int(ft)):
               if(nstav >1):
                   avt= np.nanmean(self.data[i:int(i+nstav-1)],axis=0)
               elif(nstav==1):
                   avt=self.data[i]
               AVT[a]=avt
               t=self.MJD[b]
               T.append(t)
              # print(i,b)
               i=i+nstav 
               b=b+2*nstav     
               #print(i)
               #print(i+nstav)
               a=a+1                
           return AVT, T 

        #######################################################################
    def averaget2(self):
            #print(self.nstamp,self.nstav)
            nst=int(self.nstamp)
            nstav=int(self.nstav)          
            ft=nst/nstav
           # print(ft)
            AVT=np.empty([int(ft),int(self.nchan),2],dtype=float)
            T=[]
            i=0
            a=0
            b=nstav
            for el in range(int(ft)):
                if(nstav >1):
                    avt= np.nanmean(self.data[i:int(i+nstav-1)],axis=0)
                elif(nstav==1):
                    avt= self.data[i]
                AVT[a]=avt
                t=self.MJD[i]
                T.append(t)
                i=i+nstav 
                #print(i)
                #print(i+nstav)
                a=a+1                
            return AVT, T 
##########################Average along channels/frequency ###################
    def averagec(self):
           nc=int(self.nchan)
           nchav=int(self.nchav)
           fc=nc/nchav         
           AVC=np.empty([int(self.nstamp),int(fc),4],dtype=float)
           #print AVC.shape
           F=[]
           C=[]
           b=0
           for el in range(int(self.nstamp)): 
                 i=0
                 a=0
                 for el in range(int(fc)):
                   if(nchav >1):
                       avc= np.nanmean(self.data[b,i:int(i+nchav-1),:],axis=0)
                   elif(nchav==1):
                       avc=self.data[b,i]
                   #print avc
                   AVC[b,a,:]=avc
                   i=i+nchav 
                   a=a+1
                   if(b==0):         
                      f=self.freq1+i*self.chw 
                      C.append(a)
                      F.append(f)
                 b=b+1    
           return AVC, F, C   
   
    def averagec2(self):
           nc=int(self.nchan)
           nchav=int(self.nchav)
           fc=nc/nchav         
           AVC=np.empty([int(self.nstamp),int(fc),2],dtype=float)
           #print AVC.shape
           F=[]
           C=[]
           b=0
           for el in range(int(self.nstamp)): 
                 i=0
                 a=0
                 for el in range(int(fc)):
                   if(nchav >1):  
                       avc= np.nanmean(self.data[b,i:int(i+nchav-1),:],axis=0)
                   elif(nchav==1):
                       avc=self.data[b,i]
                   #print avc
                   AVC[b,a,:]=avc
                   i=i+nchav 
                   a=a+1
                   if(b==0):         
                      f=self.freq1+i*self.chw 
                      C.append(a)
                      F.append(f)
                 b=b+1    
           return AVC, F, C  
######################### Along polarization for XX,YY ######################     
    def averagep(self):
           return np.nanmean(self.data[:,:,0:1],axis=2)
           #return (self.data[:,:,0]+self.data[:,:,1])/2
           #return self.data[:,:,0]
##############################################################################
               
       
        
