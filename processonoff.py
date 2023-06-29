from astropy.io import fits
import numpy as np
import numpy.ma as ma
import process.modules.coordinate as cor
import process.modules.calibrate1 as cal
#import modules.createcube as cube
#import modules.dopset as dpset
#import modules.baselinesub as blsub
import process.modules.rearrangedata as rarr
import process.modules.rfiremovalRM as rfi
import process.modules.continuum_sub as cnsub
import process.modules.freqtovelocity as fqv
import process.modules.dopset as dps
from scipy.optimize import curve_fit
from scipy.ndimage import gaussian_filter1d
from scipy.ndimage import median_filter
#import modules.plotspect as plt
#from modules import coordinate
#from coordinate import ra, dec
import astropy.io.ascii as ast
from astropy.time import Time
#import modules.plotspect as plt
#from modules import coordinate
#from coordinate import ra, dec
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, AutoMinorLocator,IndexLocator, FormatStrFormatter
import os
os.environ['KMP_DUPLICATE_LIB_OK']='True'
import logging
from datetime import datetime
import pandas as pd
from astropy.table import Table
from scipy.ndimage import gaussian_filter1d
#####################################################################
Gain0=np.array([15.98,	16.27,	16.48,	16.1,	16.12,	15.94,	15.98,	16.02,	15.71])
filepath='/sugon/MaGroup/chandola/pipeline3.14/FASTpipeline/cube/'
file=filepath+'FASTgains.xlsx'
gf=pd.read_excel(file,engine='openpyxl')
gf1=gf['f1050']*Gain0[0]
gf2=gf['f1100']*Gain0[1]
gf3=gf['f1150']*Gain0[2]
gf4=gf['f1200']*Gain0[3]
gf5=gf['f1250']*Gain0[4]
gf6=gf['f1300']*Gain0[5]
gf7=gf['f1350']*Gain0[6]
gf8=gf['f1400']*Gain0[7]
gf9=gf['f1450']*Gain0[8]

def appeffi(ZEN,a,b,c):
    if(ZEN < 26.4):
        eta=a*ZEN+b
    elif(ZEN > 26.4):
        eta=c*ZEN+b+26.4*(a-c)
    return eta

class process:
    def __init__(self,b,rang,beamsize,datapath,tcalpath,driftname,driftdate,cycleontime,slewtime,ncalon,nstav,nchav,flag,sfreq,efreq,flags,maskfreq,contsub):
        self.b=b
        self.rang=rang
        self.beamsize=beamsize
        self.datapath=datapath
        self.tcalpath=tcalpath
        self.driftname=driftname
        self.driftdate=driftdate
        self.cycleontime=cycleontime
        self.slewtime=slewtime
        self.ncalon=ncalon
        self.nstav=nstav
        self.nchav=nchav
        self.flag=flag
        self.sfreq=sfreq
        self.efreq=efreq
        self.flags=flags
        self.maskfreq=maskfreq
        self.contsub=contsub
       # self.lffreq=lffreq
    def process(self):
        data_path=self.datapath
        sfreq=self.sfreq
        efreq=self.efreq
        flags=self.flags
        maskfreq=self.maskfreq
        ontime=self.cycleontime
        slewtime=self.slewtime
        if(self.b < 9):
            beam='M0'+str(self.b+1)
        elif(self.b>=9):
            beam='M'+str(self.b+1)   
        print("Processing for beam:",beam,"\n")
        os.system('echo'+' file > '+beam+'tmp.txt')
        rp=0
        for date in self.driftdate:
            FI=data_path + self.driftname + '/' + date + '/' + self.driftname + '_onoff-'+beam+'_W_*'+'.fits'
            
            if(rp==0):
                os.system('ls '+ FI+' > '+beam+'tmp1.txt')
                os.system('cat '+beam+'tmp.txt '+beam+'tmp1.txt > '+beam+'files.txt')
                r=ast.read(beam+"files.txt","r")
                file1=r['file']
                File=file1
            else:
                os.system('ls '+ FI+' > '+beam+'tmp1.txt')
                os.system('cat '+beam+'tmp.txt '+beam+'tmp1.txt > '+beam+'files.txt')             
                s= ast.read(beam+"files.txt","r")
                file2=s['file']
                File=np.concatenate((File,file2),axis=0)
            rp=rp+1
       
        i=0
        nstamp=0
        for files in File:
            inpf=File[i]
            print("Processing file",inpf,"\n")
            file=fits.open(inpf)
            header=file[1].header
            #print("The header of the file is as given below")
            #print "+++++++++++++++++++++++++++++++++++++++++"
            #print(header)
            #data1=file[1].data
           # print("Shape of data structure: time stamps, channels, polarizations")
           # print(data1['DATA'].shape)
           #header=file[1].header
           #print "The header of the file is as given below"
           #print "+++++++++++++++++++++++++++++++++++++++++"
           #print header
            data=file[1].data.field('DATA')
            #print("Shape of data structure: time stamps, channels, polarizations")
            #print(data.shape)
#################### Extract the information ########################
            OBSTYPE=file[1].data.field('OBSTYPE')
            nstamp1=len(file[1].data.field('EXPOSURE'))
            freq1=np.average(file[1].data.field('FREQ'))
            exp=np.average(file[1].data.field('EXPOSURE'))
            chw=np.average(file[1].data.field('CHAN_BW'))
            nchan=np.average(file[1].data.field('NCHAN'))
#           mjd=np.average(file[1].data.field('UTOBS'))
            MJD=file[1].data.field('UTOBS')
       #     RA=file[1].data.field('OBJ_RA')
            #DEC=file[1].data.field('OBJ_DEC')
            dateobs=file[1].data.field('DATE-OBS')
            
            X=file[1].header['OBSGEO-X']
            Y=file[1].header['OBSGEO-Y']
            Z=file[1].header['OBSGEO-Z']
            freq=np.arange(nchan)*chw+freq1
            SC=np.where(( np.array(freq)> sfreq) & (np.array(freq) < efreq))
            ###################### Combine data from all FITS ##################################
            nstamp=nstamp+nstamp1
            if(i==0):
                data1=data[:,SC[0],:]
                MJD1=MJD
            elif(i>0):    
                data1=np.concatenate((data1,data[:,SC[0],:]),axis=0)
                MJD1=np.concatenate((MJD1,MJD),axis=0)
            i=i+1    
        ##########################determine the     
        i=0
        R=[]
        D=[]
        EL=[]
        AZ=[]
        ZEN=[]
        T=[]
        for el in range(len(MJD1)):
            test=cor.coordinate(self.b,self.beamsize,self.rang,X,Y,Z,MJD1[i])
            el, az, ra,dec=test.cor() 
            EL.append(el)
            AZ.append(az)
            ZEN=abs(90-np.array(el))
            
            #print(i,ra,dec)
                    #Correction for precession
                    #RA = RA(2000) + (3.075 + 1.336 * sin(RA) * tan(Dec)) * y
                    #Dec(2000) + 20.04 * cos(RA) * y 
                    #Correction for nutation
                    #delta RA = (0.9175 + 0.3978 * sin(RA) * tan(Dec)) * dL - cos(RA) * tan(Dec) * dE
                    #delta Dec = 0.3978 * cos(RA) * dL + sin(RA) * dE
                    #print ra, dec
            R.append(ra)
            D.append(dec)
            T.append(i)
            i=i+1    
        fig,(ax1,ax2)=plt.subplots(figsize=(10,7),nrows=2) 
        ax1.set_title(beam)
        ax1.plot(T,R,'o')
        ax2.plot(R,D,'o')
        plt.show()    
         
       
        #ontime=270
        #exp=1
        cycle=2*int(np.round((ontime+slewtime)/exp)) # on +off time 30 sec extra for switching over
        Ncycle=int(np.round((nstamp)/cycle))
        print("exp:",exp,"nstamp:",nstamp,"stcycle:",cycle,"Ncycle:",Ncycle)
        dataout=[]
        dataoffout=[]
        dataonout=[]
        
        for i in range(Ncycle):  
        ##############################################################################
        ##################--separate the calon and caloff--###########################
        ##############################################################################
                print("cycle no:",i)
                if((i+1)*cycle > nstamp):  
                    data2=data1[i*cycle:nstamp-1,:,:]
                    mjd=MJD1[i*cycle:nstamp-1] 
                    Nst=nstamp-1-i*cycle
                    fig,(ax1,ax2)=plt.subplots(figsize=(10,7),nrows=2) 
                    ax1.set_title(beam+'_'+str(i))
                    ax1.plot(T[i*cycle:nstamp-1],R[i*cycle:nstamp-1],'o')
                    ax2.plot(R[i*cycle:nstamp-1],D[i*cycle:nstamp-1],'o')
                    plt.show() 
                else:
                    data2=data1[i*cycle:(i+1)*cycle,:,:]
                    mjd=MJD1[i*cycle:(i+1)*cycle]
                    Nst=cycle
                    fig,(ax1,ax2)=plt.subplots(figsize=(10,7),nrows=2) 
                    ax1.set_title(beam+'_'+str(i))
                    ax1.plot(T[i*cycle:(i+1)*cycle],R[i*cycle:(i+1)*cycle],'o')
                    ax2.plot(R[i*cycle:(i+1)*cycle],D[i*cycle:(i+1)*cycle],'o')
                    plt.show() 
                print("Nst:",Nst)
                si=i*cycle
        ########################### separate cal ON and OFF data ######
                #t=rarr.rearr(int(np.round(self.ncalon/exp)),nchan,Nst,1,1,mjd,freq1,chw,data2)
                t=rarr.rearr(int(self.ncalon),si,nchan,Nst,1,1,mjd,freq1,chw,data2)
                Non,Noff=t.ONOFFt()
                calon,caloff,mjdon,mjdoff=t.dcalonoff()   
                print("calon:",len(Non),"caloff:",len(Noff))
      ####################Flagging calon and caloff data#########
              # calon=Dcalon[:,SC[0],:]
              # caloff=Dcaloff[:,SC[0],:] 
               
             #  xmajorLocator = MultipleLocator(100)
             #  ymajorLocator = MultipleLocator(100)
             #  xminorLocator = AutoMinorLocator(5)
             #  yminorLocator = AutoMinorLocator(5)


               ###################################################################################
                Pcalon=np.nanmedian(calon,axis=0)
                Pcaloff=np.nanmedian(caloff,axis=0)
                pcal=gaussian_filter1d(Pcalon-Pcaloff,201,axis=0)
                fig3,(ax1,ax2,ax3) = plt.subplots(figsize=(7, 7),nrows=3)
                ax1.set_title(beam+'_'+str(i))
                ax1.plot(freq[SC[0]],Pcalon[:,0],color='blue', label='calon pol 0')
                ax1.plot(freq[SC[0]],Pcaloff[:,0],color='red',label='caloff pol 0')
                ax2.plot(freq[SC[0]],Pcalon[:,1],color='blue',label='calon pol 1')
                ax2.plot(freq[SC[0]],Pcaloff[:,1],color='red',label='caloff pol 1')
                ax3.plot(freq[SC[0]],1.2*Pcaloff[:,0]/pcal[:,0],color='green',label='pol 0')
                ax3.plot(freq[SC[0]],1.2*Pcaloff[:,1]/pcal[:,1],color='orange',label='pol 1')
                ax1.legend(loc = 'upper left', frameon = False)
                ax2.legend(loc = 'upper left', frameon = False)
                ax3.legend(loc = 'upper left', frameon = False)
                plt.xlabel('Frequency (MHz)')
                ax1.set_ylabel('Median Power')
                ax2.set_ylabel('Median Power')
                ax3.set_ylabel('1.2*median(Pcaloff)/pcal')
            
                plt.show()
                plt.savefig(data_path + self.driftname + '/' + self.driftdate[0] + '/' + self.driftname + '_onoff-'+str(sfreq)+'_'+str(beam)+'_'+str(i)+'_power.png')
                Data=Table([freq[SC[0]],np.nanmean(Pcaloff[:,0:1],axis=1)], names=['Frequency','Power'])
                ast.write(Data,data_path + self.driftname + '/' + self.driftdate[0] + '/' + self.driftname + '_onoff-'+str(sfreq)+'_'+str(beam)+'_'+str(i)+'_caloff_power.dat',overwrite=True)
                
               ###################################################################################
              # flag=self.flag
              # if(flag=='ON'):
              #      print("RFI removal before calibration for beam:",beam)
                  # maskon,mcalon=rfi.flagbadchan2(beam,calon,len(SC[0]),np.array(freq[SC[0]]))
                  # if(len(Non)==len(Noff)):
                  #     mcaloff=ma.MaskedArray(caloff,mask=maskon,fill_value=np.nan)
                  # elif(len(Non)!= len(Noff)):
                  #     flgc=np.where(np.array(maskon[0,:,0])==True)
                      # print(flgc[0])
                  #     maskoff=caloff < -10000000000000000000000000000
                  #     maskoff[:,flgc[0],:]=True
                  #     mcaloff=ma.MaskedArray(caloff,mask=maskoff,fill_value=np.nan)
                  # elif(beam in flags['badbeam']):
             #      if(beam in flags['badbeam']):  
             #          if(i in flags['badscan']):
             #              a=np.where(flags['badbeam']==beam)
             #              b=np.where(flags['badscan'][a[0]]==i)
             #              badtime=flags['badtime'][b[0]]
             #              badtime=str(badtime[0]).split(',')
             #              badchan=flags['badchan'][b[0]]
             #              badtime=str(badtime[0]).split(',')
             #              badpol=flags['badpol'][b[0]]
             #              badtime=str(badtime[0]).split(',')
             #              maskon, mcalon=rfi.flagbaddata(calon,badtime,badchan,badpol)
                          # mcaloff=ma.MaskedArray(caloff,mask=maskon,fill_value=np.nan)
             #              maskoff,mcaloff=rfi.flagbaddata(caloff,badtime,badchan,badpol)
                         #  mcalon=ma.MaskedArray(mcalon,mask=maskoff,fill_value=np.nan)
                       
             #  elif(flag=='OFF'):
             #      maskon= calon < -100000000000000000000000000
             #      maskon[:,:,:,]=False
             #      mcalon=ma.MaskedArray(calon,mask=maskon,fill_value=np.nan)
             #      maskoff= caloff < -100000000000000000000000000
             #      maskoff[:,:,:,]=False
             #      mcaloff=ma.MaskedArray(caloff,mask=maskoff,fill_value=np.nan)
               ################## Plot flagged data  #########
               
                fig1, (ax1, ax2) = plt.subplots(figsize=(10, 7), ncols=2)
            #   ax1.xaxis.set_major_locator(xmajorLocator)
            #   ax1.yaxis.set_major_locator(ymajorLocator)
            #   ax1.xaxis.set_minor_locator(xminorLocator)
            #   ax1.yaxis.set_minor_locator(yminorLocator)

            #   ax2.xaxis.set_major_locator(xmajorLocator)
            #   ax2.yaxis.set_major_locator(ymajorLocator)
            #   ax2.xaxis.set_minor_locator(xminorLocator)
            #   ax2.yaxis.set_minor_locator(yminorLocator)

            #   ax3.xaxis.set_major_locator(xmajorLocator)
            #   ax3.yaxis.set_major_locator(ymajorLocator)
            #   ax3.xaxis.set_minor_locator(xminorLocator)
            #   ax3.yaxis.set_minor_locator(yminorLocator)

            #   ax4.xaxis.set_major_locator(xmajorLocator)
            #   ax4.yaxis.set_major_locator(ymajorLocator)
            #   ax4.xaxis.set_minor_locator(xminorLocator)
            #   ax4.yaxis.set_minor_locator(yminorLocator)

   #### XX ##################
   #XX = ax1.imshow(dataXX, cmap='RdBu', aspect='auto', interpolation='none')


   #XX = ax1.imshow(davpON, cmap='seismic', aspect='auto', interpolation='none',origin='lower')
                XX = ax1.imshow(calon[:,:,0], cmap='rainbow', aspect='auto', interpolation='none',origin='lower')
                fig1.colorbar(XX, ax=ax1,label='Power')
                ax1.set_ylabel('Time stamp')
                ax1.set_xlabel('Channel')
                string=beam+'_'+str(i)+' CALON pol0'
                ax1.set_title(string)
   #YY = ax2.imshow(dataYY, cmap='RdBu', aspect='auto', interpolation='none')
   #YY = ax2.imshow(davpOFF, cmap='seismic', aspect='auto', interpolation='none',origin='lower')
                YY = ax2.imshow(calon[:,:,1], cmap='rainbow', aspect='auto', interpolation='none',origin='lower')
                fig1.colorbar(YY, ax=ax2, label='Power')
                ax2.set_ylabel('Time stamp')
                ax2.set_xlabel('Channel')
                string=beam+'_'+str(i)+' CALON pol1'
                ax2.set_title(string)
                plt.show()

                fig2, (ax1, ax2) = plt.subplots(figsize=(10, 7), ncols=2)

                AA = ax1.imshow(caloff[:,:,0], cmap='rainbow', aspect='auto', interpolation='none',origin='lower')
   #AA = ax3.imshow(flaggedata, cmap='seismic', aspect='auto', interpolation='none',origin='lower')
                fig2.colorbar(AA, ax=ax1, label='Power')
                ax1.set_ylabel('Time stamp')
                ax1.set_xlabel('Channel')
                string=beam+'_'+str(i)+' CALOFF pol0'
                ax1.set_title(string)

                AA = ax2.imshow(caloff[:,:,1], cmap='rainbow', aspect='auto', interpolation='none',origin='lower')
   #AA = ax4.imshow(flaggedata2, cmap='seismic', aspect='auto', interpolation='none',origin='lower')
                fig2.colorbar(AA, ax=ax2, label='Power')
                ax2.set_ylabel('Time stamp')
                ax2.set_xlabel('Channel')
                string=beam+'_'+str(i)+' CALOFF pol1'
                ax2.set_title(string)

             #  plt.savefig(prefix+'colorscale.png')
                plt.show()
                #########################Separate ON and OFF position data for calon#####################################################################
                #TONcalon=int((len(Non)-(30/exp)*(1/self.ncalon))) #subtracted extra 30 stamps as stabilizing time
                #TONcalonst=int(30*(1/self.ncalon)) # started after 30 stamps
                #TONoff=TONon+1
                #TOFFcalonst=int(np.round(TONcalon+TONcalonst+(30/exp)*(1/self.ncalon))) #int((len(Non)-(15/exp))/2) # TONon + TONsec= (cycle-30)/2 + i*cycle
                #TOFFcalonend=int(np.round(2*TONcalon+(30/exp)*(1/self.ncalon)))
                
                ###################### Separate ON and OFF position data for caloff#########################################
                #TON=int((len(Noff)-(30/exp)*(self.ncalon-1)/self.ncalon)) #subtracted extra 30 stamps as stabilizing time
                #TONst=int(30*(self.ncalon-1)/self.ncalon) # started after 30 stamps
                #TONoff=TONon+1
                #TOFFst=int(np.round(TON+TONst+(30/exp)*((self.ncalon-1)/self.ncalon))) #int((len(Non)-(15/exp))/2) # TONon + TONsec= (cycle-30)/2 + i*cycle
                #TOFFend=int(np.round(2*TON+(30/exp)*((self.ncalon-1)/self.ncalon)))
                 #TOFFoff=TOFFon+1
                #print(TONst,TON,TOFFst,TOFFend)
                
               # TaON=mTa[TONst:TON-1,:,:]
               # TaOFF=mTa[TOFFst:TOFFend,:,:]
                
               #print(np.array(freq.shape,SC[0].shape,calon[:,SC[0],:].shape)
           ######################### Calibrate ON position data#########################################
                s=cal.cal(beam,self.driftdate[0],len(Non),len(Noff),len(SC[0]),np.array(freq[SC[0]]),calon,caloff,self.tcalpath)
               # s=cal.cal(beam,self.driftdate[0],TONcalon-TONcalonst,TON-TONst,len(SC[0]),np.array(freq[SC[0]]),calon[TONcalonst:TONcalon-1,:,:],caloff[TONst:TON-1,:,:],self.tcalpath)
                Ta, tcalt=s.tonoff() 
                
            ######################### Calibrate OFF position data#########################################    
               # s=cal.cal(beam,self.driftdate[0],TOFFcalonend-TOFFcalonst,TOFFend-TOFFst,len(SC[0]),np.array(freq[SC[0]]),calon[TOFFcalonst:TOFFcalonend,:,:],caloff[TOFFst:TOFFend,:,:],self.tcalpath)
               # TaOFF, tcalt=s.tonoff() 
                
               # print(TaON.shape,TaOFF.shape)
               # Ta=np.concatenate((TaON,TaOFF),axis=0)
                plt.plot(freq[SC[0]],tcalt[:,0:1],label=beam+'_'+str(i))
                plt.legend()
                plt.xlabel('Frequency (MHz)')
                plt.ylabel('Tcal (K)')
                plt.show()
               
               
               
                flag=self.flag
                if(flag=='ON'):
                   print("RFI removal after calibration for beam:",beam)
                  # mask,mTa=rfi.flagbadchan2(beam,Ta,len(SC[0]),np.array(freq[SC[0]]))
                   if(beam in flags['badbeam']):
                       if(i in flags['badscan']):
                           a=np.where(flags['badbeam']==beam)
                           b=np.where(flags['badscan'][a[0]]==i)
                           badtime=flags['badtime'][b[0]]
                           badtime=str(badtime[0]).split(',')
                           badchan=flags['badchan'][b[0]]
                           badchan=str(badchan[0]).split(',')
                           badpol=flags['badpol'][b[0]]
                           badpol=str(badpol[0]).split(',')
                           print('Flagged',badtime,badchan,badpol)
                           mask, mTa=rfi.flagbaddata(Ta,badtime,badchan,badpol)
                               
                elif(flag=='OFF'):
                   mask= Ta < -100000000000000000000000000
                   mask[:,:,:,]=False
                   mTa=np.ma.MaskedArray(Ta,mask=mask,fill_value=np.nan)
              
               ############################after calibration plots #############################################
               
              # GausssmoothF=gaussian_filter1d(Dav,300,axis=0)
              # a=Dav-GausssmoothF
                fig2, (ax1, ax2) = plt.subplots(figsize=(10, 7), ncols=2)

                AA = ax1.imshow(mTa[:,:,0], cmap='rainbow', aspect='auto', interpolation='none',origin='lower')
   #AA = ax3.imshow(flaggedata, cmap='seismic', aspect='auto', interpolation='none',origin='lower')
                fig2.colorbar(AA, ax=ax1)
                ax1.set_ylabel('Time stamp')
                ax1.set_xlabel('Channel')
                string=beam+'_'+str(i)+' Ta pol0'
                ax1.set_title(string)

                AA = ax2.imshow(mTa[:,:,1], cmap='rainbow', aspect='auto', interpolation='none',origin='lower')
   #AA = ax4.imshow(flaggedata2, cmap='seismic', aspect='auto', interpolation='none',origin='lower')
                fig2.colorbar(AA, ax=ax2)
                ax2.set_ylabel('Time stamp')
                ax2.set_xlabel('Channel')
                string=beam+'_'+str(i)+' Ta pol1'
                ax2.set_title(string)

             #  plt.savefig(prefix+'colorscale.png')
                plt.show()
                if(i==Ncycle-1):
                    TON=int((len(Noff)-(slewtime*((self.ncalon-1)/self.ncalon)/exp))/2) #30s switching time subt.
                else:
                    TON=int((len(Noff)-(2*slewtime*((self.ncalon-1)/self.ncalon)/exp))/2) #((self.ncalon-1)/self.ncalon)) #subtracted extra 30 stamps as stabilizing time
                TONst=int((slewtime/2)*((self.ncalon-1)/self.ncalon)/exp)#*(self.ncalon-1)/self.ncalon) # started after 30 stamps
                #TONoff=TONon+1
                TOFFst=int(np.round(TON+TONst+(slewtime*((self.ncalon-1)/self.ncalon))/exp)) #*((self.ncalon-1)/self.ncalon))) #int((len(Non)-(15/exp))/2) # TONon + TONsec= (cycle-30)/2 + i*cycle
                TOFFend=int(np.round(2*TON+(slewtime*((self.ncalon-1)/self.ncalon))/exp)) #/((self.ncalon-1)/self.ncalon)))
                 #TOFFoff=TOFFon+1
                print("TONst:",TONst,"TON:",TON,"TOFFst:",TOFFst,"TOFFend:",TOFFend)
                
                TaON=mTa[TONst:TON-1,:,:]
                TaOFF=mTa[TOFFst:TOFFend,:,:]
                TaONave=np.nanmedian(TaON,axis=0)
                TaOFFave=np.nanmedian(TaOFF,axis=0)
                if((beam=='M08') or (beam=='M14')):
                    Taave=TaOFFave-TaONave
                else:
                    Taave=TaONave-TaOFFave
               # print(TaON)
               # print(TaOFF)
                fig, (ax1, ax2, ax3) = plt.subplots(figsize=(10, 5), nrows=3)
                ax1.plot(freq[SC[0]],TaONave[:, 0], color = 'green', label = 'polarization_1 ON')
                
                #ax1.plot(freq[SC[0]],np.nanmean(TaONave[:, 0:1],axis=1), label = 'average in polarization')

                ax1.plot(freq[SC[0]],TaOFFave[:, 0], color = 'red', label = 'polarization_1 OFF')
                
                ax2.plot(freq[SC[0]],TaONave[:, 1], color = 'green', label = 'polarization_2 ON')
                ax2.plot(freq[SC[0]],TaOFFave[:, 1], color = 'red', label = 'polarization_2 OFF')
                #ax2.plot(freq[SC[0]],np.nanmean(TaOFFave[:, 0:1],axis=1), label = 'average in polarization')

                ax3.plot(freq[SC[0]],Taave[:, 0], color = 'red', label = 'polarization_1 ON-OFF')
                ax3.plot(freq[SC[0]],Taave[:, 1], color = 'green', label = 'polarization_2 ON-OFF')
                #ax3.plot(freq[SC[0]],np.nanmean(Taave[:, 0:1]), label = 'average in polarization')
                ax1.tick_params(top = 'True', right = 'True', which = 'both')
                ax2.tick_params(top = 'True', right = 'True', which = 'both')
                ax1.legend(loc = 'upper left', frameon = False)
                ax2.legend(loc = 'upper left', frameon = False)
                ax3.legend(loc = 'upper left', frameon = False)
                plt.show()
                
                Tout=Taave
                
                #######################################################################
                #if(i==0):
                #    dataout=np.array(Tout)
                #    mjd=MJD
                #    Tstamp=len(Tout[:,0,0])
                #    mask1=mask
                #else: 
                dataout.append(Tout)
                dataoffout.append(TaOFFave)
                dataonout.append(TaONave)
                #    dataout=np.concatenate((dataout,np.array(Tout)))
                #    mjd=np.concatenate((mjd,np.array(MJD)),axis=0)
                #    Tstamp=Tstamp+len(Tout[:,0,0]) 
                #    mask1=np.concatenate((mask1,mask),axis=0)
        print("Data shape final:",np.array(dataout).shape)
        dataoutfinal=np.nanmean(np.array(dataout),axis=0)
        dataoffoutfinal=np.nanmean(np.array(dataoffout),axis=0)
        dataonoutfinal=np.nanmean(np.array(dataonout),axis=0)
        
        fig,ax=plt.subplots(figsize=(10,5),nrows=1)    
        ax.plot(freq[SC[0]],dataoffoutfinal[:,0],label='OFF pol 0')
        ax.plot(freq[SC[0]],dataonoutfinal[:,0],label='ON pol 0')
        ax.set_xlabel('Frequency')
        ax.set_ylabel('Ta(K)')
        plt.legend()
        plt.show()  
        
        fig,ax=plt.subplots(figsize=(10,5),nrows=1)    
        ax.plot(freq[SC[0]],dataoffoutfinal[:,1],label='OFF pol 1')
        ax.plot(freq[SC[0]],dataonoutfinal[:,1],label='ON pol 1')
        ax.set_xlabel('Frequency')
        ax.set_ylabel('Ta(K)')
        plt.legend()
        plt.show()  
            
        fig,ax=plt.subplots(figsize=(10,5),nrows=1)    
        ax.plot(freq[SC[0]],dataoutfinal[:,0],label='ON-OFF pol 0')
        ax.plot(freq[SC[0]],gaussian_filter1d(dataoutfinal[:,0],5),color='red')
        ax.plot(freq[SC[0]],dataoutfinal[:,1],label='ON-OFF pol 1')
        ax.plot(freq[SC[0]],gaussian_filter1d(dataoutfinal[:,1],5),color='green')
        ax.set_xlabel('Frequency')
        ax.set_ylabel('Ta(K)')
        plt.legend()
        plt.show()
        
        
        #Mask the line and RFI frequencies if any
        #######################################################
        if(self.contsub=='ON'):
            mask= dataoutfinal < -100000000000000000000000000
            mask[:,:]=False
            #mdfout=np.ma.MaskedArray(datafinalout,mask=mask,fill_value=np.nan)
            i=0
            for el in maskfreq:
                a=np.min(np.where(np.array(freq[SC[0]]) > maskfreq[i][0]))
                b=np.max(np.where(np.array(freq[SC[0]]) < maskfreq[i][1]))
                print(a,b)
                mask[a:b,:]=True
                i=i+1
            x=np.ma.MaskedArray(freq[SC[0]],mask=mask[:,0],fill_value=np.nan)
            y=np.ma.MaskedArray(dataoutfinal,mask=mask,fill_value=np.nan)
        ############################################################################
            xData=x.filled()
            yData=y.filled() 
           # print(np.shape(xData),np.shape(yData),np.shape(mask),np.shape(y[:,0]))
        ############################################################################
           # print(np.shape(freq[SC[0]]),np.shape(dataoutfinal[:,0]))
            x0,y0,poly0,fun0,res0=cnsub.polyfitsub(np.array(xData),np.array(yData[:,0]))
            x1,y1,poly1,fun1,res1=cnsub.polyfitsub(np.array(xData),np.array(yData[:,1]))
        ###########################################################################   
            fig,ax=plt.subplots(figsize=(10,5),nrows=1)    
            ax.plot(x0,y0,label='ON-OFF pol 0')
       # ax.plot(freq[SC[0]],gaussian_filter1d(y[:,0],5),color='red')
            ax.plot(freq[SC[0]],poly0(freq[SC[0]]),color='black',lw=2)
            ax.plot(x1,y1,label='ON-OFF pol 1')
       # ax.plot(freq[SC[0]],gaussian_filter1d(y[:,1],5),color='green')
            ax.plot(freq[SC[0]],poly1(freq[SC[0]]),color='blue',lw=2)
            ax.set_xlabel('Frequency')
            ax.set_ylabel('Ta (K)')
            plt.legend()
            plt.show()
            
        ########################Remove left over ripples if any#####################################################      
            a=np.where(np.array(x0)==np.array(x0))
            fVals,amp0,Ph0=cnsub.FFT(a,x0,res0,chw)
            fVals,amp1,Ph1=cnsub.FFT(a,x1,res1,chw)
            
            #print(chw,fVals,amp1,Ph1)
            #########################################################################
            xlim=10

            fig, (ax1,ax2) = plt.subplots(nrows=2, ncols=1) #create figure handle
           # fVals=np.arange(start = -NFFT/2,stop = NFFT/2)*fs/NFFT
            ax1.plot(fVals,amp0,'b')
            ax1.plot(fVals,amp1,'b')
            ax1.set_title('Double Sided FFT - with FFTShift')
            ax1.set_xlabel('Frequency (MHz)')         
            ax1.set_ylabel('Amplitude')
            ax1.set_xlim(-xlim,xlim)
                
            ax2.plot(fVals,Ph0)
            ax2.plot(fVals,Ph1)
            ax2.set_xlim(-xlim,xlim)
            ax2.set_xlabel('Frequency (MHz)')         
            ax2.set_ylabel('Phase')
            #ax.set_xticks(np.arange(-50, 50+10,10))
            fig.show()
            
           
        #######################################################################################
            p=np.where(fVals>0)
            M0=np.max(amp0[p])
            m0=np.max(np.where(amp0==M0))
            #fp=1/(tarr[m]*1000000
            fp0=fVals[m0]#*np.power(10,6)
            C0=Ph0[m0]
            
            q=np.where(fVals>0)
            M1=np.max(amp1[q])
            m1=np.max(np.where(amp1==M1))
            #fp=1/(tarr[m]*1000000
            fp1=fVals[m1]#*np.power(10,6)
            C1=Ph0[m1]
            ##########################################################################
            
            
            #define sine function with frequency fp
            def f2(x,A,C,fp):
                    return A*np.sin(2*np.pi*fp*x+C)        
            
           # chw=1
            ini_val0=[M0,C0,fp0]
            param0, param_cov = curve_fit(f2,np.array(x0),np.array(gaussian_filter1d(res0,5)),p0=ini_val0)
            ini_val1=[M1,C1,fp1]
            param1, param_cov = curve_fit(f2,np.array(x1),np.array(gaussian_filter1d(res1,5)),p0=ini_val1)
            print(np.where(amp0==M0),np.where(amp1==M1))
            print(m0,m1)
            print(ini_val0,ini_val1)
            print(param0,param1)
            fig,(ax1,ax2)=plt.subplots(2)
            ax1.plot(x0,y0,x0,fun0)
            ax1.plot(x1,y1,x1,fun1)
            ax2.plot(x0,gaussian_filter1d(res0,5))
            ax2.plot(x1,gaussian_filter1d(res1,5)) 
            ax2.plot(x0,f2(x0,param0[0],param0[1],param0[2]))
            ax2.plot(x0,f2(x0,param1[0],param1[1],param1[2]))
           # ax3.plot(x0,res0,'D')
           # ax3.plot(x1,res1,'D') 
           # ax3.plot(x0,f2(x0,param0[0],param0[1],param0[2]))
           # ax3.plot(x1,f2(x1,param1[0],param1[1],param1[2]))  
            fig.show()
            
            
            
        #############################################################################    
            dataoutfinal[:,0]=dataoutfinal[:,0]-poly0(freq[SC[0]])-f2(freq[SC[0]],param0[0],param0[1],param0[2])
            dataoutfinal[:,1]=dataoutfinal[:,1]-poly1(freq[SC[0]])-f2(freq[SC[0]],param1[0],param1[1],param1[2])
           # Tout=cnsub.continuum_sub2(dataoutfinal,1,np.array(freq[SC[0]]),lffreq[0],lffreq[1],lffreq[2],lffreq[3])
##############################################################################
        def G(ZEN):
                if(ZEN < 26.4):
                    if((sfreq >= 1050) & (sfreq < 1100)):            # Jiang et al. (2020)
                        out=gf1[self.b]  
                    elif((sfreq >= 1100) & (sfreq < 1150)):
                        out=gf2[self.b]
                    elif((sfreq >= 1150) & (sfreq < 1200)):
                        out=gf3[self.b]
                    elif((sfreq >= 1200) & (sfreq < 1250)):
                        out=gf4[self.b]
                    elif((sfreq >= 1250) & (sfreq < 1300)):
                        out=gf5[self.b]
                    elif((sfreq >= 1300) & (sfreq < 1350)):
                        out=gf6[self.b]
                    elif((sfreq >= 1350) & (sfreq < 1400)):
                        out=gf7[self.b]
                    elif((sfreq >= 1400) & (sfreq < 1450)):
                        out=gf8[self.b]
                    elif((sfreq >= 1450)):
                        out=gf9[self.b]
    
                if(ZEN> 26.4):
        #a=np.arange(26.4,63.4,1.0)
                    out=16.46-0.02*np.power(ZEN-26.4,2)-0.12*(ZEN-26.4) # Zhang et al. 2019
                return out   
 #########################  conversion to mJy########################  
     
        flux=1000*dataoutfinal/G(ZEN)
        fig,ax=plt.subplots(figsize=(10,5),nrows=1)    
        ax.plot(freq[SC[0]],flux[:,0],label='ON-OFF pol 0')
        ax.plot(freq[SC[0]],gaussian_filter1d(flux[:,0],5),color='red')
        ax.plot(freq[SC[0]],flux[:,1],label='ON-OFF pol 1')
        ax.plot(freq[SC[0]],gaussian_filter1d(flux[:,1],5),color='green')
        ax.set_xlabel('Frequency')
        ax.set_ylabel('Flux (mJy)')
        plt.legend()
        plt.show()
        
        fig,ax=plt.subplots(figsize=(10,5),nrows=1)    
        ax.plot(freq[SC[0]],np.average(flux[:,0:1],axis=1),label='ON-OFF pol av')
        ax.plot(freq[SC[0]],gaussian_filter1d(np.average(flux[:,0:1],axis=1),5),color='red')
        ax.set_xlabel('Frequency [MHz]')
        ax.set_ylabel('Flux [mJy]')
        plt.legend()
        plt.show()
        
#################################################################################################
        F=freq[SC[0]]
        C=np.arange(len(SC[0]))+1 
        Flux=[]
        V=fqv.freqtovel(F,'optical',1420.405752)
        lon,lat,time=fqv.obsloc(np.nanmean(X),np.nanmean(Y),np.nanmean(Z),np.nanmean(MJD1))
        print('lon:',lon.deg,'lat:',lat.deg,'time:',time)
        t=Time(time,format='mjd',scale='utc')
        hV=dps.heliocentric_velocity_correction(V,lon.deg,lat.deg, t)
        r=0
        for r in range(len(R)):
            Flux.append(flux)
            r=r+1
            
        fig,ax=plt.subplots(figsize=(10,5),nrows=1)    
        ax.plot(hV,np.average(flux[:,0:1],axis=1),label='ON-OFF pol av')
        ax.plot(hV,gaussian_filter1d(np.average(flux[:,0:1],axis=1),5),color='red')
        ax.set_xlabel('Velocity-Heliocentric [km/s]')
        ax.set_ylabel('Flux [mJy]')
        plt.legend()
        plt.show()
       # print("Co-ordinate calculation done for beam:",beam,"\n")
       # return R,D,EL,AZ,RAN,Ton,Toff,Ta,Fluxon,Fluxoff,Flux,F,C,TON,TOFF, MJD, freq1,chw
        return R,D,EL,AZ,Flux,F,hV,C,T,MJD,freq1,chw
############################################################################
#-(15/(60*60*24)
#-14/(60*60*24)



