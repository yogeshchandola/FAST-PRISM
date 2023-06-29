from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, AutoMinorLocator,IndexLocator, FormatStrFormatter
import numpy as np
import astropy
from astropy.coordinates import SkyCoord, Angle, EarthLocation
from astropy.time import Time
from astropy import units as u
from astropy.utils.data import clear_download_cache
import aoflagger 

class rfiremove:
    def __init__(self,nstamp,nchan,values,npol):
        self.nstamp=nstamp
        self.nchan=nchan
        self.npol=npol
        self.values=values
    def flagger(self):
        aof=aoflagger.AOFlagger()
        #strategy = aof.load_strategy("/data/FAST_cluster_project_Yinzhe/pipeline3.1/cube/process/modules/default-2.14.0.rfis")
        strategy = aof.load_strategy("/data/FAST_cluster_project_Yinzhe/pipeline3.1/cube/process/modules/newstrategy7.rfis")
        #strategy = aof.load_strategy("/media/yogesh/5TBAIPSDATA1/FAST_cluster_project_Yinzhe/pipeline3.1/cube/process/modules/J1550flagnewstrategy.rfis")
#nstamp=10000
#nchan=256
        data = aof.make_image_set(int(self.nstamp), int(self.nchan), int(self.npol))
        print("Number of times: " + str(data.width()))
        print("Number of channels: " + str(data.height()))
        print(np.shape(data))
# When flagging multiple baselines, iterate over the baselines and
# call the following code for each baseline
# (to use multithreading, make sure to create an imageset for each
# thread)
# Make four images:  for 4 pol
        flagarray=[]
        i=0
        #imgindex=0
        for imgindex in range(int(self.npol)):
        # Initialize data
        #    values = np.zeros([int(self.nchan), int(self.nstamp)])
        #    print(np.shape(self.values))
        #    a=np.array(self.values)
            print(np.shape(self.values))
        #   a=np.array(self.values[i,:,:],dtype='double')
            a=np.array(self.values[:,:],dtype='double')
            print(np.shape(a))
            print(np.array(a).dtype)
            data.set_image_buffer(imgindex, a)
            flags = aof.run(strategy, data)
            flagvalues = flags.get_buffer()
            flagcount = sum(sum(flagvalues))
            print(flagvalues,flagcount,int(self.nchan)*int(self.nstamp))
            print("Percentage flags on data: " + str(flagcount * 100.0 / (int(self.nchan)*int(self.nstamp))) + "%")
            flagarray.append(flagvalues)
            i=i+1
        
        print(np.shape(flagarray))
       # flagarray=np.array(flagarray)
        #data=davpON
       # print(np.shape(data))
       # flags = aof.run(strategy, data)
       # flagvalues = flags.get_buffer()
       # print(np.shape(flagvalues))
        return flagarray
