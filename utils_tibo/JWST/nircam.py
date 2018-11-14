
from astropy.io import ascii
import numpy as np
import os

class nircamFilter(object):
    
    def __init__(self,filterName,throughputModel=None,nircamDir=None):
        self.filterName = filterName
        if throughputModel is None:
            self.throughputModel = '_NRC_and_OTE_ModAB_mean'
        else:
            self.throughputModel = throughputModel
        if nircamDir is None:
            self.nircamDir = 'nircam_throughputs/modAB_mean/nrc_plus_ote/'
        else:
            self.nircamDir = nircamDir
        self.filterFile = "%s%s%s%s"%(self.nircamDir, self.filterName, self.throughputModel, '.txt')
        if not os.path.exists(self.filterFile):
            print('filterFile does not exist...',self.filterFile)
        else:
            self.load()
            
    def load(self,threshold=None):

        filterTab = ascii.read(self.filterFile)
        wl = filterTab['microns'].data
        throughput = filterTab['throughput'].data
        if threshold is None:
            threshold = 1.e-4
        good = np.where(throughput > threshold)
        self.wl = wl[good] * 1.e4 # convert to Angstrom
        self.throughput = throughput[good]
        self.lambda_min = np.amin(self.wl)
        self.lambda_max = np.amax(self.wl)
        self.lambda_mean = np.mean(self.wl)
        
def plotFilter(filter):
    import matplotlib.pyplot as plt
    f,a = plt.subplots()
    a.semilogy()
    a.plot(filter.wl,filter.throughput)
