# class and functions to deal with rays (for ColumnDensity code).
from scipy.io import FortranFile 
import numpy as np
import lya_utils as lya  # all lya-specific constants and conversions.

class raylist(object):
    
    def __init__(self,icFile,resFile,load=True):
        # init object from results of a run 
        self.icFile  = icFile
        self.resFile = resFile
        if load:
            self.load_ic()
            self.load_res()

    def load_ic(self):
        # read rays IC file
        f = FortranFile(self.icFile)
        [self.nrays]  = f.read_ints()
        self.ID_ic = f.read_ints()
        self.nu_ic = f.read_reals()
        xx = f.read_reals()
        xx = xx.reshape((self.nrays,3))
        self.x_ic = xx[:,0]
        self.y_ic = xx[:,1]
        self.z_ic = xx[:,2]
        xx = f.read_reals()
        xx = xx.reshape((self.nrays,3))
        self.kx_ic = xx[:,0]
        self.ky_ic = xx[:,1]
        self.kz_ic = xx[:,2]
        f.close()

            
    def load_res(self): 
        # read results from RASCAS
        f = FortranFile(self.resFile)
        [self.nrays] = f.read_ints()
        # ID
        self.ID = f.read_ints() 
        # dist
        self.dist = f.read_reals()
        # tau
        self.dist = f.read_reals()
        f.close()


        
