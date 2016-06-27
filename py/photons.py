
import fortranfile
import numpy as np

clight   = 2.9979250e+10       #[cm/s] light speed
cmtoA    = 1.e8                # from cm to A
lambda_0 = 1215.6701           #[A] Lya wavelength
lambda_0_cm = lambda_0 / cmtoA              # cm
nu_0 = clight / lambda_0_cm                 # Hz

class photons(object):

    def __init__(self,ID,status,x,y,z,nu,wl,nscat,time):
        self.ID = ID
        self.status =status
        self.x = x
        self.y = y
        self.z = z
        self.nu = nu
        self.wl = wl
        self.nscat = nscat
        self.time = time

def from_file(file):

    print "--> reading photons in file ",file
    f = fortranfile.FortranFile(file)
    
    [nphoton] = f.readInts()
    print "--> Nphotons =",nphoton
    
    # ID
    xx = f.readInts()
    ID = xx 
    # status
    xx = f.readInts()
    status = xx
    # xlast(3)
    xx = f.readReals('d')
    print np.shape(xx)
    xx = xx.reshape((nphoton,3))
    print np.shape(xx)
    x = xx[:,0]
    y = xx[:,1]
    z = xx[:,2]
    # nu_ext
    xx = f.readReals('d')
    nu = xx
    # nb_abs
    xx = f.readInts()
    nscat = xx
    # time
    xx = f.readReals('d')
    time = xx
    
    f.close()

    wl = clight / nu * cmtoA
    
    print "--> read done"

    classP = photons(ID,status,x,y,z,nu,wl,nscat,time)

    return classP
