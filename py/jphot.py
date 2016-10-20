# nouvelle tentative de classe photons.
import fortranfile
import lya_utils as lya  # all lya-specific constants and conversions.

class photonlist(object):
    
    def __init__(self,icFile,resFile):
        self.icFile  = icFile
        self.resFile = resFile

        
    def load_ic(self):
        # read photn IC file
        f = fortranfile.FortranFile(self.icFile)
        [self.nphoton] = f.readInts()
        [self.iseed_ic]   = f.readInts()
        self.ID_ic = f.readInts()
        self.nu_ic = f.readReals('d')
        xx = f.readReals('d')
        xx = xx.reshape((self.nphoton,3))
        self.x_ic = xx[:,0]
        self.y_ic = xx[:,1]
        self.z_ic = xx[:,2]
        xx = f.readReals('d')
        xx = xx.reshape((self.nphoton,3))
        self.kx_ic = xx[:,0]
        self.ky_ic = xx[:,1]
        self.kz_ic = xx[:,2]
        self.iran_ic = f.readInts()
        f.close()

                
    def load_res(self): 
        # read results from RASCAS
        f = fortranfile.FortranFile(self.resFile)
        [self.nphoton] = f.readInts()
        # ID
        xx = f.readInts()
        self.ID = xx 
        # status
        xx = f.readInts()
        self.status = xx
        # xlast(3)
        xx = f.readReals('d')
        xx = xx.reshape((self.nphoton,3))
        self.x = xx[:,0]
        self.y = xx[:,1]
        self.z = xx[:,2]
        # nu_ext
        xx = f.readReals('d')
        self.nu = xx
        # nb_abs
        xx = f.readInts()
        self.nscat = xx
        # time
        xx = f.readReals('d')
        self.time = xx
        f.close()

        

