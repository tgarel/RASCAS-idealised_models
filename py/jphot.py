# nouvelle tentative de classe photons.
import fortranfile
import numpy as np
import lya_utils as lya  # all lya-specific constants and conversions.

class photonlist(object):
    
    def __init__(self,icFile,resFile,load=True,stars=False):
        self.icFile  = icFile
        self.resFile = resFile
        if load: 
            self.load_ic(stars=stars)
            self.load_res()

            
    def load_ic(self,stars=False):
        # read photn IC file
        f = fortranfile.FortranFile(self.icFile)
        [self.nphoton]  = f.readInts()
        [self.nRealPhotons] = f.readReals('d')
        [self.iseed_ic] = f.readInts()
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
        if stars:
            self.nu_star = f.readReals('d')
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
        # k(3)
        xx = f.readReals('d')
        xx = xx.reshape((self.nphoton,3))
        self.kx = xx[:,0]
        self.ky = xx[:,1]
        self.kz = xx[:,2]
        # nb_abs
        xx = f.readInts()
        self.nscat = xx
        # time
        xx = f.readReals('d')
        self.time = xx
        f.close()


    def extract_sample(self,indexes):
        # return a photon object, which is a subsample of self defined by indexes (i.e. self[indexes])
        p = photonlist(self.icFile,self.resFile,load=False)
        p.nphoton = len(indexes)
        # simple copy of numbers 
        p.iseed_ic = self.iseed_ic
        # selection in arrays
        p.ID_ic    = self.ID_ic[indexes]
        p.nu_ic    = self.nu_ic[indexes]
        p.x_ic     = self.x_ic[indexes]        
        p.y_ic     = self.y_ic[indexes]
        p.z_ic     = self.z_ic[indexes]
        p.kx_ic    = self.kx_ic[indexes]        
        p.ky_ic    = self.ky_ic[indexes]
        p.kz_ic    = self.kz_ic[indexes]
        p.iran_ic  = self.iran_ic[indexes]
        p.ID       = self.ID[indexes]
        p.status   = self.status[indexes]
        p.x        = self.x[indexes]        
        p.y        = self.y[indexes]
        p.z        = self.z[indexes]
        p.nu       = self.nu[indexes]
        p.nscat    = self.nscat[indexes]
        p.time     = self.time[indexes]
        p.kx       = self.kx[indexes]        
        p.ky       = self.ky[indexes]
        p.kz       = self.kz[indexes]
        return p      


    def append(self,p): 
        # append a photonlist p to self.
        self.icFile = 'mix'
        self.resFile = 'mix'
        self.nphoton = self.nphoton + p.nphoton
        # reset seed to absurd value 
        self.iseed_ic = 0
        # selection in arrays
        self.ID_ic    = np.append(self.ID_ic,p.ID_ic)
        self.nu_ic    = np.append(self.nu_ic,p.nu_ic)
        self.x_ic     = np.append(self.x_ic,p.x_ic)
        self.y_ic     = np.append(self.y_ic,p.y_ic)
        self.z_ic     = np.append(self.z_ic,p.z_ic)
        self.kx_ic    = np.append(self.kx_ic,p.kx_ic)
        self.ky_ic    = np.append(self.ky_ic,p.ky_ic)
        self.kz_ic    = np.append(self.kz_ic, p.kz_ic)
        self.iran_ic  = np.append(self.iran_ic,p.iran_ic)
        self.ID       = np.append(self.ID,p.ID)
        self.status   = np.append(self.status,p.status)
        self.x        = np.append(self.x,p.x)
        self.y        = np.append(self.y,p.y)
        self.z        = np.append(self.z,p.z)
        self.nu       = np.append(self.nu,p.nu)
        self.nscat    = np.append(self.nscat,p.nscat)
        self.time     = np.append(self.time,p.nscat)
        self.kx       = np.append(self.kx,p.kx)
        self.ky       = np.append(self.ky,p.ky)
        self.kz       = np.append(self.kz,p.kz)
        return p      

    
            
    def project_pos(self,k,thetamax):
        # compute projected positions (x,y) of photons going along k, in a plane perp. to k.
        # NB: used with thetamax=180, this selects all photons and may be used to compute
        #     a simple projection along any direction
        # inputs :
        #    - k: a 3D array [kx,ky,kz]
        #    - thetamax : max angle away from k, in deg.
        costheta = self.kx * k[0] + self.ky * k[1] + self.kz * k[2]
        # select photons along k +/- some
        ii = np.where(costheta > np.cos(thetamax*np.pi/180.)) 
        x  = np.empty_like(ii)
        y  = np.empty_like(ii)
        # define projection basis
        if (k[0] < 1.):  # then we k not colinear with x
            #-> compute u1 as cross prod of k and x. 
            u1_x = 0.
            u1_y = k[2]
            u1_z = -k[1]
            # compute u2 as cross product of k and u1.
            u2_x = k[1]*u1_z - k[2]*u1_y
            u2_y = -k[0]*u1_z
            u2_z = k[0]*u1_y
            # normalize u1 and u2
            mod_u1 = np.sqrt(u1_y*u1_y+u1_z*u1_z)
            u1_y = u1_y / mod_u1
            u1_z = u1_z / mod_u1
            mod_u2 = np.sqrt(u2_x*u2_x+u2_y*u2_y+u2_z*u2_z)
            u2_x = u2_x / mod_u2
            u2_y = u2_y / mod_u2
            u2_z = u2_z / mod_u2
        else:   # k is along x -> use u1 = y, u2 = z
            u1_x = 0.
            u1_y = 1.
            u1_z = 0.
            u2_x = 0.
            u2_y = 0.
            u2_z = 1.
        # compute projected coordinates
        x = self.x[ii] * u1_x + self.y[ii] * u1_y + self.z[ii] * u1_z
        y = self.x[ii] * u2_x + self.y[ii] * u2_y + self.z[ii] * u2_z
        return x,y
