
import numpy as np
from scipy.io import FortranFile as ff

class header(object):
    def __init__(self,file):
        f = open(file,'r')
        line         = f.readline()   # skip this line
        self.dir     = f.readline().split()[3]
        self.output  = f.readline().split()[3]
        line         = f.readline()   # skip this line
        self.boxlen  = float(f.readline().split()[2])
        self.time    = float(f.readline().split()[2])
        self.aexp    = float(f.readline().split()[2])
        self.H0      = float(f.readline().split()[2])
        self.omega_m = float(f.readline().split()[2])
        self.omega_l = float(f.readline().split()[2])
        self.omega_k = float(f.readline().split()[2])
        self.omega_b = float(f.readline().split()[2])
        self.unit_l  = float(f.readline().split()[2])
        self.unit_d  = float(f.readline().split()[2])
        self.unit_t  = float(f.readline().split()[2])
        line         = f.readline()    # skip this line
        f.close()
        print '   aexp & z                  =',self.aexp,(1./self.aexp)-1.
        print '   scale_l, scale_d, scale_t =',self.unit_l,self.unit_d,self.unit_t
        print '   simulation time [Myr]     =',self.time*self.unit_t/(365.*24.*3600.*1e6)
    

class stars(object):

    def __init__(self,file,load=True):
        self.file = file
        if load:
            self.load()
        self.units = 'pos in box units, vel in cm/s, mass in g, age in Myr, Z in absolute value'

    def load(self):
        print "   Reading Subvol file: ", self.file
        f = ff(self.file)
        [nstar] = f.read_ints()
        print '   nstar =',nstar
        xx = f.read_reals('d')
        #xx = xx.reshape((nstar,3),order='F')
        xx = xx.reshape((nstar,3))
        self.x = xx[:nstar,0]
        self.y = xx[:nstar,1]
        self.z = xx[:nstar,2]
        xx = f.read_reals('d')
        xx = xx.reshape((nstar,3))
        self.vx = xx[:nstar,0]
        self.vy = xx[:nstar,1]
        self.vz = xx[:nstar,2]
        self.mp = f.read_reals('d')
        #ap = f.read_reals('d')
        self.tp = f.read_reals('d')
        self.zp = f.read_reals('d')
        #idp = f.read_ints()
        f.close()

    def convert(self, header):

        # convert units
        # from g to solar mass
        #amass = scale_d * scale_l**3
        Msol  = 1.9891e+33
        #amass = amass/Msol
        self.mp /= Msol

        # from box unit to kpc
        size = header.boxlen*header.unit_l/3.08e21
        # print '   subvol size (kpc) ',r0,r0*size
        self.x *= size
        self.y *= size
        self.z *= size

        # from cm/s to km/s
        scale_v = 1./1.e5
        self.vx *= scale_v
        self.vy *= scale_v
        self.vz *= scale_v

        self.units = 'pos in kpc, vel in km/s, mass in Msun, age in Myr, Z in absolute value'
        print "   Units: pos in kpc, vel in km/s, masses in Msun"    

