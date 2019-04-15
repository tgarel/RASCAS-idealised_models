
# -*- coding: utf-8 -*-
'''
Modification line115 *unit_kpc to trace domain data and computation
'''

import numpy as np

class domain(object):
    
    def __init__(self, shape):
        self.shape = shape
        
    @classmethod
    def read(cls, filename):
        
        f = open(filename,'r')
        line   = f.readline()
        shape  = line.strip()
        if shape=='sphere':
            line   = f.readline()
            center = np.array(line.split(),dtype=float)
            line   = f.readline()
            radius = np.float(line)
            f.close()
            return sphere(center,radius)

        elif shape=='shell':
            line   = f.readline()
            center = np.array(line.split(),dtype=float)
            line   = f.readline()
            rin    = np.float(line)
            line   = f.readline()
            rout   = np.float(line)
            f.close()
            return shell(center,rin,rout)

        elif shape=='cube':
            line   = f.readline()
            center = np.array(line.split(),dtype=float)
            line   = f.readline()
            size   = np.float(line)
            f.close()
            return cube(center,size)

        else:
            return IOError('bouh')
        

class sphere(domain):
    
    def __init__(self, center, radius):
        domain.__init__(self, 'sphere')
        self.radius = radius
        self.center = center
        
    def area(self):
        return 4*np.pi*self.radius*self.radius

    def volume(self):
        return 4./3.*np.pi*(self.radius)**3.

    def info(self):
        print("domain INFO:")
        print("|_shape  =",self.shape)
        print("|_radius =",self.radius)
        print("|_center =",self.center)


class cube(domain):
    
    def __init__(self, center, size):
        domain.__init__(self, 'cube')
        self.center = center
        self.size   = size
        
    def area(self):
        return self.size*self.size

    def volume(self):
        return (self.size)**3.

    def info(self):
        print("domain INFO:")
        print("|_shape  =",self.shape)
        print("|_size   =",self.radius)
        print("|_center =",self.center)

class shell(domain):

    def __init__(self, center, rin, rout):
        domain.__init__(self, 'shell')
        self.center = center
        self.rin    = rin
        self.rout   = rout

    def area(self):
        return 4*np.pi*(self.rout**2-self.rin**2)

    def volume(self):
        return 4./3.*np.pi*(self.rout**3.-self.rin**3.)

    def info(self):
        print("domain INFO:")
        print("|_shape  =",self.shape)
        print("|_center =",self.center)
        print("|_rin    =",self.rin)
        print("|_rout   =",self.rout)

        
def overplot_limits(domain, unit_kpc, color=None, linestyle=None, linewidth=None):
    # OVERPLOT DOMAIN LIMITS
    from matplotlib.patches import Circle, Wedge, Polygon
    from matplotlib.collections import PatchCollection
    patches=[]
    if (domain.shape == 'sphere'):
#        wedge = Wedge((domain.center[0]*unit_kpc, domain.center[1]*unit_kpc), domain.radius*unit_kpc, 0., 360.)
	wedge = Wedge((0.0, 0.0), domain.radius*unit_kpc, 0., 360.)
    elif (meshDom.shape == 'shell'):
        wedge = Wedge((domain.center[0], domain.center[1]), domain.rout, 0., 360., width=domain.rout-domain.rin)
    patches.append(wedge)
    if color is None:
        color='k'
    if linestyle is None:
        linestyle='solid'
    if linewidth is None:
        linewidth=2
    p = PatchCollection(patches, edgecolors=color, facecolors='none', linewidths=linewidth, linestyles=[linestyle],zorder=1)
    return p
