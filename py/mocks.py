import lya_utils as ly
import jphot as jp
import numpy as np

class mockobs(object):

    #
    # This class handles mock observations output by RASCAS
    #

    def __init__(self,Dir,FileName,load_flux=False,load_spectrum=False,load_image=False,load_cube=False,unit_l_arcsec=1.0,lumdist_cm=1.0,redshift=0.0,cube_h5_file=None,idirection=1,ICFile=None):
        self.Dir      = Dir
        self.FileName = FileName 
        if ICFile is None:
            ICFile = FileName.replace('result.','')
            
        p = jp.photonlist("%s/%s"%(self.Dir,ICFile),'',load=False)
        nRealPhotons = p.get_nRealPhotons()
        nPhotons     = p.get_nphoton()
        self.nPhotPerPacket = nRealPhotons / nPhotons
        self.LumPerPacket = self.nPhotPerPacket * ly.h_cgs * ly.nu0 

        if load_flux:
            self.aperture,self.flux = self.__read_flux(index=idirection)
            #self.flux = self.flux * self.nPhotPerPacket
            self.flux = self.flux * self.LumPerPacket # [erg/s]
        if load_image:
            self.imsize_cu, self.image_npix, self.imtot = self.__read_image(index=idirection)
            p.load_ic()
            self.LumPerPacket = self.nPhotPerPacket * ly.h_cgs * np.mean(p.nu_ic)
            self.imtot = self.imtot * self.LumPerPacket # [erg/s]
        if load_cube:  # obs-frame cube
            self.cube_nlbda,self.cube_nxy,self.cube_lmin,self.cube_lmax,self.cube_imsize,self.cube = self.__read_cube(index=idirection)
            self.cube_lmin = self.cube_lmin * (1+redshift)
            self.cube_lmax = self.cube_lmax * (1+redshift)
            dl = (self.cube_lmax - self.cube_lmin) / self.cube_nlbda  # [A] 
            #l = np.arange(self.cube_lmin,self.cube_lmax,dl) + 0.5*dl  # [A]
            l = np.linspace(self.cube_lmin+ 0.5*dl, self.cube_lmax-0.5*dl, self.cube_nlbda)
            self.cube_lbda_Angstrom  = l 
            self.cube_dlbda_Angstrom = dl
            l = l * 1e-8  # [cm]
            energy = ly.h_cgs * ly.clight / l      # [erg]
            energy = energy * self.nPhotPerPacket  # [erg / s / phot packet]
            energy = energy / dl                   # [erg / s / A / phot packet]
            energy = energy / (self.cube_imsize*unit_l_arcsec/self.cube_nxy)**2 # [erg / s / A / arcsec2 / phot packet]
            energy = energy / (4. * np.pi * lumdist_cm**2)  # [erg / s/ A / arcsec2 / cm2 / phot packet] 
            self.cube = energy * self.cube   # [erg / s / A / arcsec2 / cm2 ]
            x = np.arange(-0.5*self.cube_imsize,0.5*self.cube_imsize,self.cube_imsize/self.cube_nxy) + 0.5 * self.cube_imsize/self.cube_nxy
            self.cube_x_arcsec = x * unit_l_arcsec
            self.cube_y_arcsec = x * unit_l_arcsec
            self.cube_xmin_arcsec = -0.5 * self.cube_imsize * unit_l_arcsec
            self.cube_xmax_arcsec = 0.5 * self.cube_imsize * unit_l_arcsec
            self.cube_ymin_arcsec = -0.5 * self.cube_imsize * unit_l_arcsec
            self.cube_ymax_arcsec = 0.5 * self.cube_imsize * unit_l_arcsec
            self.cube_dx_arcsec = self.cube_imsize * unit_l_arcsec /self.cube_nxy

    def __read_flux(self,index=1):
        from scipy.io import FortranFile as ff
        f = ff('%s/%s_flux.%5.5i'%(self.Dir,self.FileName,0),'r')
        for j in range(index):
            aperture,flux = f.read_reals('d')
        f.close()
        return aperture,flux

    def __read_image(self,index=1):
        from scipy.io import FortranFile as ff
        f = ff('%s/%s_image.%5.5i'%(self.Dir,self.FileName,0))
        for j in range(index):
            n        = f.read_ints()[0]
            imsize   = f.read_reals('d')[0]
            imcenter = f.read_reals('d')
            im       = f.read_reals('d')
        im = im.reshape((n,n))
        f.close()
        return imsize,n,im
    
    def show_image(self,plotFile,vmin=1e-15,vmax=10,showCB=True,smooth=False,smooth_scale_pix=2,noiseLevel=0):
        #import matplotlib
        #matplotlib.use('Agg')
        from matplotlib import pyplot as plt
        from matplotlib.colors import LogNorm

        image = self.imtot
        if smooth:
            import scipy.ndimage as ndimage
            image = ndimage.gaussian_filter(image,sigma=(smooth_scale_pix,smooth_scale_pix))

        if noiseLevel > 0: 
            image = image + np.random.normal(0.0,noiseLevel,image.shape)
        plt.figure(figsize=(10,10))
        plt.imshow(image,norm=LogNorm(),vmin=vmin,vmax=vmax)
        if showCB :
            plt.colorbar()
        plt.tight_layout()
        #plt.savefig(plotFile)

    def __read_cube(self,index=1):
        from scipy.io import FortranFile as ff
        f = ff('%s/%s_cube.%5.5i'%(self.Dir,self.FileName,0))
        for j in range(1,index+1):
            nlbda,nxy = f.read_ints()
            lmin,lmax,imsize = f.read_reals('d')
            center = f.read_reals('d')
            cube = f.read_reals('d')
        cube = cube.reshape((nxy,nxy,nlbda))
        f.close()
        return nlbda,nxy,lmin,lmax,imsize,cube

    def show_cube_image(self,ax=None,vmin=1e-21,vmax=1e-16,smooth=False,smooth_sigma_arcsec=0.285):
        from numpy import sum
        import matplotlib.pyplot as plt
        from matplotlib.colors import LogNorm
        from mpl_toolkits.axes_grid1.inset_locator import inset_axes

        #plt.ion()
        #plt.figure()
        #plt.clf()
        image = sum(self.cube,axis=2) * self.cube_dlbda_Angstrom
        if smooth:
            import scipy.ndimage as ndimage
            smooth_scale_pix = smooth_sigma_arcsec / self.cube_dx_arcsec
            image = ndimage.gaussian_filter(image,sigma=(smooth_scale_pix,smooth_scale_pix))
            cs = plt.contour(image,[1e-20,1e-19,1e-18],origin='lower',
                                 extent=(self.cube_xmin_arcsec,self.cube_xmax_arcsec,self.cube_xmin_arcsec,self.cube_xmax_arcsec),
                                 colors=['limegreen','deepskyblue','red'],linestyles=['--','-','-'],
                                 linewidths=1)
            
        if ax is None:
            ax=plt.subplot(1,1,1)
        im= plt.imshow(image,norm=LogNorm(),vmin=vmin,vmax=vmax,origin='lower', \
                    extent=[self.cube_xmin_arcsec,self.cube_xmax_arcsec,self.cube_xmin_arcsec,self.cube_xmax_arcsec],cmap='Greys')
        plt.xlabel('arcsec')
        plt.ylabel('arcsec')
        axins = inset_axes(ax,width="100%",height="3%",loc=1,borderpad=0)
        cbar = plt.colorbar(im,cax=axins,orientation='horizontal')
        if smooth:
            cbar.add_lines(cs)

    def show_cube_spec(self,smooth=False,smooth_sigma_angstrom=[0.5],unit='Angstrom'):
        from numpy import sum
        import matplotlib.pyplot as plt
        #plt.figure()
        if unit == 'Angstrom':
            lbdafact = 1.0
            lbdalab  = r'$[\AA]$'
        if unit == 'micron':
            lbdafact = 1e-4
            lbdalab  = r'$[\mu m]$'

        spectrum = sum(self.cube,axis=(0,1)) * self.cube_dx_arcsec * self.cube_dx_arcsec
        if smooth:
            import scipy.ndimage as ndimage
            for ssa in smooth_sigma_angstrom:
                smooth_scale_pix = ssa / self.cube_dlbda_Angstrom
                sspectrum = ndimage.gaussian_filter(spectrum,sigma=(smooth_scale_pix))
                plt.plot(self.cube_lbda_Angstrom*lbdafact,sspectrum,label=r'$\sigma_\lambda =$ %.1f $\AA$'%(ssa))

        plt.plot(self.cube_lbda_Angstrom*lbdafact,spectrum,color='gray',alpha=0.5,linewidth=1)
        if smooth:
            plt.legend()
        #plt.axvline(1215.67*4.20039,color='red',alpha=0.5)
        plt.xlabel(r'$\lambda_{\rm obs}$ %s'%lbdalab,fontsize=15)
        plt.ylabel(r'$F_\lambda \ [erg s^{-1} cm^{-2} \AA^{-1}]$',fontsize=15)


    def save_cube_h5(self,filename):
        import h5py 
        f = h5py.File(filename,'w')
        f.create_dataset('cube',data=self.cube)
        f.create_dataset('nlbda',data=self.cube_nlbda)
        f.create_dataset('nxy',data=self.cube_nxy)
        f.create_dataset('lbda_Angstrom',data=self.cube_lbda_Angstrom)
        f.create_dataset('dlbda_Angstrom',data=self.cube_dlbda_Angstrom)
        f.create_dataset('x_arcsec',data=self.cube_x_arcsec)
        f.create_dataset('y_arcsec',data=self.cube_y_arcsec)
        f.create_dataset('xmin_arcsec',data=self.cube_xmin_arcsec)
        f.create_dataset('xmax_arcsec',data=self.cube_xmax_arcsec)
        f.create_dataset('ymin_arcsec',data=self.cube_ymin_arcsec)
        f.create_dataset('ymax_arcsec',data=self.cube_ymax_arcsec)
        f.create_dataset('dx_arcsec',data=self.cube_dx_arcsec)
        f.close()



####
# Classes to help prepare mock observations
####

class specParams(object):
    def __init__(self,npix,aperture,lmin,lmax):
        self.npix     = npix
        self.aperture = aperture
        self.lmin     = lmin
        self.lmax     = lmax
    def write(self,f):
        f.write("%i %.10e %.10e %.10e \n "%(self.npix,self.aperture,self.lmin,self.lmax))
    def toString(self):
        return "%i %.10e %.10e %.10e \n "%(self.npix,self.aperture,self.lmin,self.lmax)
        
class imParams(object):
    def __init__(self,npix,side):
        self.npix = npix
        self.side = side
    def write(self,f):
        f.write("%i %.10e \n "%(self.npix,self.side))
    def toString(self):
        return "%i %.10e \n "%(self.npix,self.side)

class cubeParams(object):
    def __init__(self,cube_lbda_npix,cube_image_npix,cube_lmin,cube_lmax,cube_side):
        self.cube_lbda_npix  = cube_lbda_npix  
        self.cube_image_npix = cube_image_npix 
        self.cube_lmin       = cube_lmin       
        self.cube_lmax       = cube_lmax       
        self.cube_side       = cube_side       
    def write(self,f):
        f.write("%i %i %.10e %.10e %.10e \n "%(self.cube_lbda_npix,self.cube_image_npix,self.cube_lmin,self.cube_lmax,self.cube_side))
    def toString(self):
        return "%i %i %.10e %.10e %.10e \n "%(self.cube_lbda_npix,self.cube_image_npix,self.cube_lmin,self.cube_lmax,self.cube_side)


# simple functions to handle the 1000 different names of center and radius for different domains ...
# -> these assume a single spherical domain... 
def define_computational_domain(center,radius):
    from collections import OrderedDict
    return OrderedDict([('comput_dom_type','sphere'),('comput_dom_pos',center),('comput_dom_rsp',radius)])

def define_emission_domain(center,radius):
    from collections import OrderedDict
    return OrderedDict([('emission_dom_type','sphere'),('emission_dom_pos',center),('emission_dom_rsp',radius)])

def define_stellar_emission_domain(center,radius):
    from collections import OrderedDict
    return OrderedDict([('star_dom_type','sphere'),('star_dom_pos',center),('star_dom_rsp',radius)])
                            
def define_domain_decomposition(center,radius):
    from collections import OrderedDict
    return OrderedDict([('decomp_dom_type','sphere'),('decomp_dom_ndomain',1),\
                            ('decomp_dom_xc',center[0]),('decomp_dom_yc',center[1]),('decomp_dom_zc',center[2]),\
                            ('decomp_dom_rsp',radius)])

def define_photonsFromStarsParams(sedDir,spec_type='Mono',spec_mono_lambda0=1500.,spec_table_lmin_Ang=1170,spec_table_lmax_Ang=1260):
    from collections import OrderedDict
    PhotonsFromStarsParams = OrderedDict([ ('spec_type',spec_type),('spec_SSPdir',sedDir)])
    if spec_type == 'Mono':
        PhotonsFromStarsParams.update({'spec_mono_lambda0':spec_mono_lambda0})
    if spec_type == 'Table':
        PhotonsFromStarsParams.update({'spec_table_lmin_Ang':spec_table_lmin_Ang})
        PhotonsFromStarsParams.update({'spec_table_lmax_Ang':spec_table_lmax_Ang})
    return PhotonsFromStarsParams

def get_simulation_units(dir,ts):
    file  = '%s/output_%5.5i/info_%5.5i.txt'%(dir,ts,ts)
    f     = open(file,'r')
    lines = f.readlines()
    f.close()
    for l in lines:
        if 'unit_l' in l:
            ul = float(l.split()[2])
        if 'unit_d' in l:
            ud = float(l.split()[2])
        if 'unit_t' in l:
            ut = float(l.split()[2])
        if 'boxlen' in l:
            bl = float(l.split()[2])
    return ul,ud,ut,bl


def define_mock_params(center,mockSpec=False,mockImage=False,mockCube=False,\
                           flux_aperture_radius=0.,
                           spec_lbda_npix=0.,spec_aperture_radius=0.,spec_lbdmin=0,spec_lbdmax=0,\
                           image_npix=0,image_side=0,\
                           cube_image_npix=0,cube_lbda_npix=0,cube_lbdmin=0,cube_lbdmax=0,cubeSide_cu=0):
                           
    from collections import OrderedDict

    # define mock directions from theta = 0 to pi/2 (phi = 0) 
    nDirections = 6
    thetas = np.linspace(0,np.pi/2,nDirections)
    kxList,kyList,kzList = np.sin(thetas),np.zeros_like(thetas),np.cos(thetas)

    if mockSpec: 
        specpar = specParams(spec_lbda_npix,spec_aperture_radius,spec_lbdmin,spec_lbdmax)
    else:
        specpar = specParams(0,0,0,0)

    if mockImage:
        impar   = imParams(image_npix,image_side)
    else:
        impar   = imParams(0,0)

    if mockCube:
        cubepar = cubeParams(cube_lbda_npix,cube_image_npix,cube_lbdmin,cube_lbdmax,cubeSide_cu)
    else:
        cubepar = cubeParams(0,0,0,0,0)
        
    MockParameterFileContents = ""
    for kx,ky,kz in zip(kxList,kyList,kzList):
        MockParameterFileContents = MockParameterFileContents+"%.10e %.10e %.10e \n "%(kx,ky,kz)
        MockParameterFileContents = MockParameterFileContents+"%.10e %.10e %.10e \n "%(center[0],center[1],center[2])
        MockParameterFileContents = MockParameterFileContents+"%.10e \n"%(flux_aperture_radius)
        MockParameterFileContents = MockParameterFileContents+specpar.toString()
        MockParameterFileContents = MockParameterFileContents+impar.toString()
        MockParameterFileContents = MockParameterFileContents+cubepar.toString()
    mock_params = OrderedDict([ ('nDirections',nDirections),
                                ('MockParameterFileContents',MockParameterFileContents)])
    return mock_params

