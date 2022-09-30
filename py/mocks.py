import lya_utils as ly
import jphot as jp
import numpy as np

class mockobs(object):

    #
    # This class handles mock observations output by RASCAS
    #

    def __init__(self,Dir,FileName,ICFile,load_flux=False,load_spectrum=False,load_image=False,load_cube=False,unit_l_arcsec=1.0,lumdist_cm=1.0,redshift=0.0,idirection=1):
        self.Dir      = Dir
        self.FileName = FileName 
            
        p = jp.photonlist("%s/%s"%(self.Dir,ICFile),'',load=False)
        nRealPhotons = p.get_nRealPhotons() 
        nPhotons     = p.get_nphoton()
        self.nPhotPerPacket = nRealPhotons / nPhotons
        
        if load_flux:
            self.aperture,self.flux,self.hnu = self.__read_flux(index=idirection)
            self.LumPerPacket = self.nPhotPerPacket * self.hnu
            self.flux = self.flux * self.LumPerPacket # [erg/s]

        if load_spectrum:  # restframe spectrum
            self.spec_lmin,self.spec_lmax,self.spec_npix,self.spec = self.__read_spectrum(index=idirection)
            dl = (self.spec_lmax - self.spec_lmin) / self.spec_npix  # [A] 
            l = np.arange(self.spec_lmin,self.spec_lmax,dl) + 0.5*dl # [A]
            self.spec_lbda_Angstrom  = l
            self.spec_dlbda_Angstrom = dl
            l = l * 1e-8  # [cm]
            energy = ly.h_cgs * ly.clight / l      # [erg]
            energy = energy * self.nPhotPerPacket  # [erg / s / phot packet]
            energy = energy / dl                   # [erg / s / A / phot packet]
            self.spec = energy * self.spec         # [erg / s/ A]

        if load_image:
            self.imsize_cu, self.image_npix, self.imtot, self.imhnu = self.__read_image(index=idirection)
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
            energy = energy * (1. + redshift) # because we computed obs-frame energy, so correct for a factor (1+z) taken into account by lum_dist. 
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
        f = ff('%s/%s.flux'%(self.Dir,self.FileName),'r')
        for j in range(index):
            aperture,flux,hnu = f.read_reals('d')
        f.close()
        return aperture,flux,hnu

    def __read_spectrum(self,index=1):
        from scipy.io import FortranFile as ff
        f = ff('%s/%s.spectrum'%(self.Dir,self.FileName))
        for j in range(index):
            n = f.read_ints()[0]
            aperture,lmin,lmax = f.read_reals('d')
            spec = f.read_reals('d')
        f.close()
        return lmin,lmax,n,spec

    def __read_image(self,index=1):
        from scipy.io import FortranFile as ff
        f = ff('%s/%s.image'%(self.Dir,self.FileName))
        for j in range(index):
            n        = f.read_ints()[0]
            imsize   = f.read_reals('d')[0]
            imcenter = f.read_reals('d')
            im       = f.read_reals('d')
            imhnu   = f.read_reals('d')[0]
        im = im.reshape((n,n))
        f.close()
        return imsize,n,im,imhnu

    def __read_cube(self,index=1):
        from scipy.io import FortranFile as ff
        f = ff('%s/%s.cube'%(self.Dir,self.FileName))
        for j in range(index):
            nlbda,nxy = f.read_ints()
            lmin,lmax,imsize = f.read_reals('d')
            center = f.read_reals('d')
            cube = f.read_reals('d')
        cube = cube.reshape((nxy,nxy,nlbda))
        f.close()
        return nlbda,nxy,lmin,lmax,imsize,cube


    def show_spectrum(self,plotFile=None):
        import lya_utils as ly
        import matplotlib
        if plotFile is not None:
            matplotlib.use('Agg')
        from matplotlib import pyplot as plt
        from numpy import linspace
        plt.figure()
        l = linspace(self.spec_lmin,self.spec_lmax,num=self.spec_npix)  # [A]
        plt.plot(l,self.spec)
        if plotFile is not None:
            plt.savefig(plotFile)
        else:
            plt.show()
    
    
    def show_image(self,plotFile=None,vmin=None,vmax=None,showCB=True,smooth=False,smooth_scale_pix=2,noiseLevel=0):
        from matplotlib import pyplot as plt
        from matplotlib.colors import LogNorm
        image = self.imtot
        if smooth:
            import scipy.ndimage as ndimage
            image = ndimage.gaussian_filter(image,sigma=(smooth_scale_pix,smooth_scale_pix))
        if noiseLevel > 0: 
            image = image + np.random.normal(0.0,noiseLevel,image.shape)
        plt.figure(figsize=(6,5))
        plt.imshow(image,norm=LogNorm(vmin=vmin,vmax=vmax))
        if showCB :
            plt.colorbar()
        plt.tight_layout()
        if plotFile is not None:
            plt.savefig(plotFile)
        else:
            plt.show()


    def show_cube_image(self,ax=None,vmin=None,vmax=None,smooth=False,smooth_sigma_arcsec=0.285):
        # todo: use previous show_image 
        from numpy import sum
        import matplotlib.pyplot as plt
        from matplotlib.colors import LogNorm
        from mpl_toolkits.axes_grid1.inset_locator import inset_axes
        image = sum(self.cube,axis=2) * self.cube_dlbda_Angstrom
        if smooth:
            import scipy.ndimage as ndimage
            smooth_scale_pix = smooth_sigma_arcsec / self.cube_dx_arcsec
            image = ndimage.gaussian_filter(image,sigma=(smooth_scale_pix,smooth_scale_pix))
            cs = plt.contour(image,levels=[1e-12,1.5e-12,2.e-12],origin='lower',
                                 extent=(self.cube_xmin_arcsec,self.cube_xmax_arcsec,self.cube_xmin_arcsec,self.cube_xmax_arcsec),
                                 colors=['limegreen','deepskyblue','red'],linestyles=['--','-','-'],
                                 linewidths=1)
        if ax is None:
            ax=plt.subplot(1,1,1)
        im= plt.imshow(image,norm=LogNorm(vmin=vmin,vmax=vmax),origin='lower', interpolation='nearest', \
                    extent=[self.cube_xmin_arcsec,self.cube_xmax_arcsec,self.cube_xmin_arcsec,self.cube_xmax_arcsec],cmap='Greys')
        plt.xlabel('arcsec')
        plt.ylabel('arcsec')
        axins = inset_axes(ax,width="3%",height="100%",loc=4,borderpad=0)
        cbar = plt.colorbar(im,cax=axins)
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

    # TODO: save_cube_fits


