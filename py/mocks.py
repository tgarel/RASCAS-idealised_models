import lya_utils as ly
import jphot as jp
import numpy as np

class mockobs(object):
    def __init__(self,Dir,FileName,load_flux=False,load_spectrum=False,load_image=False,load_cube=False,unit_l_arcsec=1.0,lumdist_cm=1.0,redshift=0.0,cube_h5_file=None,idirection=1):
        self.Dir      = Dir
        self.FileName = FileName 
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
            self.imtot = self.imtot * self.LumPerPacket # [erg/s]
        if load_cube:  # obs-frame cube
            self.cube_nlbda,self.cube_nxy,self.cube_lmin,self.cube_lmax,self.cube_imsize,self.cube = self.__read_cube(index=idirection)
            self.cube_lmin = self.cube_lmin * (1+redshift)
            self.cube_lmax = self.cube_lmax * (1+redshift)
            dl = (self.cube_lmax - self.cube_lmin) / self.cube_nlbda  # [A] 
            l = np.arange(self.cube_lmin,self.cube_lmax,dl) + 0.5*dl  # [A]
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
    
    def show_image(self,plotFile,vmin=1e-15,vmax=10,showCB=True):
        import matplotlib
        matplotlib.use('Agg')
        from matplotlib import pyplot as plt
        from matplotlib.colors import LogNorm
        plt.figure(figsize=(10,10))
        plt.imshow(self.imtot,norm=LogNorm(),vmin=vmin,vmax=vmax)
        if showCB :
            plt.colorbar()    
        plt.savefig(plotFile)

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

    def show_cube_spec(self,smooth=False,smooth_sigma_angstrom=0.5):
        from numpy import sum
        import matplotlib.pyplot as plt
        #plt.figure()
        spectrum = sum(self.cube,axis=(0,1)) * self.cube_dx_arcsec * self.cube_dx_arcsec
        if smooth:
            import scipy.ndimage as ndimage
            smooth_scale_pix = smooth_sigma_angstrom / self.cube_dlbda_Angstrom
            spectrum = ndimage.gaussian_filter(spectrum,sigma=(smooth_scale_pix))

        plt.plot(self.cube_lbda_Angstrom,spectrum)
        plt.axvline(1215.67*4.20039,color='red',alpha=0.5)
        plt.xlabel(r'$\lambda_{\rm obs} \ [\AA]$',fontsize=15)
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
