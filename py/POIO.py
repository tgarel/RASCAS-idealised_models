import lya_utils as ly
import jphot as jp
import numpy as np

class mockobs(object):
    def __init__(self,Dir,FileName,ncpu,load_spectrum=False,load_image=False,load_cube=False,unit_l_arcsec=1.0,lumdist_cm=1.0,redshift=0.0,cube_h5_file=None,idirection=1):
        self.Dir      = Dir
        self.FileName = FileName 
        self.ncpu     = ncpu
        ICFile = FileName.replace('result.','')
        p = jp.photonlist("%s/%s"%(self.Dir,ICFile),'',load=False)
        nRealPhotons = p.get_nRealPhotons()
        nPhotons     = p.get_nphoton()
        self.nPhotPerPacket = nRealPhotons / nPhotons
        self.LumPerPacket = self.nPhotPerPacket * ly.h_cgs * ly.nu0 

        if load_image:
            self.imsize_cu, self.image_npix, self.imtot = self.__read_image()
            self.imtot = self.imtot * self.LumPerPacket # [erg/s]

        if load_spectrum:  # restframe spectrum
            self.spec_lmin,self.spec_lmax,self.spec_npix,self.spec = self.__read_spectrum()
            dl = (self.spec_lmax - self.spec_lmin) / self.spec_npix  # [A] 
            l = np.arange(self.spec_lmin,self.spec_lmax,dl) + 0.5*dl # [A]
            self.spec_lbda_Angstrom  = l
            self.spec_dlbda_Angstrom = dl
            l = l * 1e-8  # [cm]
            energy = ly.h_cgs * ly.clight / l      # [erg]
            energy = energy * self.nPhotPerPacket  # [erg / s / phot packet]
            energy = energy / dl                   # [erg / s / A / phot packet]
            self.spec = energy * self.spec         # [erg / s/ A]

        if load_cube:  # obs-frame cube
            if cube_h5_file is not None:
                self.load_cube_h5(cube_h5_file)
            else:
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
            
    def __read_cube(self,index=1):
        from scipy.io import FortranFile as ff
        first = True
        for i in range(1,self.ncpu):
            f = ff('%s/%s_cube.%5.5i'%(self.Dir,self.FileName,i))
            for j in range(1,index+1): 
                nlbda,nxy = f.read_ints()
                lmin,lmax,imsize = f.read_reals('d')
                center = f.read_reals('d')
                cube = f.read_reals('d')
            cube = cube.reshape((nxy,nxy,nlbda))
            f.close()
            if first:
                cubetot = cube
                first=False
            else:
                cubetot = cubetot + cube
        return nlbda,nxy,lmin,lmax,imsize,cubetot

    
    def __read_image(self):
        from scipy.io import FortranFile as ff
        first = True
        for i in range(1,self.ncpu):
            f = ff('%s/%s_image.%5.5i'%(self.Dir,self.FileName,i))
            n = f.read_ints()[0]
            imsize = f.read_reals('d')[0]
            imcenter = f.read_reals('d')
            im = f.read_reals('d')
            im=im.reshape((n,n))
            f.close()
            if first:
                imtot = im
                first=False
            else:
                imtot = imtot + im
        return imsize,n,imtot 


    def __read_spectrum(self):
        from scipy.io import FortranFile as ff
        first = True
        for i in range(1,self.ncpu):
            f = ff('%s/%s_spectrum.%5.5i'%(self.Dir,self.FileName,i))
            n = f.read_ints()[0]
            aperture,lmin,lmax = f.read_reals('d')
            spec = f.read_reals('d')
            f.close()
            if first:
                spectot = spec
                first=False
            else:
                spectot = spectot + spec
        return lmin,lmax,n,spectot

    
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

        
    def show_spectrum(self,plotFile):
        import lya_utils as ly
        import matplotlib
        matplotlib.use('Agg')
        from matplotlib import pyplot as plt
        from numpy import linspace
        plt.figure()
        l = linspace(self.spec_lmin,self.spec_lmax,num=self.spec_npix)  # [A]
        plt.plot(l,self.spec)
        plt.savefig(plotFile)


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

    def load_cube_h5(self,filename):
        import h5py 
        f = h5py.File(filename,'r')
        self.cube=f['cube'][:,:,:]
        self.cube_nlbda = f['nlbda'].value
        self.cube_nxy = f['nxy'].value
        self.cube_lbda_Angstrom = f['lbda_Angstrom'][:]
        self.cube_dlbda_Angstrom = f['dlbda_Angstrom'].value
        self.cube_x_arcsec = f['x_arcsec'][:]
        self.cube_y_arcsec = ['y_arcsec'][:]
        self.cube_xmin_arcsec = f['xmin_arcsec'].value
        self.cube_xmax_arcsec = f['xmax_arcsec'].value
        self.cube_ymin_arcsec = f['ymin_arcsec'].value
        self.cube_ymax_arcsec = f['ymax_arcsec'].value
        self.cube_dx_arcsec = f['dx_arcsec'].value
        f.close()

    def write_fits_cube(self,filename):
        ## DOES NOT WORK ... 
        import astropy.io.fits as fits
        import astropy.wcs as pywcs
        from datetime import datetime
        # write a datacube into a FITS file readable by HSIM
        # primary header + 1 DATA extension
        primary_header=fits.Header()
        data_header=fits.Header()
        wcs = pywcs.WCS(naxis=3)
        wcs.wcs.crpix = np.array([(self.cube_nxy+1)/2.,(self.cube_nxy+1)/2.,1.0])
        wcs.wcs.crval = np.array([0.,0.,self.cube_lmin])
        wcs.wcs.ctype = ['x', 'y','wavelength']
        wcs.wcs.cunit = ['arcsec', 'arcsec','Angstrom']
        wcs.wcs.crota=[0.0,0.0,0.0]
        wcs.wcs.cdelt=[self.cube_dx_arcsec,self.cube_dx_arcsec,self.cube_dlbda_Angstrom]
        hdrwcs = wcs.to_header()
        primary_header['date'] = (str(datetime.now()), 'creation date')
        primary_header['author'] = ('HSIMREADY', 'origin of the file')
        hdulist = fits.HDUList([fits.PrimaryHDU(header=primary_header)])
        keys = set(data_header.keys()) - set(hdrwcs.keys())
        for card in data_header.cards:
           if card.keyword not in keys:
               continue
           hdrwcs[card.keyword] = (card.value, card.comment)

        hdrwcs['FUNITS']='erg/s/cm2/A/arcsec2'
        hdulist.append(fits.ImageHDU(name='DATA',data=self.cube,header=hdrwcs))
        hdulist.writeto(filename, clobber=True, output_verify='silentfix')


