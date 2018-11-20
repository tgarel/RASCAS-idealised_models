class mockobs(object):
    def __init__(self,Dir,FileName,ncpu,ICFile=None):
        self.Dir      = Dir
        self.FileName = FileName 
        self.ncpu     = ncpu

        import jphot as jp
        if ICFile is None:
            ICFile = FileName.replace('result.','')
            p = jp.photonlist("%s/%s"%(self.Dir,ICFile),'',load=False)
        else:
            p = jp.photonlist(ICFile,'',load=False)
        p.load_ic()
        self.nPhotPerPacket = p.nRealPhotons / p.nphoton
        import lya_utils as ly
        self.LumPerPacket = self.nPhotPerPacket * ly.h_cgs * ly.nu0 

        self.imsize_cu, self.image_npix, self.imtot = self.__read_image()
        self.imtot = self.imtot * self.LumPerPacket # [erg/s]
        self.lmin,self.lmax,self.spec_npix,self.spectot = self.__read_spectrum()

        
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
            lmin,lmax = f.read_reals('d')
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
        l = linspace(self.lmin,self.lmax,num=self.spec_npix)  # [A]
        plt.plot(l,self.spectot)
        plt.savefig(plotFile)



