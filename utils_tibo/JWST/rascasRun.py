 
#---------------------------------------------------------------------------------------------------
# RASCAS generic run class 
#---------------------------------------------------------------------------------------------------

# to do's :
# - add power law options in __generate_PhotTables and __setup_PhotonsFromStars
# - __setup_PhotonsFromStars assumes a spherical domain ... fix that ? (at least an error message ?)
#

# wish list :
# - add a feature to re-start a run ... (and check whether a run is completed or not). 
# 

class RascasSurvey(object):

    def __init__(self,surveyName,rascasDir,DomDumpDir,ramsesDir,ramsesTimestep,rascas_f90):
        """
        - Generate directory structure for rascas run (if not already there)
        - whatelse...
        """
        self.rascasDir      = rascasDir      # where rascas runs
        self.rascas_f90     = rascas_f90     # where rascas code/exec are
        self.surveyName     = surveyName     # name of survey
        self.surveyDir      = "%s/%s"%(self.rascasDir,self.surveyName) # outputs of RASCAS are in rascasDir/surveyName/
        self.DomDumpDir     = "%s/%s"%(self.rascasDir,DomDumpDir)      # outputs of CreateDomDump are in rascasDir/DomDumpDir/
        self.ramsesDir      = ramsesDir      # where ramses outputs are 
        self.ramsesTimestep = ramsesTimestep # number of simulation snapshot to process

        

    def setup_Lya_survey(self,ComputationalDomain, DomainDecomposition, EmissionDomain,\
                             nphotons = 10000, seedStep=100, doRecombs=True, doColls=False, ramses_params=None, \
                             gas_composition_params=None, HI_model_params=None, dust_model_params=None):
        
        """
        Setup a Lya survey ... 
        NB: do either recombinations or collisions. The surveyName should be different for each. 

        """

        # create directory structure for rascas run if it does not exist yet
        self.__make_directories()

        # create a CreateDomDump parameter file if it does not exist yet
        self.__setup_CreateDomDump(ComputationalDomain,DomainDecomposition,\
                                       ramses_params,gas_composition_params,\
                                       reading_method='hilbert')

        # Create Photon ICs with LyaPhotonsFromGas 
        runNum, PhotonICFile = self.__setup_LyaPhotonsFromGas(seedStep,EmissionDomain,nphotons,ramses_params,\
                                      doRecombs,doColls)

        # Create Rascas parameter file 
        self.__setup_rascas(runNum, PhotonICFile,\
                           gas_composition_params,HI_model_params,dust_model_params)






        
        
    def setup_broad_band_survey(self, PhotometricTableParams, ComputationalDomain, DomainDecomposition, StellarEmissionDomain,
                                    nphotons = 10000, seedStep=100,
                                    ramses_params=None, gas_composition_params=None,
                                    HI_model_params=None, dust_model_params=None):

        """

        Set up a rascas survey. Used to create a new survey or to complement an existing one. 

        Input parameters:
        -----------------

        ComputationalDomain : a dictionary which defines the computational domain
        DomainDecomposition : a dictionary which defines the sub-domains
        StellarEmissionDomain 

        Relative to creation of photometric tables : PhotometricTableParams is a dictionary containing:
        - sedDir : directory where SED models are.
        - sedModel : which SED library to use 
        - lbdamin_Angstrom : the min wavelength [A]
        - lbdamax_Angstrom : the max wavelength [A]
        - photTableDir : absolute path to directory where the photometric tables will be written.
        - method : either 'TabulatedSpectra' or 'PowerLaw'
        - dlambda : used for method=='TabulatedSpectra'

        Relative to PhotonsFromStars:
        - star_dom_type ... has to be a sphere for now. 
        - star_dom_pos 
        - star_dom_rsp
        - seedStep : the random seed will be defined as the minimum seed found in existing 
                     parameter files (0 if no file is found) minus seedStep. 
        - nphotons : the number of MC photons to define

        

        """

        # ----------------------------------------------------------------------
        # create directory structure for rascas run if it does not exist yet
        # ----------------------------------------------------------------------
        self.__make_directories()

        # ----------------------------------------------------------------------            
        # create a CreateDomDump parameter file if it does not exist yet
        # ----------------------------------------------------------------------
        self.__setup_CreateDomDump(ComputationalDomain,DomainDecomposition,\
                                       ramses_params,gas_composition_params,\
                                       reading_method='hilbert')

        # ----------------------------------------------------------------------
        # Generate photometric tables if they do not exist
        # ----------------------------------------------------------------------
        PhotTableFile = self.__generate_PhotTables(PhotometricTableParams)
        
        # ----------------------------------------------------------------------
        # Create PhotonsFromStars parameter file 
        # ----------------------------------------------------------------------
        runNum, PhotonICFile = self.__setup_PhotonsFromStars(seedStep,StellarEmissionDomain,\
                                                                 PhotTableFile,PhotometricTableParams,\
                                                                 nphotons,ramses_params)

        # ----------------------------------------------------------------------
        # Create Rascas parameter file 
        # ----------------------------------------------------------------------
        self.__setup_rascas(runNum, PhotonICFile,\
                           gas_composition_params,HI_model_params,dust_model_params)


        
    def load(self,savePhots=False,saveOnlyEscaped=True,maxRuns=100,startRun=1):

        import os
        import numpy as np
        # rascas import :
        import jphot as jp

        self.noResults = False
        self.noICs     = False
        # check which result files are available
        resNums,icNums = [],[]
        for irun in range(startRun,maxRuns+1):
            RascasOutputFile = '%s/%3.3i.res'%(self.surveyDir,irun)
            print(RascasOutputFile) ###
            if os.path.exists(RascasOutputFile):
                resNums.append(irun)
            PhotonICFile = '%s/%3.3i.IC'%(self.surveyDir,irun)
            if os.path.exists(PhotonICFile):
                icNums.append(irun)
        if len(resNums) == 0: # no result file
            self.noResults = True

        if len(icNums) == 0:
            self.noICs = True

        if self.noResults and self.noICs :
            return

        # read total luminosity (i.e. nb of real photons sent).
        irun = icNums[0]
        PhotonICFile = '%s/%3.3i.IC'%(self.surveyDir,irun)
        p = jp.photonlist(PhotonICFile,'',load=False)
        self.nRealPhotons = p.get_nRealPhotons()

        # read result files 
        if not self.noResults:
            if savePhots:
                first = True
                for irun in resNums:
                    PhotonICFile = '%s/%3.3i.IC'%(self.surveyDir,irun)
                    RascasOutputFile = '%s/%3.3i.res'%(self.surveyDir,irun)
                    if first:
                        p = jp.photonlist(PhotonICFile,RascasOutputFile)
                        first = False
                    else:
                        p.append(jp.photonlist(PhotonICFile,RascasOutputFile))          
                escaped   = np.where(p.status == 1)[0]
                self.p    = p
                self.pesc = p.extract_sample(escaped)
                nphotTot  = self.p.nphoton

            elif saveOnlyEscaped:
                first = True
                nphotTot = 0
                for irun in resNums:
                    PhotonICFile = '%s/%3.3i.IC'%(self.surveyDir,irun)
                    RascasOutputFile = '%s/%3.3i.res'%(self.surveyDir,irun)
                    p = jp.photonlist(PhotonICFile,RascasOutputFile)
                    nphotTot = nphotTot + p.nphoton
                    escaped = np.where(p.status == 1)[0]
                    if first:
                        pesc = p.extract_sample(escaped)                    
                        first = False
                    else:
                        pesc.append(p.extract_sample(escaped))            
                self.p    = None
                self.pesc = pesc
            self.nphotTot = nphotTot


    def make_mock_cube(self,k,snap,thetamax=20.,nxybins=100,xc=None,yc=None,zc=None,pixsize=0.2,lbdapix=1.25,lbdamin_restframe=1200,lbdamax_restframe=1230):
        # Inputs : 
        # - k : direction of observation (3d vector)
        # - snap : snapshot object (contains redshift, unit_l_arcsec, ...)
        # - thetamax : selection angle around k [deg]
        # - photFlux : flux carried by each photon packet [erg/s/cm2]
        # - pixsize (optional) : size of spatial pixels [arcsec]
        # - lbdapix (optional) : size of spectral pixels [A]

        import numpy as np
        import lya_utils as lya
        import jphot as jp

        # select photons going along k, +/- thetamax
        costheta = self.pesc.kx * k[0] + self.pesc.ky * k[1] + self.pesc.kz * k[2]
        ii       = np.where(costheta > np.cos(thetamax*np.pi/180.))[0]
        p        = self.pesc.extract_sample(ii)
        domega   = 2.*np.pi*(1-np.cos(thetamax*np.pi/180.))
        percm2   = 1. / (domega*snap.lumDist_cm**2)  # /cm2
        
        # define observer-frame wavelengths
        lbda_rest = lya.clight/p.nu*1e8  # [A]
        lbda      = lbda_rest * (1.+snap.redshift)
        lbdamin   = lbdamin_restframe * (1.+snap.redshift)
        lbdamax   = lbdamax_restframe * (1.+snap.redshift)
        
        # define projected coordinates in arcsec
        if xc is None: xc = np.median(p.x)
        if yc is None: yc = np.median(p.y)
        if zc is None: zc = np.median(p.z)
        print("xc,yc,zc = %f, %f, %f"%(xc,yc,zc))
        x,y = jp.projected_positions(k,p.x,p.y,p.z,xc,yc,zc)   # now centered around 0
        x   = x * snap.unit_l_arcsec
        y   = y * snap.unit_l_arcsec
        
        # define sampling  
        # define effective extent of the map in arcsec.
        xmin = - 0.5 *nxybins * pixsize
        xmax = xmin + nxybins*pixsize
        ymin = - 0.5 *nxybins * pixsize
        ymax = ymin + nxybins*pixsize
        # define effective extent of spectral direction
        dl = lbdamax - lbdamin
        nlbins = int(dl / lbdapix)
        lmin = lbdamin
        lmax = lbdamin + nlbins*lbdapix
        
        # construct cube
        ener = (p.nRealPhotons/float(p.nphoton) * lya.h_cgs * lya.clight*1e8)/lbda # energy of each photon packet, observer-frame [erg/s]
        cube,xx = np.histogramdd(np.array([lbda,x,y]).T, normed=False,
                                     weights=ener,
                                     bins=(nlbins,nxybins,nxybins),
                                     range=((lmin,lmax),(xmin,xmax),(ymin,ymax)))
        # normalize
        pixarea_arcsec2 = pixsize*pixsize  # [arcsec2]
        cube = cube / pixarea_arcsec2 / lbdapix * percm2 # erg/s/cm2 / arcsec2 / A
        
        return cube,lmin,nlbins,xmin,ymin,nxybins
        
        
    def mock_image(self,k=[0,0,1],thetamax=20,nbins=100,xc=None,yc=None,zc=None,rmax=None,lumDist_cm=1.,escaped=True):

        import numpy as np
        import lya_utils as lya
        import jphot as jp

        # select photons going along k within +/- thetamax
        if escaped:
            costheta = self.pesc.kx * k[0] + self.pesc.ky * k[1] + self.pesc.kz * k[2]
            ii       = np.where(costheta > np.cos(thetamax*np.pi/180.))[0]
            p        = self.pesc.extract_sample(ii)
        else:
            costheta = self.p.kx * k[0] + self.p.ky * k[1] + self.p.kz * k[2]
            ii       = np.where(costheta > np.cos(thetamax*np.pi/180.))[0]
            p        = self.p.extract_sample(ii)

        # define image center (in 3D)
        if xc is None: xc = np.median(p.x)
        if yc is None: yc = np.median(p.y)
        if zc is None: zc = np.median(p.z)
        if rmax is None: rmax = max([xc - p.x.min(),p.x.max()-xc,yc - p.y.min(), p.y.max() - yc,zc - p.z.min(), p.z.max() - zc])

        print(xc,yc,zc,rmax)
            
        # get coords in direction of observation
        x,y = jp.projected_positions(k,p.x,p.y,p.z,xc,yc,zc)   # now centered around 0

        # compute the map
        nphot_per_packet = self.pesc.nRealPhotons / self.pesc.nphoton
        ener    = nphot_per_packet * lya.h_cgs * p.nu  # [erg/s/MC phot]
        h,hx,hy = np.histogram2d(x,y,bins=(nbins,nbins),weights=ener,range=[(-rmax,rmax),(-rmax,rmax)])  # erg/s in each pixel

        # solid angle of selection : 
        domega   = 2.*np.pi*(1-np.cos(thetamax*np.pi/180.))
        h = h / (domega*lumDist_cm**2)  # erg/s/cm2
        
        # -> divide by pixel size to get a surface brightness. 

        return h,hx,hy 

    
    def imagexy(self,nbins=100,xc=None,yc=None,zc=None,rmax=None,escaped=True,initial_pos=False):

        import numpy as np
        import lya_utils as lya

        if escaped:
            p = self.pesc
        else:
            p = self.p

        if initial_pos:
            p.x = p.x_ic
            p.y = p.y_ic
            p.z = p.z_ic
        # define image center (in 3D)
        if xc is None: xc = np.median(p.x)
        if yc is None: yc = np.median(p.y)
        if zc is None: zc = np.median(p.z)
        if rmax is None: rmax = max([xc - p.x.min(),p.x.max()-xc,yc - p.y.min(), p.y.max() - yc,zc - p.z.min(), p.z.max() - zc])

        # get coords in direction of observation
        x,y = p.x,p.y

        # compute the map
        nphot_per_packet = p.nRealPhotons / p.nphoton
        ener    = nphot_per_packet * lya.h_cgs * p.nu  # [erg/s/MC phot]
        h,hx,hy = np.histogram2d(x,y,bins=(nbins,nbins),weights=ener,range=[(xc-rmax,xc+rmax),(yc-rmax,yc+rmax)])  # erg/s in each pixel
        # -> divide by pixel size to get a surface brightness. 

        return h,hx,hy 

    
    def abs_mag(self,no_dust=False):

        # compute absolute magnitude for photons in a monochromatic survey

        import numpy as np
        import lya_utils as lya

        Tenpc_cm = 3.086e19   # 10pc in cm... 
        #m = self.nRealPhotons  # [erg/s/A]  !!! This is due to a bug in CPT ... 
        m  = self.nRealPhotons * np.mean(self.pesc.nu) * lya.h_cgs # erg/s/A 
        m = m * 1e8  # [erg/s/cm]
        wavelength_Angstrom = lya.clight / np.mean(self.pesc.nu) * 1e8
        m = m * (wavelength_Angstrom*1e-8)**2 / lya.clight  # [erg/s/Hz]
        m = m / (4.*np.pi*Tenpc_cm**2) # [erg/s/Hz/cm2]
        if no_dust:
            ext_factor = 1.
        else:
            ext_factor = float(self.pesc.nphoton) / float(self.nphotTot)
        m = m * ext_factor # [erg/s/Hz/cm2]
        m = -2.5 * np.log10(m) - 48.6 # AB mag
        print("fesc = ")
        print(ext_factor)
        return m

    
    def abs_mag_angular_map(self,nside=24):

        # compute absolute magnitude for photons in a monochromatic survey

        import numpy as np
        import lya_utils as lya
        import healpy as hp

        Tenpc_cm = 3.086e19   # 10pc in cm... 
        m = self.nRealPhotons  # [erg/s/A]  !!! This is due to a bug in CPT ... 
        m = m * 1e8  # [erg/s/cm]
        wavelength_Angstrom = lya.clight / np.mean(self.pesc.nu) * 1e8
        m = m * (wavelength_Angstrom*1e-8)**2 / lya.clight  # [erg/s/Hz]
        m = m / (4.*np.pi*Tenpc_cm**2) # [erg/s/Hz/cm2]
        nphot_map = self.MCPhotCount_angular_map(nside)
        npix = float(hp.nside2npix(nside))
        ext_factor = nphot_map / float(self.nphotTot) * npix
        m = m * ext_factor # [erg/s/Hz/cm2]
        m = -2.5 * np.log10(m) - 48.6
        
        return m

    
    def Luminosity(self,no_dust=False,mono=False):
        # compute Luminosity in a Line survey

        import numpy as np
        import lya_utils as lya

        if mono: # (bugged thing ...)
            L = self.nRealPhotons * (lya.clight/np.mean(self.pesc.nu)*1e8) # erg/s/A -> erg/s
            if not no_dust:
                L = L * (float(self.pesc.nphoton) / self.nphotTot)

        else:  # general case ... 
            nphot_per_MC = self.nRealPhotons / self.nphotTot
            if no_dust:
                L = nphot_per_MC * lya.h_cgs * np.sum(self.p.nu)    # erg/s
            else:
                L = nphot_per_MC * lya.h_cgs * np.sum(self.pesc.nu)

        return L

    def Luminosity_within_shell(self,xc,rin,rout,intrinsic_emission=True,initial_positions=True,mono=False):

        import numpy as np
        import lya_utils as lya
        # count photons within sphere
        if intrinsic_emission:
            if initial_positions:
                d2 = (self.p.x_ic - xc[0])**2 + (self.p.y_ic - xc[1])**2 + (self.p.z_ic - xc[2])**2
            if not initial_positions:
                d2 = (self.p.x - xc[0])**2 + (self.p.y - xc[1])**2 + (self.p.z - xc[2])**2
        if not intrinsic_emission:
            if initial_positions:
                d2 = (self.pesc.x_ic - xc[0])**2 + (self.pesc.y_ic - xc[1])**2 + (self.pesc.z_ic - xc[2])**2
            if not initial_positions:
                d2 = (self.pesc.x - xc[0])**2 + (self.pesc.y - xc[1])**2 + (self.pesc.z - xc[2])**2
        ii = np.where((d2 < rout*rout)&(d2>rin*rin))[0]

        # compute Luminosity in a Line survey
        if mono: # (bugged thing ...)
            L = self.nRealPhotons * (lya.clight/np.mean(self.pesc.nu)*1e8) # erg/s/A -> erg/s
            if not no_dust:
                L = L * (float(self.pesc.nphoton) / self.nphotTot)
            L = L * float(len(ii))/len(d2)
        else:  # general case ... 
            nphot_per_MC = self.nRealPhotons / self.nphotTot
            if intrinsic_emission:
                L = nphot_per_MC * lya.h_cgs * np.sum(self.p.nu[ii])    # erg/s
            else:
                L = nphot_per_MC * lya.h_cgs * np.sum(self.pesc.nu[ii])

        return L

    def MCPhotCount_angular_map(self,nside):
        import numpy as np
        import lya_utils as lya
        import healpy as hp
        ppix_line = hp.pixelfunc.vec2pix(nside,self.pesc.kx,self.pesc.ky,self.pesc.kz,nest=False)
        npix  = hp.nside2npix(nside) 
        carte = np.zeros(npix)
        for j in ppix_line:
            carte[j] = carte[j]+1.
        return carte

    def fesc_angular_map(self,nside):
        import numpy as np
        import lya_utils as lya
        import healpy as hp
        ppix_line = hp.pixelfunc.vec2pix(nside,self.pesc.kx,self.pesc.ky,self.pesc.kz,nest=False)
        npix  = hp.nside2npix(nside) 
        carte = np.zeros(npix)
        for j in ppix_line:
            carte[j] = carte[j]+1. 
        carte = carte / self.nphotTot * float(npix) * 100.  # escape fraction in %
        return carte

    def fesc_percentiles(self,nside=4,percentiles=[25,50,75]):
        import numpy as np
        carte = self.fesc_angular_map(nside)
        return np.percentile(carte,percentiles)
    

    def __make_directories(self):

        """

        Make the rascas directories. The imposed structure is as follows: 
        --- /rascasDir/surveyName -> contains parameter files for a rascas run, and will contain photon ICs and results
        --- /rascasDir/DomDumpDir -> contains parameter files for a CreateDomDump run and the results of that. 

        """
        
        import os
        if not os.path.exists(self.rascasDir):
            os.makedirs(self.rascasDir)
        if not os.path.exists(self.DomDumpDir):
            os.makedirs(self.DomDumpDir)
        if not os.path.exists(self.surveyDir):
            os.makedirs(self.surveyDir)


        
        
    def __setup_CreateDomDump(self,ComputationalDomain,DomainDecomposition,ramses_params,
                                  gas_composition_params,reading_method='hilbert'):

        """

        Create a CreateDomDump parameter file if it does not exist yet

        """

        import os
        from collections import OrderedDict
        import setupUtils as sU  # read/write parameter files in RASCAS format

        cddParamFile = "%s/CDD_params.conf"%(self.DomDumpDir)
        if not os.path.exists(cddParamFile):
            # Define CreateDomDump parameters
            center = ComputationalDomain['comput_dom_pos']
            if DomainDecomposition['decomp_dom_type'] == 'shell':
                cdd_params = OrderedDict([('DomDumpDir',self.DomDumpDir),('repository',self.ramsesDir),('snapnum',self.ramsesTimestep),
                                            ('comput_dom_type',ComputationalDomain['comput_dom_type']),
                                            ('comput_dom_pos',"%.16e, %.16e, %.16e"%(center[0],center[1],center[2])),
                                            ('comput_dom_rsp',"%.16e"%(ComputationalDomain['comput_dom_rsp'])),
                                            ('decomp_dom_type',DomainDecomposition['decomp_dom_type']),
                                            ('decomp_dom_ndomain',DomainDecomposition['decomp_dom_ndomain']),
                                            ('decomp_dom_xc',DomainDecomposition['decomp_dom_xc']),
                                            ('decomp_dom_yc',DomainDecomposition['decomp_dom_yc']),
                                            ('decomp_dom_zc',DomainDecomposition['decomp_dom_zc']),
                                            ('decomp_dom_rin',DomainDecomposition['decomp_dom_rin']),
                                            ('decomp_dom_rout',DomainDecomposition['decomp_dom_rout']),
                                            ('verbose','T'),('reading_method',reading_method)])
            if DomainDecomposition['decomp_dom_type'] == 'sphere':
                cdd_params = OrderedDict([('DomDumpDir',self.DomDumpDir),('repository',self.ramsesDir),('snapnum',self.ramsesTimestep),
                                            ('comput_dom_type',ComputationalDomain['comput_dom_type']),
                                            ('comput_dom_pos',"%.16e, %.16e, %.16e"%(center[0],center[1],center[2])),
                                            ('comput_dom_rsp',"%.16e"%(ComputationalDomain['comput_dom_rsp'])),
                                            ('decomp_dom_type',DomainDecomposition['decomp_dom_type']),
                                            ('decomp_dom_ndomain',DomainDecomposition['decomp_dom_ndomain']),
                                            ('decomp_dom_xc',DomainDecomposition['decomp_dom_xc']),
                                            ('decomp_dom_yc',DomainDecomposition['decomp_dom_yc']),
                                            ('decomp_dom_zc',DomainDecomposition['decomp_dom_zc']),
                                            ('decomp_dom_rsp',DomainDecomposition['decomp_dom_rsp']),
                                            ('verbose','T'),('reading_method',reading_method)])
            parameters = OrderedDict([('CreateDomDump',cdd_params),
                        ('gas_composition',gas_composition_params),
                        ('ramses',ramses_params)])
            sU.write_parameter_file(cddParamFile,parameters)

            #tibo
            fff = "%s.pbscript"%(cddParamFile)
            f = open(fff,'w')
            f.write("#PBS -S /bin/bash \n")
            f.write("#PBS -j oe   \n")
            f.write("#PBS -l select=1:ncpus=16 \n")
            f.write("#PBS -l place=scatter:excl   \n")
            f.write("#PBS -l walltime=24:00:00   \n")
            f.write("#PBS -q workq \n")
            f.write("cd $PBS_O_WORKDIR   \n")
            f.write("source /usr/share/modules/init/bash   \n")
            f.write("module load intel-icc-14/14.0.3.174   \n")
            f.write("module load intel-ifc-14/14.0.3.174   \n")
            f.write("module load intel-mkl-11/11.1.3   \n")
            f.write("module load mpt/2.10   \n")
            f.write("export OMP_NUM_THREADS=16 \n")
            f.write("set -x   \n")
            #f.write("/scratch/garel/rascas_sphinx/f90/CreateDomDump %s >& %s.log \n"%(cddParamFile,cddParamFile))
            f.write("%sCreateDomDump %s >& %s.log \n"%(self.rascas_f90,cddParamFile,cddParamFile))

            f.close()
            

            
    def __generate_PhotTables(self,PhotometricTableParams):
        """
        
        Generate Photometric tables if they are not found.

        """

        import os
        from collections import OrderedDict
        import CreatePhotTables as CPT  # Utils to generate photometric tables

        if not os.path.exists(PhotometricTableParams['photTableDir']):
            os.makedirs(PhotometricTableParams['photTableDir'])
        PhotTableFile = '%s/%s_PhotTable_%s.dat'%(PhotometricTableParams['photTableDir'],
                                                      self.surveyName,PhotometricTableParams['sedModel'])
        if not os.path.exists(PhotTableFile):
            ssp = CPT.readRamsesSEDs(sedDir="%s/%s"%(PhotometricTableParams['sedDir'],PhotometricTableParams['sedModel']))  # read SED models
            if PhotometricTableParams['method'] == 'TabulatedSpectra':
                CPT.gen_tabulated_spectra(PhotometricTableParams['lbdamin_Angstrom'],
                                              PhotometricTableParams['lbdamax_Angstrom'],
                                              ssp,PhotTableFile,
                                              dlambda=PhotometricTableParams['dlambda'],
                                              comment='Generated with %s/%s'%(PhotometricTableParams['sedDir'],\
                                              PhotometricTableParams['sedModel']))
            elif PhotometricTableParams['method'] == 'Monochromatic':
                CPT.gen_monochromatic_tables(PhotometricTableParams['lbda0_Angstrom'],ssp,PhotTableFile,\
                                             comment='Generated with %s/%s'%(PhotometricTableParams['sedDir'],\
                                             PhotometricTableParams['sedModel']))
            #elif PhotometricTableParams['method'] == 'PowerLaw':
            # --- TO BE DONE
                
        return PhotTableFile

    

    def __setup_PhotonsFromStars(self,seedStep,StellarEmissionDomain,PhotTableFile,PhotometricTableParams,nphotons,ramses_params):
        
        """
        Create PhotonsFromStars parameter file 
        
        """
        import os
        from collections import OrderedDict
        import setupUtils as sU  # read/write parameter files in RASCAS format

        # determine new file number and seed for random number generator
        OK = False
        runNum  = 1
        minseed = 0
        while not OK:
            pfsParamFile = "%s/%3.3i.PFS.conf"%(self.surveyDir,runNum)
            if os.path.exists(pfsParamFile):
                f = open(pfsParamFile,'r')
                lines = f.readlines()
                for l in lines:
                    if 'ranseed' in l:
                        seed = int(l.split()[2])
                        if seed < minseed:
                            minseed = seed
                        break
                runNum = runNum + 1
            else:
                OK = True
        runSeed = minseed - seedStep
        
        # define PhotonsFromStars parameters
        PhotonICFile = '%s/%3.3i.IC'%(self.surveyDir,runNum)
        sdp=StellarEmissionDomain['star_dom_pos']
        if PhotometricTableParams['method'] == 'TabulatedSpectra':
            pfs_params = OrderedDict([('outputfile',PhotonICFile),
                                        ('repository',self.ramsesDir),
                                        ('snapnum',self.ramsesTimestep),
                                        ('star_dom_type',StellarEmissionDomain['star_dom_type']),
                                        ('star_dom_pos',"%.16e %.16e %.16e"%(sdp[0],sdp[1],sdp[2])),
                                        ('star_dom_rsp','%.16e'%(StellarEmissionDomain['star_dom_rsp'])),
                                        ('weight_type','Table'),
                                        ('weight_input_file',PhotTableFile),
                                        ('spec_type','Table'),
                                        ('nphot',nphotons),('ranseed',runSeed),('verbose','T')])
        elif PhotometricTableParams['method'] == 'Monochromatic':
            pfs_params = OrderedDict([('outputfile',PhotonICFile),
                                        ('repository',self.ramsesDir),
                                        ('snapnum',self.ramsesTimestep),
                                        ('star_dom_type',StellarEmissionDomain['star_dom_type']),
                                        ('star_dom_pos',"%.16e %.16e %.16e"%(sdp[0],sdp[1],sdp[2])),
                                        ('star_dom_rsp','%.16e'%(StellarEmissionDomain['star_dom_rsp'])),
                                        ('weight_type','Mono'),
                                        ('weight_input_file',PhotTableFile),
                                        ('spec_type','Mono'),
                                        ('nphot',nphotons),('ranseed',runSeed),('verbose','T')])

        parameters = OrderedDict([('PhotonsFromStars',pfs_params),('ramses',ramses_params)])
        sU.write_parameter_file(pfsParamFile,parameters)

        # tibo
        fff = "%s.pbscript"%(pfsParamFile)
        f = open(fff,'w')
        f.write("#PBS -S /bin/bash \n")
        f.write("#PBS -j oe   \n")
        f.write("#PBS -l select=1:ncpus=1 \n")
        f.write("#PBS -l place=scatter:excl   \n")
        f.write("#PBS -l walltime=24:00:00   \n")
        f.write("#PBS -q workq \n")
        f.write("cd $PBS_O_WORKDIR   \n")
        f.write("source /usr/share/modules/init/bash   \n")
        f.write("module load intel-icc-14/14.0.3.174   \n")
        f.write("module load intel-ifc-14/14.0.3.174   \n")
        f.write("module load intel-mkl-11/11.1.3   \n")
        f.write("module load mpt/2.10   \n")
        f.write("set -x   \n")
        #f.write("/scratch/garel/rascas_sphinx/f90/PhotonsFromStars %s >& %s.log \n"%(pfsParamFile,pfsParamFile))
        f.write("%sPhotonsFromStars %s >& %s.log \n"%(self.rascas_f90,pfsParamFile,pfsParamFile))

        f.close()

        
        return runNum, PhotonICFile


    def __setup_LyaPhotonsFromGas(self,seedStep,EmissionDomain,nphotons,ramses_params,doRecombs,doColls,tcool_resolution=10.):

        import os
        from collections import OrderedDict
        import setupUtils as sU  # read/write parameter files in RASCAS format

        if (doRecombs and doColls) or ((not doRecombs) and (not doColls)):
            raise NameError("ERROR : __setup_LyaPhotonsFromGas needs exactly one of doRecombs and  doColls to be True")
        
        # determine new file number and seed for random number generator
        OK = False
        runNum  = 1
        minseed = 0
        while not OK:
            lpfgParamFile = "%s/%3.3i.LPFG.conf"%(self.surveyDir,runNum)
            if os.path.exists(lpfgParamFile):
                f = open(lpfgParamFile,'r')
                lines = f.readlines()
                for l in lines:
                    if 'ranseed' in l:
                        seed = int(l.split()[2])
                        if seed < minseed:
                            minseed = seed
                        break
                runNum = runNum + 1
            else:
                OK = True
        runSeed = minseed - seedStep
        
        # define PhotonsFromStars parameters
        PhotonICFile = '%s/%3.3i.IC'%(self.surveyDir,runNum)
        edp=EmissionDomain['emission_dom_pos']
        if doRecombs:
            PhotonICFileR = PhotonICFile
            PhotonICFileC = ''
        else:
            PhotonICFileR = ''
            PhotonICFileC = PhotonICFile
            
        lpfg_params  = OrderedDict([('outputfileRec',PhotonICFileR),('outputfileCol',PhotonICFileC),\
                                        ('repository',self.ramsesDir),('snapnum',self.ramsesTimestep),\
                                        ('emission_dom_type',EmissionDomain['emission_dom_type']),\
                                        ('emission_dom_pos',"%.16e %.16e %.16e"%(edp[0],edp[1],edp[2])),\
                                        ('emission_dom_rsp','%.16e'%EmissionDomain['emission_dom_rsp']),\
                                        ('nphotons',nphotons),\
                                        ('ranseed',runSeed),('doRecombs',doRecombs),('doColls',doColls),\
                                        ('tcool_resolution',"%.16e"%(tcool_resolution)),('verbose','T')])
        parameters = OrderedDict([('LyaPhotonsFromGas',lpfg_params),('ramses',ramses_params)])
        sU.write_parameter_file(lpfgParamFile,parameters)

        return runNum, PhotonICFile


        
        
    def __setup_rascas(self,runNum,PhotonICFile,\
                           gas_composition_params,HI_model_params,dust_model_params):

        """ 

        setup a rascas run ... 

        """
        from collections import OrderedDict
        import setupUtils as sU  # read/write parameter files in RASCAS format

        RascasOutputFile = '%s/%3.3i.res'%(self.surveyDir,runNum)
        RascasBackupFile = '%s/%3.3i.bak'%(self.surveyDir,runNum)
        rascas_params = OrderedDict([('DomDumpDir',self.DomDumpDir),('DomDumpFile','domain_decomposition_params.dat'),
                                         ('PhotonICFile',PhotonICFile),
                                         ('fileout',RascasOutputFile),
                                         ('nbuffer',10),('verbose','T')])
        worker_params = OrderedDict([('verbose','F')])
        master_params = OrderedDict([('verbose','F'),('restart','F'),
                                         ('PhotonBakFile',RascasBackupFile),
                                         ('dt_backup',18000.)])
        parameters = OrderedDict([('RASCAS',rascas_params),('worker',worker_params),('master',master_params),
                                      ('gas_composition',gas_composition_params),('HI',HI_model_params),
                                      ('dust',dust_model_params)])
        rascasParamFile = "%s/%3.3i.RASCAS.conf"%(self.surveyDir,runNum)
        sU.write_parameter_file(rascasParamFile,parameters)

         #tibo
        fff = "%s.pbscript"%(rascasParamFile)
        f = open(fff,'w')
        f.write("#PBS -S /bin/bash \n")
        f.write("#PBS -j oe   \n")
        f.write("#PBS -l select=4:ncpus=16:mpiprocs=16 \n")
        f.write("#PBS -l place=scatter:excl   \n")
        f.write("#PBS -l walltime=24:00:00   \n")
        f.write("#PBS -q workq \n")
        f.write("cd $PBS_O_WORKDIR   \n")
        f.write("source /usr/share/modules/init/bash   \n")
        f.write("module load intel-icc-14/14.0.3.174   \n")
        f.write("module load intel-ifc-14/14.0.3.174   \n")
        f.write("module load intel-mkl-11/11.1.3   \n")
        f.write("module load mpt/2.10   \n")
        f.write("set -x   \n")
        #f.write("time mpiexec_mpt -np 64 dplace -s1 -c0-15 /scratch/garel/rascas_sphinx/f90/rascas %s >& %s.log \n"%(rascasParamFile,rascasParamFile))
        f.write("time mpiexec_mpt -np 64 dplace -s1 -c0-15 %srascas %s >& %s.log \n"%(self.rascas_f90,rascasParamFile,rascasParamFile))

        f.close()
            

        
# Gather a collection of rascas runs (within a same directory, but e.g. different bands)
class RascasRunCollection(object):
    def __init__(self,rascasDir,SurveyFiles,savePhots=False,saveOnlyEscaped=True,maxRuns=5,startRun=1):
        try:
            surveys = {}
            for k in SurveyFiles.keys():
                surveys[k] = RascasRun(rascasDir,SurveyFiles[k],
                                        savePhots=savePhots,
                                        saveOnlyEscaped=saveOnlyEscaped,
                                        maxRuns=maxRuns,startRun=startRun)
        except:
            raise NameError('-> ERROR in RascasRunCollection constructor.')

        self.surveys = surveys
        
