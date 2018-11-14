
################################################################################

from minirats.HaloFinder.py import haloCatalog as hC
import rascasRun as RS
import numpy as np
from collections import OrderedDict
import nircam

################################################################################

HI_model_params   = OrderedDict([('isotropic','F'),('recoil','T')])
dust_model_params = OrderedDict([('albedo','%.16e'%(0.32)),('g_dust','%.16e'%(0.73)),('dust_model','SMC')])

gas_composition_params  = OrderedDict([ ('f_ion','%.16e'%(0.01)),('Zref','%.16e'%(0.005)),
                                     ('gas_overwrite','F'),('verbose','T')])

################################################################################
# RAMSES-SIMULATION STUFF 
################################################################################
# RAMSES OUTPUT
ramsesDir      = '/cral2/sphinx/05_F1000/02_IC20_BP/'
ramsesTimestep = 183 

# CreateDomDump extra parameters
CreateDomDumpOptions = OrderedDict([('reading_method','hilbert'),('verbose','T')])

ramses_params = OrderedDict([('self_shielding','T'),('ramses_rt','T'),('verbose','F'),
                                 ('use_initial_mass','T'),('cosmo','T'),('use_proper_time','T'),
                                 ('read_rt_variables','T')])

################################################################################

hcat = hC.haloCatalog(ramsesDir,ramsesTimestep,zoom=False)
hcat.load_catalog()
hcat.convert_distances()
mstar = hcat.get_Mstar()

redshift = hcat.info['redshift']

ids = np.where(mstar > 1.e-3)

for i in range(len(mstar[ids])):
    print(' ')
    print('=============')
    print('halo %i'%hcat.hnum[ids][i])
    print('Mstar = %.8e'%(mstar[ids][i]*1.e11))
    print('coordinates = ', (hcat.x_cu[ids][i]), (hcat.y_cu[ids][i]), (hcat.z_cu[ids][i]))
    print('Rvir = ',hcat.rvir_cu[ids][i])

    xh = hcat.x_cu[ids][i]
    yh = hcat.y_cu[ids][i]
    zh = hcat.z_cu[ids][i]
    rh = hcat.rvir_cu[ids][i]
    
    ################################################################################
    ### HALO-dependent STUFF ###
    ################################################################################
    
    # computational domain 
    comput_dom_type = 'sphere' 
    comput_dom_pos  = xh, yh, zh
    comput_dom_rsp  = rh
    ComputationalDomain = OrderedDict([('comput_dom_type',comput_dom_type),
                                        ('comput_dom_pos',comput_dom_pos),
                                        ('comput_dom_rsp',comput_dom_rsp)])
    # domain decomposition
    decomp_dom_type    = 'sphere'
    decomp_dom_ndomain =     1 
    decomp_dom_xc      = xh
    decomp_dom_yc      = yh
    decomp_dom_zc      = zh
    decomp_dom_rsp     = rh*1.10
    DomainDecomposition = OrderedDict([('decomp_dom_type',decomp_dom_type),
                        ('decomp_dom_ndomain',decomp_dom_ndomain),('decomp_dom_xc',decomp_dom_xc),
                        ('decomp_dom_yc',decomp_dom_yc),('decomp_dom_zc',decomp_dom_zc),
                        ('decomp_dom_rsp',decomp_dom_rsp)])
    
    # stellar emission domain
    star_dom_type = 'sphere' 
    star_dom_pos  = xh, yh, zh
    star_dom_rsp  = rh
    StellarEmissionDomain = OrderedDict([('star_dom_type',star_dom_type),
                                        ('star_dom_pos',star_dom_pos),
                                        ('star_dom_rsp',star_dom_rsp)])
    
    ################################################################################
    ################################################################################
    # RASCAS PARAMETERS 
    ################################################################################
    HI_model_params   = OrderedDict([('isotropic','F'),('recoil','T')])
    #dust_model_params = OrderedDict([('albedo','%.16e'%(0.32)),('g_dust','%.16e'%(0.73)),('dust_model','SMC')])
    gas_composition_params  = OrderedDict([ ('f_ion','%.16e'%(0.01)),('Zref','%.16e'%(0.005)),
                                        ('gas_overwrite','F'),('verbose','T')])
    
    # RASCAS DIRECTORY 
    rascasDir  = 'halo%i'%hcat.hnum[ids][i]
    DomDumpDir = 'CDD_HI_dust'   # directory inside rascasDir to contain all CDD outputs.
    
    # PHOTOMETRY parameters
    sedDir = '/cral/leo/seds/'
    photTableDir = 'photTables'
    
    # define NIRCam surveys
    # get NIRCam filters
    surveyName = ('F070W', 'F090W', 'F115W', 'F150W', 'F200W', 'F277W', 'F356W', 'F444W')

    # according to Li&Draine Table6
    lambda_model = np.array([912., 1000., 1216., 1500., 2200., 3000., 3650. , 4400., 5500., 7000., 9000., 12200., 16300., 22000., 34500., 36000.])
    albedo_model = np.array([0.24, 0.27, 0.32, 0.38, 0.42, 0.58, 0.62, 0.65, 0.67, 0.66, 0.63, 0.58, 0.51, 0.43, 0.28, 0.26])
    g_dust_model = np.array([0.73, 0.72, 0.73, 0.70, 0.56, 0.57, 0.58, 0.57, 0.54, 0.48, 0.40, 0.29, 0.21, 0.13, 0.005, -0.004])

    for j in range(len(surveyName)):
        print('---> ',surveyName[j])
        a = RS.RascasSurvey(surveyName[j],rascasDir,DomDumpDir,ramsesDir,ramsesTimestep)
        
        # get NIRCam filter according to surveyName
        f = nircam.nircamFilter(surveyName[j])
        #
        if surveyName[j] == 'F444W':
            f.load(threshold=1.e-2)
            print('threshold used = 0.01')
        lambda_min = f.lambda_min / (1.+redshift)
        lambda_max = f.lambda_max / (1.+redshift)
        
        if (lambda_min < 1215.):
            print('need Lyman-alpha for this filter',surveyName[j])
        else:
            PhotometricTableParams=OrderedDict([('sedDir',sedDir),
                                                ('sedModel','bpass100'),
                                                ('lbdamin_Angstrom',lambda_min),
                                                ('lbdamax_Angstrom',lambda_max),
                                                ('photTableDir',photTableDir),
                                                ('method','TabulatedSpectra'),
                                                ('dlambda',0.1)])
            # get pivot wavelength of the filter and interpolate albedo & g_dust (constant over the whole filter)
            albedo = np.interp(f.lambda_mean, lambda_model, albedo_model)
            g_dust = np.interp(f.lambda_mean, lambda_model, g_dust_model)
            
            dust_model_params = OrderedDict([('albedo','%.16e'%(albedo)),('g_dust','%.16e'%(g_dust)),('dust_model','SMC')])
            
            a.setup_broad_band_survey(PhotometricTableParams,ComputationalDomain,DomainDecomposition,StellarEmissionDomain,
                                          nphotons=10000000,ramses_params=ramses_params,gas_composition_params=gas_composition_params,
                                          HI_model_params=HI_model_params,dust_model_params=dust_model_params)

################################################################################

