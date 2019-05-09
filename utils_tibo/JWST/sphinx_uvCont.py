
################################################################################
import matplotlib
matplotlib.use('Agg')

import numpy as np
from collections import OrderedDict
import os
import sys

#from minirats.HaloFinder.py import haloCatalog as hC
import nircam
import rascasRun as RS
from minirats.HaloFinder.py.haloCatalog import haloCatalog as hC
from minirats.HaloFinder.py.galaxyCatalog import galaxyCatalog as gC


HI_model_params   = OrderedDict([('isotropic','F'),('recoil','T')])
dust_model_params = OrderedDict([('albedo','%.16e'%(0.32)),('g_dust','%.16e'%(0.73)),('dust_model','SMC')])

gas_composition_params  = OrderedDict([ ('f_ion','%.16e'%(0.01)),('Zref','%.16e'%(0.005)),
                                     ('gas_overwrite','F'),('verbose','T')])


################################################################################
# RAMSES-SIMULATION STUFF 
################################################################################
# RAMSES OUTPUT
#ramsesDir      = '/scratch/blaizot/sphinx/05_F1000/02_IC20_BP/' # 183
ramsesDir      = '/scratch/garel/SPHINX_outputs/05_F1000/02_IC20_BP_test_HaloFinder/'
#ramsesTimestep = [] 147 # 183 
list_timesteps = [183,147] #,173,164,155,147,133,121,101,85,74,63]

# CreateDomDump extra parameters
CreateDomDumpOptions = OrderedDict([('reading_method','hilbert'),('verbose','T')])

ramses_params = OrderedDict([('self_shielding','F'),('ramses_rt','T'),('verbose','T'),
                                 ('use_initial_mass','T'),('cosmo','T'),('use_proper_time','T'),
                                 ('read_rt_variables','F')])


for its in range(len(list_timesteps)):
    
    ramsesTimestep = list_timesteps[its]

    ################################################################################
    # RAScas dir. 
    rascas_directory = '/scratch/garel/rascas_sphinx/output/05_F1000/02_IC20_BP_test_HaloFinder/GF_rho_1000_alphap1_npart100/'
    rascas_f90       = '/scratch/garel/rascas_sphinx/f90/'
    runcmd = 'mkdir '+rascas_directory
    os.system(runcmd)
    
    gcat = gC(ramsesDir,ramsesTimestep,HaloDir='GF_rho_1000_alphap1_npart100/',load=True)
    gcat.add_subhalo_mass()
    g_mstar, g_ages, g_metallicities= gcat.get_starprops()
    mstar = g_mstar*1e11
    mstar     = mstar[:gcat.nhalos]   

    redshift = gcat.info['redshift']

    haloid_list = [] #np.array([8,8,8,8,8,8,8,8,8,8])
    print('Halo list')
    print(haloid_list)

    # Compute cell size in code units for domain decompositon
    ###############################################################################
    
    ncell          = 512.0
    coarse_cell_MaxLength_comoving_cu = 3.**(0.5) * 1.0 / 512.0

    ###############################################################################

    print(' ith timestep = ',its)

    for i in range(len(mstar)):
        print(' ')
        print('=============')
        print('igal = ',i, ' / ',len(mstar))

        haloid_list.append(gcat.hnum[:gcat.nhalos][i])
        
        xh = gcat.x_cu[:gcat.nhalos][i]
        yh = gcat.y_cu[:gcat.nhalos][i]
        zh = gcat.z_cu[:gcat.nhalos][i]
        rh = gcat.rvir_cu[:gcat.nhalos][i]
    
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

        # Use rh*factor_rvir + coarser cell max length (in comoving code units) to ensure small haloes can access neighbor cells info
        factor_rvir        = 1.0 
        decomp_dom_rsp     = rh*factor_rvir + coarse_cell_MaxLength_comoving_cu
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
        
        # RASCAS DIRECTORY 
        rascasHaloDir  = '%s/%5.5i/halo%i'%(rascas_directory,ramsesTimestep,gcat.hnum[:gcat.nhalos][i])
        DomDumpDir = 'CDD_HI_dust'   # directory inside rascasHaloDir to contain all CDD outputs.
        
        # PHOTOMETRY parameters
        sedDir       = '/home/garel/seds/'
        photTableDir = '%s/photTables'%(rascas_directory)
        sedModel     = 'bpass100'
        spec_type    = 'Monochromatic'
        nphot        = 200000
        dust_model   = 'SMC'
    
        # define surveys
        # get filters
        surveyName = ['1500A_rf']
        lambda_model = np.array([1500.])
        albedo_model = np.array([0.38])
        g_dust_model = np.array([0.70])
        
        lambda_model = np.array([1500.])
        albedo_model = np.array([0.38])
        g_dust_model = np.array([0.70])
        
        for j in range(len(surveyName)):
            print('---> ',surveyName[j])
            a = RS.RascasSurvey(surveyName[j],rascasHaloDir,DomDumpDir,ramsesDir,ramsesTimestep,rascas_f90)
            PhotometricTableParams=OrderedDict([('sedDir',sedDir),
                                                ('sedModel',sedModel),
                                                ('lbda0_Angstrom',lambda_model[j]),
                                                ('photTableDir',photTableDir),
                                                ('method',spec_type),
                                                ])
            # get pivot wavelength of the filter and interpolate albedo & g_dust (constant over the whole filter)
            albedo = albedo_model[j]
            g_dust = g_dust_model[j]
            
            dust_model_params = OrderedDict([('albedo','%.16e'%(albedo)),('g_dust','%.16e'%(g_dust)),('dust_model',dust_model)])
            
            a.setup_broad_band_survey(PhotometricTableParams,ComputationalDomain,DomainDecomposition,StellarEmissionDomain,
                                        nphotons=nphot,ramses_params=ramses_params,gas_composition_params=gas_composition_params,
                                        HI_model_params=HI_model_params,dust_model_params=dust_model_params)

    # TIBO - write haloid list file
    fff = "%s/%5.5i/haloid_list.dat"%(rascas_directory,ramsesTimestep)
    f = open(fff,'w')
    f.write("Halo IDs \n")
    for j in range(len(haloid_list)):
        f.write("%i \n"%(haloid_list[j]))
    f.close()
    #OBIT

    print(' ')
    print('===============================')
    print(' ')


################################################################################

