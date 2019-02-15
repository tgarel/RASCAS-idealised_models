
################################################################################
import matplotlib
matplotlib.use('Agg')

from minirats.HaloFinder.py import haloCatalog as hC
import rascasRun as RS
import numpy as np
from collections import OrderedDict
import nircam


################################################################################
# RAScas dir. 
rascas_directory = '/scratch/garel/rascas_sphinx/output/sfr_gt_1e-4_corrBuffer_domain'
rascas_f90       = '/scratch/garel/rascas_sphinx/f90/'

################################################################################

HI_model_params   = OrderedDict([('isotropic','F'),('recoil','T')])
dust_model_params = OrderedDict([('albedo','%.16e'%(0.32)),('g_dust','%.16e'%(0.73)),('dust_model','SMC')])

gas_composition_params  = OrderedDict([ ('f_ion','%.16e'%(0.01)),('Zref','%.16e'%(0.005)),
                                     ('gas_overwrite','F'),('verbose','T')])

################################################################################
# RAMSES-SIMULATION STUFF 
################################################################################
# RAMSES OUTPUT
#ramsesDir      = '/scratch/blaizot/sphinx/05_F1000/02_IC20_BP/' # 183
ramsesDir      = '/scratch/garel/SPHINX_outputs/02_IC20_BP/'
ramsesTimestep = 147 # 183 

# CreateDomDump extra parameters
CreateDomDumpOptions = OrderedDict([('reading_method','hilbert'),('verbose','T')])

ramses_params = OrderedDict([('self_shielding','F'),('ramses_rt','T'),('verbose','T'),
                                 ('use_initial_mass','T'),('cosmo','T'),('use_proper_time','T'),
                                 ('read_rt_variables','F')])

################################################################################

hcat = hC.haloCatalog(ramsesDir,ramsesTimestep,zoom=False)
hcat.load_catalog()
hcat.convert_distances()
mstar = hcat.get_Mstar() # in 1e11 Msun 
sfr100 = hcat.get_halo_SFRs() # Msun/yr
mvir   = hcat.mvir * 1.0e11

redshift = hcat.info['redshift']

# Compute cell size in code units for domain decompositon
###############################################################################

ncell          = 512.0

#boxsize_cmpc   = 10.0
#cell_size_cmpc = boxsize_cmpc / ncell
#coarse_cell_MaxLength_cmpc = 3.**(0.5) * cell_size_cmpc # comoving Mpc
#coarse_cell_MaxLength_pmpc = coarse_cell_MaxLength_cmpc / (1.0+redshift) # physical Mpc
#pc2cm       = 3.086e18
#cm2pc       = 1./pc2cm
#unit_l_pmpc = hcat.info['unit_l'] * cm2pc / 1.e6 # converts comoving code units to phys. Mpc
#coarse_cell_MaxLength_comoving_cu = coarse_cell_MaxLength_pmpc / unit_l_pmpc # comoving code units

# actually, boxsize=1 in cu, so i can do: coarse_cell_MaxLength_comoving_cu = 1 / 512...
coarse_cell_MaxLength_comoving_cu = 3.**(0.5) * 1.0 / 512.0

###############################################################################



#ids = np.where(mstar > 1.e-3) 
ids = np.where(sfr100 > 1.e-4)
#ids = np.where(sfr100 > 1.e-6)
#ids = np.where(mvir > 1.e7)
#ids = np.where((mvir > 1.e7) & (mstar > 0.0e0))

print(len(mstar),len(mvir),len(mstar[ids]))

#TIBO------------
# Create haloid list array
haloid_list = [] #np.array([8,8,8,8,8,8,8,8,8,8])
print('Halo list')
print(haloid_list)
#OBIT------------

for i in range(len(mstar[ids])):
    print(' ')
    print('=============')
    print('halo %i'%hcat.hnum[ids][i])
    # TIBO - Write to haloid list
    print(i)
    haloid_list.append(hcat.hnum[ids][i])
    #haloid_list[i] = hcat.hnum[ids][i]
    # OBIT
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

    # Use rh*1.10 + coarser cell max length (in comoving code units) to ensure small haloes can access neighbor cells info
    decomp_dom_rsp     = rh*1.10 + coarse_cell_MaxLength_comoving_cu
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
    rascasDir  = '%s/%5.5i/halo%i'%(rascas_directory,ramsesTimestep,hcat.hnum[ids][i])
    DomDumpDir = 'CDD_HI_dust'   # directory inside rascasDir to contain all CDD outputs.
    
    # PHOTOMETRY parameters
    sedDir       = '/home/garel/seds/'
    photTableDir = '%s/photTables'%(rascas_directory)
    sedModel     = 'bpass100'
    spec_type    = 'Monochromatic'
    nphot        = 1000000
    dust_model   = 'SMC'

    # define surveys
    # get filters
    surveyName = ['1500A_rf']
    lambda_model = np.array([1500.])
    albedo_model = np.array([0.38])
    g_dust_model = np.array([0.70])

    # according to Li&Draine Table6
    # lambda_model = np.array([912., 1000., 1216., 1500., 2200., 3000., 3650. , 4400., 5500., 7000., 9000., 12200., 16300., 22000., 34500., 36000.])
    # albedo_model = np.array([0.24, 0.27, 0.32, 0.38, 0.42, 0.58, 0.62, 0.65, 0.67, 0.66, 0.63, 0.58, 0.51, 0.43, 0.28, 0.26])
    # g_dust_model = np.array([0.73, 0.72, 0.73, 0.70, 0.56, 0.57, 0.58, 0.57, 0.54, 0.48, 0.40, 0.29, 0.21, 0.13, 0.005, -0.004])

    lambda_model = np.array([1500.])
    albedo_model = np.array([0.38])
    g_dust_model = np.array([0.70])

    for j in range(len(surveyName)):
        print('---> ',surveyName[j])
        a = RS.RascasSurvey(surveyName[j],rascasDir,DomDumpDir,ramsesDir,ramsesTimestep,rascas_f90)
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
 
################################################################################

