# demo of usage of rascasRun.

import rascasRun as RS ## anticipating renaming to rascasSurvey

from collections import OrderedDict



HI_model_params   = OrderedDict([('isotropic','F'),('recoil','T')])
dust_model_params = OrderedDict([('albedo','%.16e'%(0.32)),('g_dust','%.16e'%(0.73)),('dust_model','SMC')])


gas_composition_params  = OrderedDict([ ('f_ion','%.16e'%(0.01)),('Zref','%.16e'%(0.005)),
                                     ('gas_overwrite','F'),('verbose','T')])



################################################################################
### HALO-dependent STUFF ###
################################################################################

# computational domain 
comput_dom_type = 'sphere' 
comput_dom_pos  = 6.0334515571594238e-01, 2.2293025255203247e-01, 1.4043170213699341e-01 
comput_dom_rsp  = 1.1605447158217430e-02 
ComputationalDomain = OrderedDict([('comput_dom_type',comput_dom_type),
                                       ('comput_dom_pos',comput_dom_pos),
                                       ('comput_dom_rsp',comput_dom_rsp)])

# domain decomposition
decomp_dom_type    = 'sphere'
decomp_dom_ndomain =     1 
decomp_dom_xc      = 6.0334515571594238e-01 
decomp_dom_yc      = 2.2293025255203247e-01 
decomp_dom_zc      = 1.4043170213699341e-01 
decomp_dom_rsp     = 1.2765991874039173e-02 
DomainDecomposition = OrderedDict([('decomp_dom_type',decomp_dom_type),
                    ('decomp_dom_ndomain',decomp_dom_ndomain),('decomp_dom_xc',decomp_dom_xc),
                    ('decomp_dom_yc',decomp_dom_yc),('decomp_dom_zc',decomp_dom_zc),
                    ('decomp_dom_rsp',decomp_dom_rsp)])

# stellar emission domain
star_dom_type = 'sphere' 
star_dom_pos  = 6.0334515571594238e-01, 2.2293025255203247e-01, 1.4043170213699341e-01 
star_dom_rsp  = 1.1605447158217430e-02 
StellarEmissionDomain = OrderedDict([('star_dom_type',star_dom_type),
                                       ('star_dom_pos',star_dom_pos),
                                       ('star_dom_rsp',star_dom_rsp)])

################################################################################



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


################################################################################
# RASCAS PARAMETERS 
################################################################################
HI_model_params   = OrderedDict([('isotropic','F'),('recoil','T')])
#dust_model_params = OrderedDict([('albedo','%.16e'%(0.32)),('g_dust','%.16e'%(0.73)),('dust_model','SMC')])
gas_composition_params  = OrderedDict([ ('f_ion','%.16e'%(0.01)),('Zref','%.16e'%(0.005)),
                                     ('gas_overwrite','F'),('verbose','T')])



# RASCAS DIRECTORY 
rascasDir  = 'demo-ubvri-survey/'
DomDumpDir = 'CDD_HI_dust'   # directory inside rascasDir to contain all CDD outputs.


# PHOTOMETRY parameters 
#sedDir = '/Users/blaizot/Documents/Astro/seds/'
sedDir = '/cral/leo/seds/'
photTableDir = rascasDir

# define UBVRI surveys
surveyName = ('3600mono', '4400mono', '5500mono', '7000mono', '9000mono')
lambda0    = (3600, 4400, 5500, 7000, 9000)
# according to Li&Draine
albedo     = (0.62, 0.65, 0.67, 0.66, 0.63)
g_dust     = (0.58, 0.57, 0.54, 0.48, 0.40)

for i in range(len(surveyName)):
    print('---> ',surveyName[i])
    a = RS.RascasSurvey(surveyName[i],rascasDir,DomDumpDir,ramsesDir,ramsesTimestep)
    PhotometricTableParams=OrderedDict([('sedDir',sedDir),
                                            ('sedModel','bc16_Basel'),
                                            ('lbda0_Angstrom',lambda0[i]),
                                            ('photTableDir',photTableDir),
                                            ('method','Monochromatic')])
    dust_model_params = OrderedDict([('albedo','%.16e'%(albedo[i])),('g_dust','%.16e'%(g_dust[i])),('dust_model','SMC')])
    a.setup_broad_band_survey(PhotometricTableParams,ComputationalDomain,DomainDecomposition,StellarEmissionDomain,
                                nphotons=11000000,ramses_params=ramses_params,gas_composition_params=gas_composition_params,
                                HI_model_params=HI_model_params,dust_model_params=dust_model_params)



#a.load()



