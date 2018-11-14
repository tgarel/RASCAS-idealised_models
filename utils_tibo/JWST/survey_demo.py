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
dust_model_params = OrderedDict([('albedo','%.16e'%(0.32)),('g_dust','%.16e'%(0.73)),('dust_model','SMC')])
gas_composition_params  = OrderedDict([ ('f_ion','%.16e'%(0.01)),('Zref','%.16e'%(0.005)),
                                     ('gas_overwrite','F'),('verbose','T')])



# RASCAS DIRECTORY 
rascasDir  = '/cral2/blaizot/RASCAS_test/'
DomDumpDir = 'CDD_HI_dust'   # directory inside rascasDir to contain all CDD outputs.
surveyName = 'mySurvey'


# PHOTOMETRY parameters 
#sedDir = '/Users/blaizot/Documents/Astro/seds/'
sedDir = '/home/cral/blaizot/seds/'
photTableDir = rascasDir
PhotometricTableParams=OrderedDict([('sedDir',sedDir),
                                         ('sedModel','bc16_Basel'),
                                         ('lbdamin_Angstrom',2000),
                                         ('lbdamax_Angstrom',2200),
                                         ('photTableDir',photTableDir),
                                         ('method','TabulatedSpectra'),
                                         ('dlambda',0.1)])

a = RS.RascasSurvey(surveyName,rascasDir,DomDumpDir,ramsesDir,ramsesTimestep)
a.setup_broad_band_survey(PhotometricTableParams,ComputationalDomain,DomainDecomposition,StellarEmissionDomain,
                              nphotons=1000000,ramses_params=ramses_params,gas_composition_params=gas_composition_params,
                              HI_model_params=HI_model_params,dust_model_params=dust_model_params)


# setup a Lya run
# emission domain
emission_dom_type = 'sphere' 
emission_dom_pos  = 6.0334515571594238e-01, 2.2293025255203247e-01, 1.4043170213699341e-01 
emission_dom_rsp  = 1.1605447158217430e-02 
EmissionDomain = OrderedDict([('emission_dom_type',emission_dom_type),
                                       ('emission_dom_pos',emission_dom_pos),
                                       ('emission_dom_rsp',emission_dom_rsp)])

surveyName = 'LyaCol'
a = RS.RascasSurvey(surveyName,rascasDir,DomDumpDir,ramsesDir,ramsesTimestep)
a.setup_Lya_survey(ComputationalDomain, DomainDecomposition, EmissionDomain,\
                       nphotons = 10000, seedStep=100, doRecombs=False, doColls=True, ramses_params=ramses_params, \
                       gas_composition_params=gas_composition_params, HI_model_params=HI_model_params, \
                       dust_model_params=dust_model_params)

surveyName = 'LyaRec'
a = RS.RascasSurvey(surveyName,rascasDir,DomDumpDir,ramsesDir,ramsesTimestep)
a.setup_Lya_survey(ComputationalDomain, DomainDecomposition, EmissionDomain,\
                       nphotons = 10000, seedStep=100, doRecombs=True, doColls=False, ramses_params=ramses_params, \
                       gas_composition_params=gas_composition_params, HI_model_params=HI_model_params, \
                       dust_model_params=dust_model_params)



