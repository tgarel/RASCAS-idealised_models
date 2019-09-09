# example use of the class rascasSurvey (from rascasRun.py) (and some functions in mocks.py).
from collections import OrderedDict
from rascasRun import RascasSurvey 
import mocks as mo

#### define some paths
## Which SED model we use and where it is: 
sedModel = 'bpass100'
sedDir   = '/Users/blaizot/Documents/Astro/seds/%s/'%(sedModel)
## Where are RAMSES outputs
ramDir = '/Volumes/SIMULATIONS/P13-20h-1Mpc-MUSIC/Zoom-7-10508/SPHINX_run/'
## Where RASCAS will run (see note below)
RascasRootDir = "%s/%s/"%(ramDir,'RASCAS_test')



#### define RASCAS parameters describing the general properties of the simulation to process
## indexes of useful variables in RAMSES outputs 
ipre    = 5
imet    = 6
ixhii   = 7
ixheii  = 8
ixheiii = 9
## ramses parameters:
ramses_params    = OrderedDict([('self_shielding','T'),('ramses_rt','T'),('verbose','F'),\
                                    ('use_initial_mass','T'),('cosmo','T'),('use_proper_time','T'),\
                                    ('read_rt_variables','T'),('QuadHilbert','T'),\
                                    ('itemp',"%i"%ipre),('imetal',"%i"%imet),('ihii',"%i"%ixhii),('iheii',"%i"%ixheii),('iheiii',"%i"%ixheiii),\
                                    ('U_precision','8'),('RT_precision','8')])
                                    
#### define parameters for RASCAS gas mix: The mix itself, and the parameters for each scatterer. 
DomDumpDir = 'HI_D_Dust'   # This is the name of the directory where the domain dumps will be stored (NB: relative to rascasDir defined below). 
Deuterium_params = OrderedDict([('isotropic','F'),('recoil','T')])
HI_model_params  = OrderedDict([('isotropic','F'),('recoil','T'),('HI_core_skip','T'),('xcritmax','%.10e'%(100.0))])
dust_params      = OrderedDict([('albedo','%.16e'%(0.32)),('g_dust','%.16e'%(0.73)),('dust_model','SMC')])
gas_comp_params  = OrderedDict([('deut2H_nb_ratio','%.16e'%(3.000E-05)),('f_ion','%.16e'%(0.01)), \
                                   ('Zref','%.16e'%(0.005)),('gas_overwrite','F'),('verbose','T')])


#### What is below could be within a loop on halos or timesteps or both ... 
#### -> the general strtucture is to have one directory "rascasDir" _per target_
#### Then, within this directory, a number of sub-directories will be created to contain surveys.
#### ---- NB: i generally create a RASCAS directory next to the outputs of ramses, and then a sub-dir per target
#### ---- (e.g. for rascas runs of halo 101 of timestep 203 of a simulation stored at /some/path/to/simulation/,
#### ---- I would define the RASCAS directory as /some/path/to/simulation/RASCAS/203/101/ and all the RASCAS stuff will be there.)
#### ---- NB: directories are created by rascasRun.py -> no need to create them before hand. 
my_halo_num = 101
timestep    = 164 
rascasDir   = "%s/%5.5i/%5.5i/"%(RascasRootDir,my_halo_num,timestep)

## define all domains with same center and extent (should be replaced by coords of halos,...)
CENTER, RADIUS = [0.5,0.5,0.5], 0.12
REDSHIFT = 3.0  # should depend on timestep for cosmo runs 

# use trivial functions to convert center/radius into domains 
ComputationalDomain = mo.define_computational_domain(CENTER,RADIUS)
DomainDecomposition = mo.define_domain_decomposition(CENTER,RADIUS)
EmissionDomain      = mo.define_emission_domain(CENTER,RADIUS)

# define properties of mock observations (here to produce a data-cube as with MUSE). 
Lya_mock_params = mo.define_mock_params(CENTER,mockSpec=False,mockImage=False,mockCube=True,\
                                        cube_image_npix=200,cube_lbda_npix=100,cube_lbdmin=1190,cube_lbdmax=1240,cubeSide_cu=2.*RADIUS)


# Lya from recombinations (with core-skipping with xcritmax=100 -> "CS100"),
# and out to the virial radius (-> that would depend on the value of RADIUS defined above). 
surveyName = 'RecLya_CS100_Rvir'
r = RascasSurvey(surveyName,rascasDir,DomDumpDir,ramDir,timestep)
r.setup_Lya_survey(ComputationalDomain, DomainDecomposition, EmissionDomain,\
                       nphotons = 100000, seedStep=100, doRecombs=True, doColls=False,\
                       ramses_params=ramses_params, gas_composition_params=gas_comp_params,\
                       HI_model_params=HI_model_params, dust_model_params=dust_params,\
                       mock_params=Lya_mock_params,reading_method='fullbox_omp')

# Lya from collisions
surveyName = 'ColLya_CS100_Rvir'
r = RascasSurvey(surveyName,rascasDir,DomDumpDir,ramDir,timestep)    
r.setup_Lya_survey(ComputationalDomain, DomainDecomposition, EmissionDomain,\
                       nphotons = 100000, seedStep=100, doRecombs=False, doColls=True,\
                       ramses_params=ramses_params, gas_composition_params=gas_comp_params,\
                       HI_model_params=HI_model_params, dust_model_params=dust_params,\
                       mock_params=Lya_mock_params,reading_method='fullbox_omp')

# Continuum around Lya
surveyName = '1190-1240A_Rvir'
r = RascasSurvey(surveyName,rascasDir,DomDumpDir,ramDir,timestep)    
PhotonsFromStarsParams = mo.define_photonsFromStarsParams(sedDir,spec_type='Table',spec_table_lmin_Ang=1190.,spec_table_lmax_Ang=1240.)
StellarEmissionDomain  = mo.define_stellar_emission_domain(ComputationalDomain['comput_dom_pos'],ComputationalDomain['comput_dom_rsp'])
r.setup_broad_band_survey(PhotonsFromStarsParams, ComputationalDomain, DomainDecomposition, StellarEmissionDomain,\
                              nphotons = 1000000, seedStep=100,\
                              ramses_params=ramses_params, gas_composition_params=gas_comp_params,\
                              dust_model_params=dust_params,reading_method='fullbox_omp',\
                              mock_params=Lya_mock_params,HI_model_params=HI_model_params)
 




