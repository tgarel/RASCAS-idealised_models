
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
rascas_directory = '/scratch/garel/rascas_sphinx/output/test_Mstar_gt1e-3/'
rascas_f90       = '/scratch/garel/rascas_sphinx/f90/'

################################################################################

# RAMSES OUTPUT
ramsesDir      = '/scratch/blaizot/sphinx/05_F1000/02_IC20_BP/'
ramsesTimestep = 183
DomDumpDir     = 'CDD_HI_dust'

hcat = hC.haloCatalog(ramsesDir,ramsesTimestep,zoom=False)
hcat.load_catalog()
hcat.convert_distances()
mstar = hcat.get_Mstar()

redshift = hcat.info['redshift']

ids = np.where(mstar > 1.e-3)

m1500 = np.array(len(mstar[ids]))
surveyName = ['1500A_rf']

for i in range(len(mstar[ids])):
    print(' ')
    print('=============')
    print('halo %i'%hcat.hnum[ids][i])
    print('Mstar = %.8e'%(mstar[ids][i]*1.e11))
    print('coordinates = ', (hcat.x_cu[ids][i]), (hcat.y_cu[ids][i]), (hcat.z_cu[ids][i]))
    print('Rvir = ',hcat.rvir_cu[ids][i])

    print( hcat.hnum[ids][i])
    if hcat.hnum[ids][i] != 6114:
        rascasDir  = '%s/%5.5i/halo%i'%(rascas_directory,ramsesTimestep,hcat.hnum[ids][i])
        print(rascasDir)
        print(DomDumpDir)
        print(ramsesDir)
        
        lambda_model = np.array([1500.])
        albedo_model = np.array([0.38])
        g_dust_model = np.array([0.70])
        
        
        for j in range(len(surveyName)):
            print('---> ',surveyName[j])
            a = RS.RascasSurvey(surveyName[j],rascasDir,DomDumpDir,ramsesDir,ramsesTimestep,rascas_f90)
            a.load(savePhots=True,maxRuns=1)
            print(a.noResults)
            print(a.abs_mag())
            print(a.nRealPhotons)
           
################################################################################

