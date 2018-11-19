
import os

################################################################################

#halo_ids = ['11518','12763','18546','23233','24695','41361','47123','49444','49905','6114']
#timestep = '00183'

#path2confFiles = '/scratch/garel/rascas_sphinx/output/test_Mstar_gt1e-3/'

halo_ids = []
ramsesTimestep = 183

rascas_directory = '/scratch/garel/rascas_sphinx/output/sfr_gt_1e-4'

fff = "%s/%5.5i/haloid_list.dat"%(rascas_directory,ramsesTimestep)
f = open(fff,'r')
header1 = f.readline()

for line in f:
    p = line.split()
    halo_ids.append(p[0])

f.close()


for j in range(len(halo_ids)):
    # CDD 
    cdd_dir = "%s/%5.5i%s%s%s"%(rascas_directory,ramsesTimestep,'/halo',halo_ids[j],'/CDD_HI_dust/')
    cdd_cmd = "%s%s%s"%('qsub ',cdd_dir,'CDD_params.conf.pbscript')
    print(cdd_cmd)
    os.system(cdd_cmd)
   

################################################################################

