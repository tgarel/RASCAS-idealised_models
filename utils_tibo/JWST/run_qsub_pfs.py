
import os

################################################################################

#halo_ids = ['11518','12763','18546','23233','24695','41361','47123','49444','49905','6114']
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
    
print(halo_ids)

for j in range(len(halo_ids)):
    # PFS 
    pfs_dir = "%s/%5.5i%s%s%s"%(rascas_directory,ramsesTimestep,'/halo',halo_ids[j],'/1500A_rf/')
    pfs_cmd = "%s%s%s"%('qsub ',pfs_dir,'001.PFS.conf.pbscript')
    print(pfs_cmd)
    os.system(pfs_cmd)
  

################################################################################

