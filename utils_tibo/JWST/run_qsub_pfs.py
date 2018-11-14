
import os

################################################################################

halo_ids = ['11518','12763','18546','23233','24695','41361','47123','49444','49905','6114']
timestep = '00183'

path2confFiles = '/scratch/garel/rascas_sphinx/output/test_Mstar_gt1e-3/'
    
for j in range(len(halo_ids)):
    # PFS 
    pfs_dir = "%s%s%s%s%s"%(path2confFiles,timestep,'/halo',halo_ids[j],'/1500A_rf/')
    pfs_cmd = "%s%s%s"%('qsub ',pfs_dir,'001.PFS.conf.pbscript')
    print(pfs_cmd)
    os.system(pfs_cmd)
  

################################################################################

