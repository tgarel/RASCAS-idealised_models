
import os

################################################################################

halo_ids = ['11518','12763','18546','23233','24695','41361','47123','49444','49905','6114']

path2confFiles = '/scratch/garel/test_rascas_sphinx/00183/'

# Clean CDD files if needed

for j in range(len(halo_ids)):
    # CDD 
    cdd_dir = "%s%s%s%s"%(path2confFiles,'halo',halo_ids[j],'/CDD_HI_dust/')  
    cdd_cmd = "%s%s%s"%('rm ',cdd_dir,'/CDD*')  
    os.system(cdd_cmd)
    
  

################################################################################

