
import os

################################################################################


halo_ids = []
ramsesTimestep = 183

rascas_directory = '/scratch/garel/rascas_sphinx/output/sfr_gt_1e-4_corrBuffer_domain'

fff = "%s/%5.5i/haloid_list.dat"%(rascas_directory,ramsesTimestep)
f = open(fff,'r')
header1 = f.readline()

for line in f:
    p = line.split()
    halo_ids.append(p[0])

f.close()

# Clean CDD files if needed

for j in range(len(halo_ids)):
    # CDD 
    cdd_dir = "%s/%5.5i%s%s%s"%(rascas_directory,ramsesTimestep,'/halo',halo_ids[j],'/CDD_HI_dust/')
    cdd_cmd = "%s%s%s"%('rm ',cdd_dir,'/CDD*')  
    os.system(cdd_cmd)
    
  

################################################################################

