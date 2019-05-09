
import os

################################################################################

list_timesteps = [183,147] #,173,164,155,147,133,121,101,85,74,63]

rascas_directory = '/scratch/garel/rascas_sphinx/output/05_F1000/02_IC20_BP_test_HaloFinder/GF_rho_1000_alphap1_npart100/'

for its in range(len(list_timesteps)):
    
    ramsesTimestep = list_timesteps[its]

    fff = "%s/%5.5i/haloid_list.dat"%(rascas_directory,ramsesTimestep)
    f = open(fff,'r')
    header1 = f.readline()

    halo_ids = []
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

