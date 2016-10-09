#############################################################################
# Fortran compiler options and directives

#SYSTYPE="CCF"
#SYSTYPE="lyocral"
SYSTYPE="mac"
#SYSTYPE="Pmac"


OPTIONS = -cpp
#OPTIONS += -DSWITCH_OFF_UPARALLEL
#OPTIONS += -DDEBUG


ifeq ($(SYSTYPE),"lyocral")
F90         = ifort
FFLAGS      = -g -traceback -fpp -check all -debug -warn all -ftrapuv 
##FFLAGS      = -O3 -fpp -ftz  -xavx -sox -ipo -fno-alias -fno-fnalias -no-prec-div -no-prec-sqrt -align all
##FFLAGS      = -O3 -fpp -ftz -ipo -fno-alias -fno-fnalias -no-prec-div -no-prec-sqrt -align all
##FFLAGS      = -p -g
FFLAGS += $(OPTIONS)
LDFLAGS     =  
endif


ifeq ($(SYSTYPE),"CCF")
F90         = ifort
#FFLAGS      = -g -fpp -check all -debug -warn all -ftrapuv -traceback
##FFLAGS      = -O3 -fpp -ftz  -xavx -sox -ipo -fno-alias -fno-fnalias -no-prec-div -no-prec-sqrt -align all
FFLAGS      = -O3 -fpp -ftz  -xavx -ipo -fno-alias -fno-fnalias -no-prec-div -no-prec-sqrt -align all
##FFLAGS      = -p -g
FFLAGS += $(OPTIONS)
LDFLAGS     =  -lmpi
endif

ifeq ($(SYSTYPE),"mac")
F90         = gfortran-mp-4.7 -ffree-line-length-none -ffree-form  ###-fopenmp
FFLAGS      = -g -fcheck=all -ftrapv -fbacktrace -Wall
FFLAGS      = -O2
FFLAGS += $(OPTIONS)
###FFLAGS      = -O3
LDFLAGS     = 
endif


ifeq ($(SYSTYPE),"Pmac")
F90         = mpif90-openmpi-gcc48  # mpif90-openmpi-mp
FFLAGS      = -g -fcheck=all -ffree-line-length-none -ffree-form -ftrapv -fbounds-check -fbacktrace
#FFLAGS      = -O3
FFLAGS += $(OPTIONS)
LDFLAGS     = #-lmpi
endif


#############################################################################
# All objects

COMPOBJS = module_dust_model.o module_D_model.o module_HI_model.o module_gas_composition_HI_D_dust.o
#COMPOBJS = module_HI_model.o module_gas_composition_HI.o
#COMPOBJS = module_gas_ramses.o

#GASOBJS = module_constants.o module_ramses.o module_random.o module_utils.o module_uparallel.o $(COMPOBJS)
GASOBJS = module_domain.o module_constants.o module_ramses.o module_random.o module_utils.o module_uparallel.o $(COMPOBJS)

#OBJSCOM = $(GASOBJS) module_select.o module_domain.o module_mesh.o
OBJSCOM = $(GASOBJS) module_select.o module_mesh.o

OBJS1 = $(OBJSCOM) output2subvol.o

OBJS2 = $(OBJSCOM) CreateDomDump.o

#OBJSMC = $(GASOBJS) module_domain.o module_mesh.o module_photon.o
OBJSMC = $(GASOBJS) module_mesh.o module_photon.o 

OBJSMPI = $(OBJSMC) module_parallel_mpi.o module_worker.o module_master.o MCLya.o

#############################################################################
#.SUFFIXES: .o .f90
#.f90.o:
#	$(F90) -c  $(FFLAGS) $<

%.o:%.f90
	$(F90) $(FFLAGS) -c $^ -o $@

#############################################################################
all: clean CreateDomDump PhotonsMonoPointSource serial PhotonsFromStars ## MCLya

PhotonsMonoPointSource: $(OBJSMC) PhotonsMonoPointSource.o
	$(F90) $(FFLAGS) $(LDFLAGS) -o $@ $(OBJSMC) PhotonsMonoPointSource.o

PhotonsFromStars: $(OBJSMC) PhotonsFromStars.o
	$(F90) $(FFLAGS) $(LDFLAGS) -o $@ $(OBJSMC) PhotonsFromStars.o

serial: $(OBJSMC) serial.o
	$(F90) $(FFLAGS) $(LDFLAGS) -o $@ $(OBJSMC) serial.o

MCLya: $(OBJSMPI) MCLya.o
	$(F90) $(FFLAGS) $(LDFLAGS) -o $@ $(OBJSMPI)

output2subvol: $(OBJS1) output2subvol.o
	$(F90) $(FFLAGS) $(LDFLAGS) -o $@ $(OBJS1)

CreateDomDump: $(OBJS2) CreateDomDump.o
	$(F90) $(FFLAGS) $(LDFLAGS) -o $@ $(OBJS2)

clean:
	rm -f $(OBJS1) $(OBJS2) $(OBJSMC) $(OBJSMPI) *.o *.mod
#############################################################################
