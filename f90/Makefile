#############################################################################
# Fortran compiler options and directives

#SYSTYPE="CCF"
#SYSTYPE="lyocral"
#SYSTYPE="lyocralP"
SYSTYPE="mac"
#SYSTYPE="macP"
#SYSTYPE="ccage"

OPTIONS = -cpp
#OPTIONS += -DSWITCH_OFF_UPARALLEL
#OPTIONS += -DDEBUG


ifeq ($(SYSTYPE),"lyocral")
F90         = ifort
#FFLAGS      = -g -traceback -fpp -check all -debug -warn all -ftrapuv 
##FFLAGS      = -O3 -fpp -ftz  -xavx -sox -ipo -fno-alias -fno-fnalias -no-prec-div -no-prec-sqrt -align all
FFLAGS      = -O3 -fpp -ftz -ipo -fno-alias -fno-fnalias -no-prec-div -no-prec-sqrt -align all
##FFLAGS      = -p -g
FFLAGS += $(OPTIONS)
LDFLAGS     = -lmpi
endif

ifeq ($(SYSTYPE),"lyocralP")
F90         = mpif90 -ffree-line-length-none -ffree-form
FFLAGS      = #-g -fcheck=all -ftrapv -fbacktrace -Wall
FFLAGS      = -O2
FFLAGS += $(OPTIONS)
LDFLAGS     = -lmpi
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
#F90         = gfortran-mp-4.7 -ffree-line-length-none -ffree-form  ###-fopenmp
F90         = gfortran-mp-6 -ffree-line-length-none -ffree-form  ###-fopenmp
#FFLAGS      = -g -fcheck=all -ftrapv -fbacktrace -Wall -fbounds-check
FFLAGS      = -O2
FFLAGS += $(OPTIONS)
###FFLAGS      = -O3
LDFLAGS     = 
endif


ifeq ($(SYSTYPE),"macP")
F90         = mpif90-openmpi-mp
FFLAGS      = -g -fcheck=all -ffree-line-length-none -ffree-form -ftrapv -fbounds-check -fbacktrace -Wall
#FFLAGS      = -O3
FFLAGS += $(OPTIONS)
LDFLAGS     = #-lmpi
endif


ifeq ($(SYSTYPE),"ccage")
F90         = mpif90
#FFLAGS      = -g -fcheck=all -ffree-line-length-none -ffree-form -ftrapv -fbounds-check -fbacktrace
FFLAGS      = -O3 -ffree-line-length-none -ffree-form
FFLAGS += $(OPTIONS)
LDFLAGS     = 
endif


#############################################################################
# All objects

# choose the gas mix you want to use
# 1/ HI, D, and dust
#COMPOBJS = module_dust_model.o module_D_model.o module_HI_model.o module_gas_composition_HI_D_dust.o
# 2/ HI and dust
#COMPOBJS = module_dust_model.o module_HI_model.o module_gas_composition_HI_dust.o
# 3/ pure HI
#COMPOBJS = module_HI_model.o module_gas_composition_HI.o
# 4/ pure RAMSES hydro variables
#COMPOBJS = module_gas_composition_ramses.o
# 5/ SiII
COMPOBJS = module_SiII_1190_model.o module_SiII_1193_model.o module_gas_composition_SiII.o


GASOBJS = module_domain.o module_constants.o module_ramses.o module_random.o module_utils.o module_uparallel.o $(COMPOBJS)

OBJSDOM = $(GASOBJS) module_select.o module_mesh.o

OBJSPHO = $(GASOBJS) module_mesh.o module_photon.o 

OBJSMPI = $(OBJSPHO) module_parallel_mpi.o module_worker.o module_master.o

#############################################################################
#.SUFFIXES: .o .f90
#.f90.o:
#	$(F90) -c  $(FFLAGS) $<

%.o:%.f90
	$(F90) $(FFLAGS) -c $^ -o $@

#############################################################################
all: clean CreateDomDump PhotonsFromSourceModel PhotonsFromStars rascas-serial

PhotonsFromSourceModel: $(OBJSPHO) PhotonsFromSourceModel.o
	$(F90) $(FFLAGS) $(LDFLAGS) -o $@ $(OBJSPHO) PhotonsFromSourceModel.o

PhotonsFromStars: $(OBJSPHO) PhotonsFromStars.o
	$(F90) $(FFLAGS) $(LDFLAGS) -o $@ $(OBJSPHO) PhotonsFromStars.o

CreateDomDump: $(OBJSDOM) CreateDomDump.o
	$(F90) $(FFLAGS) $(LDFLAGS) -o $@ $(OBJSDOM) CreateDomDump.o

ExtractSubvol: $(OBJSDOM) ExtractSubvol.o
	$(F90) $(FFLAGS) $(LDFLAGS) -o $@ $(OBJSDOM) ExtractSubvol.o

rascas-serial: $(OBJSPHO) rascas-serial.o
	$(F90) $(FFLAGS) $(LDFLAGS) -o $@ $(OBJSPHO) rascas-serial.o

rascas: $(OBJSMPI) rascas.o
	$(F90) $(FFLAGS) $(LDFLAGS) -o $@ $(OBJSMPI) rascas.o

clean:
	rm -f $(OBJSDOM) $(OBJSPHO) $(OBJSMPI) *.o *.mod
#############################################################################
