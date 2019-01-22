#############################################################################
# Fortran compiler options and directives

#SYSTYPE="CCF"
#SYSTYPE="lyocral"
#SYSTYPE="lyocralP"
SYSTYPE="macP"
#SYSTYPE="mac"
#SYSTYPE="ccage"

OPTIONS = -cpp


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
F90         = gfortran -ffree-line-length-none -ffree-form -fopenmp
FFLAGS      = -g -fcheck=all -ftrapv -fbacktrace -fbounds-check
#FFLAGS      = -O2
FFLAGS += $(OPTIONS)
LDFLAGS     = 
endif


ifeq ($(SYSTYPE),"macP")
#F90         = mpif90-openmpi-mp
F90	= mpif90 -ffree-line-length-none -ffree-form 
#FFLAGS      = -g -fcheck=all -ftrapv -fbounds-check -fbacktrace
FFLAGS      = -O3
FFLAGS += $(OPTIONS)
LDFLAGS     = -lmpi
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
# 5/ SiII 1190/1193
#COMPOBJS = module_SiII_1190_model.o module_SiII_1193_model.o module_gas_composition_SiII_1190_1193.o
# 6/ SiII 1190/1193 and dust
#COMPOBJS = module_dust_model.o module_SiII_1190_model.o module_SiII_1193_model.o module_gas_composition_SiII_1190_1193_dust.o
# 7/ MgII
#COMPOBJS = module_MgII_2796_model.o module_MgII_2804_model.o module_gas_composition_MgII.o
# 8/ MgII and dust
#COMPOBJS = module_dust_model.o module_MgII_2796_model.o module_MgII_2804_model.o module_gas_composition_MgII_dust.o
# 9/ SiII_1260 and dust
COMPOBJS = module_dust_model.o module_SiII_1260_model.o module_gas_composition_SiII_1260_dust.o
# 10/ dust only
#COMPOBJS = module_dust_model.o module_gas_composition_dust.o
# 11/ HI and HII
#COMPOBJS = module_HI_model.o module_gas_composition_HI_HII.o
# 12/ OI_1302 and dust
#COMPOBJS = module_dust_model.o module_OI_1302_model.o module_gas_composition_OI1302_dust.o


GASOBJS = module_constants.o module_random.o module_utils.o module_voigt.o module_domain.o coolrates_module.o module_ramses.o module_uparallel.o $(COMPOBJS)

OBJSDOM = $(GASOBJS) module_select.o module_mesh.o

OBJSPHO = $(GASOBJS) module_mesh.o module_mock.o module_photon.o module_spectra.o module_CD.o #module_opt_ray.o

# GATHER -- 
#OBJSMPI = $(OBJSPHO) module_parallel_mpi.o module_worker.o module_master.o
OBJSMPI = module_parallel_mpi.o $(OBJSPHO) module_worker.o module_master.o
# -- GATHER

#############################################################################
#.SUFFIXES: .o .f90
#.f90.o:
#	$(F90) -c  $(FFLAGS) $<

%.o:%.f90
	$(F90) $(FFLAGS) -c $^ -o $@

#############################################################################
all: clean CreateDomDump PhotonsFromStars rascas

PhotonsFromSourceModel: $(OBJSPHO) PhotonsFromSourceModel.o
	$(F90) $(FFLAGS) $(LDFLAGS) -o $@ $(OBJSPHO) PhotonsFromSourceModel.o

PhotonsFromStars: $(OBJSPHO) PhotonsFromStars.o
	$(F90) $(FFLAGS) $(LDFLAGS) -o $@ $(OBJSPHO) PhotonsFromStars.o

star_luminosities: $(OBJSPHO) star_luminosities.o
	$(F90) $(FFLAGS) $(LDFLAGS) -o $@ $(OBJSPHO) star_luminosities.o

CreateDomDump: $(OBJSDOM) CreateDomDump.o
	$(F90) $(FFLAGS) $(LDFLAGS) -o $@ $(OBJSDOM) CreateDomDump.o

ExtractSubvol: $(OBJSDOM) ExtractSubvol.o
	$(F90) $(FFLAGS) $(LDFLAGS) -o $@ $(OBJSDOM) ExtractSubvol.o

rascas-serial: $(OBJSPHO) rascas-serial.o
	$(F90) $(FFLAGS) $(LDFLAGS) -o $@ $(OBJSPHO) rascas-serial.o

rascas: $(OBJSMPI) rascas.o
	$(F90) $(FFLAGS) $(LDFLAGS) -o $@ $(OBJSMPI) rascas.o

optical_bin: $(OBJSPHO) optical_bin.o
	$(F90) $(FFLAGS) $(LDFLAGS) -o $@ $(OBJSPHO) optical_bin.o

IC_ColumnDensity: $(OBJSPHO) IC_ColumnDensity.o
	$(F90) $(FFLAGS) $(LDFLAGS) -o $@ $(OBJSPHO) IC_ColumnDensity.o

ColumnDensity: $(OBJSPHO) ColumnDensity.o
	$(F90) $(FFLAGS) $(LDFLAGS) -o $@ $(OBJSPHO) ColumnDensity.o

clean:
	rm -f $(OBJSDOM) $(OBJSPHO) $(OBJSMPI) *.o *.mod
#############################################################################
