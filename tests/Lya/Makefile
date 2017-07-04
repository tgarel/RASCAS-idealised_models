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
F90         = gfortran-mp-4.7 ###-fopenmp
FFLAGS      = -g -fcheck=all -ffree-line-length-none -ffree-form -ftrapv -fbacktrace
##FFLAGS      = -O2
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

COMPOBJS = ../../f90/module_dust_model.o ../../f90/module_D_model.o ../../f90/module_HI_model.o
GASOBJS = ../../f90/module_constants.o ../../f90/module_random.o ../../f90/module_utils.o ../../f90/module_uparallel.o $(COMPOBJS)

#############################################################################
#.SUFFIXES: .o .f90
#.f90.o:
#	$(F90) -c  $(FFLAGS) $<

%.o:%.f90
	$(F90) $(FFLAGS) -c $^ -o $@

#############################################################################
all: clean2 tests clean

tests: $(GASOBJS) tests.o 
	$(F90) $(FFLAGS) -o $@ $(GASOBJS) tests.o

clean:
	rm -f $(GASOBJS) *.o *.mod

clean2:
	rm -f $(GASOBJS) *.o *.mod
	rm -f *.o
	rm -f *.mod
	rm -f *~
	rm -f tests 	
#############################################################################
