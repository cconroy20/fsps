
#======================
# Fortran 90 compiler
# (Uncomment only one)
#======================
# GNU
F90 = gfortran
#---------------------
# Intel
#F90 = ifort 
#---------------------

#========================
# Compiler Optimizations
#  (Uncomment only one)
#========================
# GNU
F90FLAGS := -O3 -cpp -fPIC $(F90FLAGS) # prepend so a user can set -Ofast from outside the Makefile
#-fopenmp -Wl,-stack_size,4000000
#---------------------
# Intel
#F90FLAGS = -O3 -assume byterecl
# Intel (EM64T/AMD64 processor)
#F90FLAGS = -O3 -funroll-loops -xP
#---------------------

PROGS = sps_vars.o sps_utils.o simple

COMMON = sps_vars.o sps_utils.o sps_setup.o compsp.o \
    csp_gen.o ssp_gen.o ztinterp.o \
	mod_hb.o mod_gb.o \
	add_nebular.o add_bs.o add_xrb.o add_remnants.o \
	add_dust.o add_agb_dust.o agn_dust.o attn_curve.o \
	get_tuniv.o imf.o imf_weight.o \
	sfhstat.o sfh_weight.o intsfwght.o sfhlimit.o sfhinfo.o setup_tabular_sfh.o \
	getspec.o getmags.o getindx.o sbf.o \
	smoothspec.o vacairconv.o write_isochrone.o \
	igm_absorb.o get_lumdist.o \
	funcint.o locate.o linterp.o tsum.o linterparr.o

all : $(PROGS)

clean :
	rm -rf *.o *.mod *.MOD *~

autosps : autosps.o $(COMMON)
	$(F90) $(F90FLAGS) -o autosps.exe autosps.o $(COMMON)

simple : simple.o $(COMMON)
	$(F90) $(F90FLAGS) -o simple.exe simple.o $(COMMON)

lesssimple : lesssimple.o $(COMMON)
	$(F90) $(F90FLAGS) -o lesssimple.exe lesssimple.o $(COMMON)

spec_bin : spec_bin.o sps_vars.o
	$(F90) $(F90FLAGS) -o spec_bin.exe spec_bin.o sps_vars.o

sps_vars.o : sps_vars.f90
	$(F90) $(F90FLAGS) sps_vars.f90 -c

sps_utils.o : sps_utils.f90 sps_vars.o

%.o : %.f90
	$(F90) $(F90FLAGS) -o  $@ -c $<



