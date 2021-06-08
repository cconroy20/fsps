# NB: this code is fairly memory intensive.  If you run into memory errors
# when building this code, you may have to increase the default memory 
# allocation.

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
# PGI
#F90 = pgf90
# PGI (franklin.nersc.gov)
#F90 = ftn
#---------------------
# Pathscale
#F90 = pathf90
#---------------------
# IBM
#F90 = xlf90_r
#---------------------

#========================
# Compiler Optimizations
#  (Uncomment only one)
#========================
F90FLAGS := -O3 -cpp -fPIC $(F90FLAGS) # prepend so a user can set -Ofast from outside the Makefile
#---------------------
# Intel (Itanium 2 processor)
#F90FLAGS = -O3 -mcpu=itanium2 -cpp
# Intel (EM64T/AMD64 processor)
#F90FLAGS = -O3 -funroll-loops -xP -cpp
#---------------------

PROGS = simple lesssimple autosps spec_bin

COMMON = sps_vars.o sps_utils.o compsp.o csp_gen.o ssp_gen.o \
	getmags.o locate.o funcint.o sps_setup.o pz_convol.o \
	get_tuniv.o intsfwght.o imf.o imf_weight.o add_dust.o \
	getspec.o sbf.o add_bs.o mod_hb.o add_remnants.o getindx.o \
	smoothspec.o mod_gb.o add_nebular.o write_isochrone.o \
	sfhstat.o linterp.o tsum.o add_agb_dust.o linterparr.o \
	ztinterp.o vacairconv.o igm_absorb.o get_lumdist.o attn_curve.o \
	sfh_weight.o sfhlimit.o sfhinfo.o setup_tabular_sfh.o agn_dust.o

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
	$(F90) $(F90FLAGS) sps_utils.f90 -c

%.o : %.f90 sps_vars.o sps_utils.o
	$(F90) $(F90FLAGS) -o $@ -c $<

