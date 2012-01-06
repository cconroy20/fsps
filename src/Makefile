# Modified 2008-10-24 by Ted Kisner <tskisner@lbl.gov>

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
# GNU (version < 4.3.0)
F90FLAGS = -O -cpp
# GNU (version >= 4.3.0)
#F90FLAGS = -O3 -march=native -cpp
#---------------------
# Intel (Itanium 2 processor)
#F90FLAGS = -O3 -mcpu=itanium2 -cpp
# Intel (EM64T/AMD64 processor)
#F90FLAGS = -O3 -funroll-loops -xP -cpp
#---------------------
# PGI
#F90FLAGS = -O3 -fastsse -Mcache_align -cpp
#---------------------

PROGS = simple lesssimple autosps 

COMMON = sps_vars.o nrtype.o sps_utils.o nr.o nrutil.o compsp.o	 \
ssp_gen.o getmags.o ran.o locate.o qromb.o polint.o trapzd.o tridag.o \
sps_setup.o pz_convol.o get_tuniv.o imf.o imf_weight.o add_dust.o getspec.o \
add_bs.o mod_hb.o add_remnants.o getindx.o velbroad.o mod_gb.o \
write_isochrone.o sfhstat.o linterp.o


all : $(PROGS)

clean :
	rm -rf *.o *.mod *.MOD *~

autosps : autosps.o $(COMMON)
	$(F90) -o autosps.exe autosps.o $(COMMON)

simple : simple.o $(COMMON)
	$(F90) -o simple.exe simple.o $(COMMON)

lesssimple : lesssimple.o $(COMMON)
	$(F90) -o lesssimple.exe lesssimple.o $(COMMON)

autosps.o : sps_vars.o nrtype.o nr.o sps_utils.o

simple.o : sps_vars.o nrtype.o sps_utils.o

lesssimple.o : sps_vars.o nrtype.o sps_utils.o

sps_utils.o : nrutil.o nr.o sps_vars.o


%.o : %.f90
	$(F90) $(F90FLAGS) -o $@ -c $<

%.o : nr/%.f90
	$(F90) $(F90FLAGS) -o $@ -c $<


