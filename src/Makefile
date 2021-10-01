
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
F90FLAGS = -O -cpp #-march=native -mcmodel=large
#-fopenmp -Wl,-stack_size,4000000 
#---------------------
# Intel 
#F90FLAGS = -O3 -assume byterecl
# Intel (EM64T/AMD64 processor)
#F90FLAGS = -O3 -funroll-loops -xP
#---------------------

PROGS = runsps simple lesssimple autosps metalml \
	spec_bin bytype colorteff bpass_bin

COMMON =  nrtype.o sps_vars.o sps_utils.o nr.o \
	nrutil.o compsp.o csp_gen.o ssp_gen.o getmags.o gasdev.o \
	ran_state.o locate.o qromb.o polint.o trapzd.o sps_setup.o \
	add_agb_dust.o ran1.o intsfwght.o ztinterp.o \
	get_tuniv.o imf.o imf_weight.o linterparr.o add_dust.o getspec.o  \
	mod_hb.o add_remnants.o getindx.o smoothspec.o mod_gb.o \
	write_isochrone.o sfhstat.o linterp.o add_bs.o \
	tsum.o add_nebular.o sbf.o vacairconv.o igm_absorb.o get_lumdist.o \
	sfh_weight.o sfhlimit.o sfhinfo.o setup_tabular_sfh.o attn_curve.o \
	agn_dust.o funcint.o

all : $(PROGS)

clean :
	rm -rf *.o *.mod *.MOD *~

autosps : autosps.o $(COMMON)
	$(F90) -o autosps.exe autosps.o $(COMMON)

bytype : bytype.o $(COMMON)
	$(F90) -o bytype.exe bytype.o $(COMMON)

spec_bin : spec_bin.o $(COMMON)
	$(F90) -o spec_bin.exe spec_bin.o $(COMMON)

bpass_bin : bpass_bin.o $(COMMON)
	$(F90) -o bpass_bin.exe bpass_bin.o $(COMMON)

runsps : runsps.o $(COMMON)
	$(F90) -o runsps.exe runsps.o $(COMMON)

colorteff : colorteff.o $(COMMON)
	$(F90) -o colorteff.exe colorteff.o $(COMMON)

simple : simple.o $(COMMON)
	$(F90) -o simple.exe simple.o $(COMMON)

lesssimple : lesssimple.o $(COMMON)
	$(F90) -o lesssimple.exe lesssimple.o $(COMMON)

gencat : gencat.o $(COMMON)
	$(F90) -o gencat.exe gencat.o $(COMMON)

metalml : metalml.o $(COMMON)
	$(F90) -o metalml.exe metalml.o $(COMMON)

#------------------------------------------------------------------------#

autosps.o : sps_vars.o nrtype.o nr.o sps_utils.o

bytype.o : sps_vars.o nrtype.o nr.o sps_utils.o

spec_bin.o : sps_vars.o nrtype.o nr.o sps_utils.o

bpass_bin.o : sps_vars.o nrtype.o nr.o sps_utils.o

colorteff.o : sps_vars.o nrtype.o sps_utils.o

runsps.o : sps_vars.o nrtype.o sps_utils.o

simple.o : sps_vars.o nrtype.o sps_utils.o

lesssimple.o : sps_vars.o nrtype.o sps_utils.o 

gencat.o : sps_vars.o nrtype.o sps_utils.o

metalml.o : sps_vars.o nrtype.o sps_utils.o

sps_utils.o : nrutil.o nrtype.o nr.o sps_vars.o

ran_state.o: nrutil.o

%.o : %.f90
	$(F90) $(F90FLAGS) -o  $@ -c $<

%.o : nr/%.f90
	$(F90) $(F90FLAGS) -o $@ -c $<


