# ===================================
# Compiler Configuration
# ===================================

# Default to gfortran, allow for an override (e.g. make FC=ifort)
FC = gfortran

# Default flags
FFLAGS ?= -O3 -cpp -fPIC

# Directory configuration
SRC_DIR := src
TEST_DIR := tests
BUILD_DIR := build

# Module directory flags:
# Gfortran uses -J to specify where to put/find .mod files
# Intel (ifort) uses -module. We default to Gfortran syntax here.
# You can override this: make MOD_FLAG=-module
MOD_FLAG ?= -J

# Combine flags to include the build dir for .mod search
FCFLAGS := $(FFLAGS) $(MOD_FLAG)$(BUILD_DIR) -I$(BUILD_DIR)

# ===================================
# Source and Object Definitions
# ===================================

# Tell make to look for source files in these directories
VPATH = $(SRC_DIR):$(TEST_DIR)

# The list of programs to build (executables)
PROGS = simple lesssimple autosps spec_bin

# The common object files required by the programs
# We wrap them in addprefix to place them inside the build directory
COMMON_NAMES = sps_vars.o sps_utils.o compsp.o csp_gen.o ssp_gen.o \
    getmags.o locate.o funcint.o sps_setup.o pz_convol.o \
    get_tuniv.o intsfwght.o imf.o imf_weight.o add_dust.o \
    getspec.o sbf.o add_bs.o mod_hb.o add_remnants.o getindx.o \
    smoothspec.o mod_gb.o add_nebular.o add_xrb.o write_isochrone.o \
    sfhstat.o linterp.o tsum.o add_agb_dust.o linterparr.o \
    ztinterp.o vacairconv.o igm_absorb.o get_lumdist.o attn_curve.o \
    sfh_weight.o sfhlimit.o sfhinfo.o setup_tabular_sfh.o agn_dust.o

COMMON_OBJS = $(addprefix $(BUILD_DIR)/, $(COMMON_NAMES))

# ===================================
# Rules
# ===================================

.PHONY: all clean test

all: $(PROGS)

# --- Compilation Rules ---

# Ensure build directory exists before compiling
$(BUILD_DIR):
	@mkdir -p $(BUILD_DIR)

# Pattern rule: Compile any .f90 found in VPATH to .o in BUILD_DIR
$(BUILD_DIR)/%.o: %.f90 | $(BUILD_DIR)
	$(FC) $(FCFLAGS) -c $< -o $@

# Specific dependencies to enforce compilation order

# sps_utils.o specifically depends on sps_vars.o
$(BUILD_DIR)/sps_utils.o: $(BUILD_DIR)/sps_vars.o

# All other common objects depend on both vars and utils.
REST_OF_COMMON = $(filter-out $(BUILD_DIR)/sps_vars.o $(BUILD_DIR)/sps_utils.o, $(COMMON_OBJS))

$(REST_OF_COMMON): $(BUILD_DIR)/sps_vars.o $(BUILD_DIR)/sps_utils.o

# Main program objects also wait for modules
$(BUILD_DIR)/simple.o $(BUILD_DIR)/lesssimple.o $(BUILD_DIR)/autosps.o $(BUILD_DIR)/spec_bin.o: $(BUILD_DIR)/sps_vars.o $(BUILD_DIR)/sps_utils.o

# Dependencies for test objects
$(BUILD_DIR)/generate_test_data.o: $(BUILD_DIR)/sps_vars.o $(BUILD_DIR)/sps_utils.o
$(BUILD_DIR)/test_runner.o: $(BUILD_DIR)/sps_vars.o $(BUILD_DIR)/sps_utils.o

# --- Linking Rules ---

autosps: $(BUILD_DIR)/autosps.o $(COMMON_OBJS)
	$(FC) $(FCFLAGS) -o $@ $^

simple: $(BUILD_DIR)/simple.o $(COMMON_OBJS)
	$(FC) $(FCFLAGS) -o $@ $^

lesssimple: $(BUILD_DIR)/lesssimple.o $(COMMON_OBJS)
	$(FC) $(FCFLAGS) -o $@ $^

spec_bin: $(BUILD_DIR)/spec_bin.o $(BUILD_DIR)/sps_vars.o
	$(FC) $(FCFLAGS) -o $@ $^

# --- Test Targets ---

generate_test_data: $(BUILD_DIR)/generate_test_data.o $(COMMON_OBJS)
	$(FC) $(FCFLAGS) -o $@ $^

test_runner: $(BUILD_DIR)/test_runner.o $(COMMON_OBJS)
	$(FC) $(FCFLAGS) -o $@ $^

test: test_runner

# --- Utilities ---

clean:
	rm -rf $(BUILD_DIR) $(PROGS) generate_test_data test_runner

