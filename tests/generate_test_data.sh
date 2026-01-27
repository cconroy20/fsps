#!/bin/bash

# Generate FSPS test reference data for multiple library configurations.
# This script must be run from the tests/ directory.
# Usage: ./generate_test_data.sh

set -e

# Validate we're in the correct directory
if [[ ! -f "generate_test_data.sh" ]]; then
    echo "ERROR: This script must be run from the tests/ directory"
    exit 1
fi

# Define the library combinations
declare -a configurations=(
    "-DMILES=1 -DMIST=1"
    "-DMILES=0 -DMIST=0 -DBASEL=1 -DPADOVA=1"
    "-DMILES=1 -DMIST=1 -DTHEMIS=1 -DDL07=0"
    "-DMILES=0 -DC3K=1 -DMIST=1"
    "-DMIST=0 -DBPASS=1"
)

# Create the data directory if it doesn't exist
mkdir -p data

# Move into src directory to run Make
cd ../src

# Main Loop
for flags in "${configurations[@]}"; do
    echo "=========================================================="
    echo "Running configuration: $flags"
    echo "=========================================================="

    # Clean previous build artifacts
    make clean > /dev/null 2>&1

    # Build ONLY the generator target with specific flags
    # This automatically builds the necessary library objects (sps_vars.o, etc.)
    # We pass the flags explicitly to override the Makefile defaults
    echo "Compiling FSPS objects..."
    make generate_test_data F90FLAGS="-cpp -fPIC -O3 $flags"

    # Run the generator
    if [ -f ./generate_test_data ]; then
        echo "Running generator..."
        ./generate_test_data
        
        # Move and rename output to the tests/data folder
        safe_name=$(echo "$flags" | sed 's/-D//g' | sed 's/ /_/g' | sed 's/=/-/g')
        mv sps_test_output.bin "../tests/data/sps_ref_${safe_name}.bin"
        
        echo "Created: tests/data/sps_ref_${safe_name}.bin"
    else
        echo "ERROR: Compilation failed for $flags"
        exit 1
    fi
    
    echo ""
done

# Cleanup the executable from src
rm -f generate_test_data
echo "Done."