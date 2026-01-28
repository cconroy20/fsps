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
# Format: "Legacy_File_Suffix|Runtime_Args"
declare -a configurations=(
    "MILES-1_MIST-1|--isoc mist --spec miles"
    "MILES-0_MIST-0_BASEL-1_PADOVA-1|--isoc pdva --spec basel"
    "MILES-1_MIST-1_THEMIS-1_DL07-0|--isoc mist --spec miles --dust themis"
    "MILES-0_C3K-1_MIST-1|--isoc mist --spec c3k_afe+0.0"
    "MIST-0_BPASS-1|--isoc bpss --spec bpass"
)

# Create the data directory if it doesn't exist
mkdir -p data

# Move into src directory to run Make
cd ../src

echo "=========================================================="
echo "Compiling FSPS objects and generator..."
echo "=========================================================="

# Clean previous build artifacts
make clean > /dev/null 2>&1

# Build the generator target ONCE
make generate_test_data F90FLAGS="-cpp -fPIC -O3"

if [ ! -f ./generate_test_data ]; then
    echo "ERROR: Compilation failed."
    exit 1
fi

# Main Loop
for config in "${configurations[@]}"; do
    IFS="|" read -r suffix args <<< "$config"
    
    echo "=========================================================="
    echo "Running configuration: $args"
    echo "Output: tests/data/sps_ref_${suffix}.bin"
    echo "=========================================================="

    # Run the generator
    ./generate_test_data $args
    
    # Move and rename output to the tests/data folder
    if [ -f "sps_test_output.bin" ]; then
        mv sps_test_output.bin "../tests/data/sps_ref_${suffix}.bin"
        echo "Created: tests/data/sps_ref_${suffix}.bin"
    else
        echo "ERROR: Output file not generated for $args"
        exit 1
    fi
    
    echo ""
done

# Cleanup the executable from src
rm -f generate_test_data
echo "Done."
