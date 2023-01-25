#!/bin/bash

# Description:
#
# The top-level run script to run the SpectralUnmixing fractional cover PGE on SISTER (Space-based Imaging Spectroscopy
# and Thermal PathfindER).
#
# File inputs:
#   - corrected_reflectance_dataset
#
# Config inputs:
#   - n_cores
#   - refl_nodata
#   - refl_scale
#   - normalization
#
# Positional inputs: (None)

# Activate conda environment from docker image
source activate spectral-unmixing

# Get repository directory and export the JULIA_PROJECT env variable
REPO_DIR=$(cd "$(dirname "$0")"; pwd -P)

# Generate runconfig
python ${REPO_DIR}/generate_runconfig.py inputs.json

# Execute isofit
python ${REPO_DIR}/sister_fractional_cover.py runconfig.json