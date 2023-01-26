#!/bin/bash

# Install script to 1) install Julia via conda, and 2) install the Julia dependencies for this project

set -x

# Define some useful directories
pge_dir=$( cd "$(dirname "$0")" ; pwd -P )
app_dir=$(dirname ${pge_dir})

# Get SpectralUnmixing repo
cd $app_dir
git clone https://github.com/emit-sds/SpectralUnmixing.git -b v0.2.1
specun_dir="$app_dir/SpectralUnmixing"

# Install Julia and then install Julia dependencies
conda create -n spectral-unmixing -y -c conda-forge julia=1.7 python=3.9 gdal
source activate spectral-unmixing
pip install Pillow

pushd $specun_dir
export JULIA_SSL_CA_ROOTS_PATH=""
julia -e 'using Pkg; Pkg.activate("."); Pkg.add(path="https://github.com/kmsquire/ArgParse2.jl"); Pkg.instantiate()'
export JULIA_PROJECT=$specun_dir
popd
