#!/bin/bash

# Install script to 1) install Julia via conda, and 2) install the Julia dependencies for this project

set -x

# Define some useful directories
pge_dir=$( cd "$(dirname "$0")" ; pwd -P )
app_dir=$(dirname ${pge_dir})

cd $app_dir

# Get SpectralUnmixing repo
git clone https://github.com/EnSpec/SpectralUnmixing.git -b v0.2.1-sister
specun_dir="$app_dir/SpectralUnmixing"

# Install Julia and then install Julia dependencies
conda create -n spectral-unmixing -y -c conda-forge julia=1.7 python=3.8 gdal=3.6.0 pandas=2.0.1 awscli
source activate spectral-unmixing
pip install Pillow

# Download snow climatology dataset
aws s3 cp s3://sister-ops-registry/packages/LIN10A1_snow_climatology_13day.tif .

git clone https://github.com/EnSpec/hytools.git
cd hytools
pip install .
cd ..

pushd $specun_dir
export JULIA_SSL_CA_ROOTS_PATH=""
julia -e 'using Pkg; Pkg.activate("."); Pkg.add(path="https://github.com/kmsquire/ArgParse2.jl"); Pkg.instantiate()'
export JULIA_PROJECT=$specun_dir
popd

cd $app_dir
