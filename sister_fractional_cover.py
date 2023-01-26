#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SISTER
Space-based Imaging Spectroscopy and Thermal PathfindER
Author: Winston Olson-Duvall
"""

import json
import os
import subprocess
import shutil
import sys
try:
    from osgeo import gdal
except:
    import gdal
import numpy as np
from PIL import Image


def get_frcov_basename(corfl_basename, crid):
    # Replace product type
    tmp_basename = corfl_basename.replace("L2A_CORFL", "L2A_FRCOV")
    # Split, remove old CRID, and add new one
    tokens = tmp_basename.split("_")[:-1] + [crid]
    return "_".join(tokens)


def generate_metadata(run_config, frcov_met_json_path):
    # Create .met.json file from runconfig for fractional cover
    metadata = run_config["metadata"]
    metadata["product"] = "FRCOV"
    metadata["processing_level"] = "L2A"
    metadata["description"] = "Fractional cover (soil, vegetation, water, snow)"
    with open(frcov_met_json_path, "w") as f:
        json.dump(metadata, f, indent=4)

def main():
    """
        This function takes as input the path to an inputs.json file and exports a run config json
        containing the arguments needed to run the SISTER ISOFIT PGE.

    """
    in_file = sys.argv[1]

    # Read in runconfig
    print("Reading in runconfig")
    with open(in_file, "r") as f:
        run_config = json.load(f)

    # Make work dir
    if not os.path.exists("work"):
        print("Making work directory")
        os.mkdir('work')

    # Make output dir
    if not os.path.exists("output"):
        print("Making output directory")
        os.mkdir('output')

    # Define paths and variables
    sister_frcov_dir = os.path.abspath(os.path.dirname(__file__))
    specun_dir = os.path.join(os.path.dirname(sister_frcov_dir), "SpectralUnmixing")

    corfl_basename = None
    for file in run_config["inputs"]["file"]:
        if "corrected_reflectance_dataset" in file:
            corfl_basename = os.path.basename(file["corrected_reflectance_dataset"])
    frcov_basename = get_frcov_basename(corfl_basename, run_config["inputs"]["config"]["crid"])

    corfl_img_path = f"work/{corfl_basename}"
    corfl_hdr_path = f"work/{corfl_basename}.hdr"
    frcov_img_path = f"work/{frcov_basename}"

    # Copy the input files into the work directory (don't use .bin)
    shutil.copyfile(f"input/{corfl_basename}/{corfl_basename}.bin",corfl_img_path)
    shutil.copyfile(f"input/{corfl_basename}/{corfl_basename}.hdr",corfl_hdr_path)

    # Build command and run unmix.jl
    unmix_exe = f"{specun_dir}/unmix.jl"
    endmember_lib_path = f"{sister_frcov_dir}/data/veg_soil_water_snow_endmembers.csv"
    log_path = f"output/{frcov_basename}.log"
    cmd = ["julia"]

    # Add parallelization if n_cores > 1
    if run_config["inputs"]["config"]["n_cores"] > 1:
        cmd += ["-p", str(run_config["inputs"]["config"]["n_cores"])]
    # Add required args
    cmd += [
        unmix_exe,
        corfl_img_path,
        endmember_lib_path,
        "class",
        frcov_img_path,
        "--mode=sma",
        f"--log_file={log_path}"
    ]
    # Add the optional args
    if run_config["inputs"]["config"]["refl_nodata"] != "None":
        cmd += [f"--refl_nodata={run_config['inputs']['config']['refl_nodata']}"]
    if run_config["inputs"]["config"]["refl_scale"] != "None":
        cmd += [f"--refl_scale={run_config['inputs']['config']['refl_scale']}"]
    if run_config["inputs"]["config"]["normalization"] != "None":
        cmd += [f"--normalization={run_config['inputs']['config']['normalization']}"]

    print("Running unmix.jl command: " + " ".join(cmd))
    subprocess.run(" ".join(cmd), shell=True)

    # Generate metadata in .met.json file
    frcov_met_json_path = f"output/{frcov_basename}.met.json"
    print(f"Generating metadata from runconfig to {frcov_met_json_path}")
    generate_metadata(run_config, frcov_met_json_path)

    # Generate quicklook
    frcov_ql_path = f"output/{frcov_basename}.png"
    print(f"Generating quicklook to {frcov_ql_path}")

    frcov_img_file = f'{frcov_img_path}_fractional_cover'
    frcov_gdal = gdal.Open(frcov_img_file)

    bands = []

    for band_num in range(1,4):
        band = frcov_gdal.GetRasterBand(band_num)
        bands.append(band.ReadAsArray())

    no_data = band.GetNoDataValue()
    rgb=  np.array(bands)
    rgb[rgb == no_data] = np.nan

    rgb = np.moveaxis(rgb,0,-1).astype(float)
    bottom = np.nanpercentile(rgb,5,axis = (0,1))
    top = np.nanpercentile(rgb,95,axis = (0,1))
    rgb = np.clip(rgb,bottom,top)
    rgb = (rgb-np.nanmin(rgb,axis=(0,1)))/(np.nanmax(rgb,axis= (0,1))-np.nanmin(rgb,axis= (0,1)))
    rgb = (rgb*255).astype(np.uint8)

    im = Image.fromarray(rgb)
    im.save(frcov_ql_path)

    #Convert to COG
    temp_file =  f'{frcov_img_file}.tif'
    out_file =  f"output/{frcov_basename}.tif"

    print(f"Creating COG {out_file}")

    band_names = ['soil',
                  'vegetation',
                  'water',
                  'snow/ice',
                  'brightness']

    units = ['PERCENT',
                  'PERCENT',
                  'PERCENT',
                  'PERCENT',
                  'UNITELESS']

    descriptions=  ['SOIL PERCENT COVER',
                     'VEGETATION PERCENT COVER',
                     'WATER PERCENTCOVER',
                     'SNOW/ICE PERCENT COVER',
                     'BRIGHTNESS']

    # Set the output raster transform and projection properties
    driver = gdal.GetDriverByName("GTIFF")
    tiff = driver.Create(temp_file,
                         frcov_gdal.RasterXSize,
                         frcov_gdal.RasterYSize,
                         frcov_gdal.RasterCount,
                         gdal.GDT_Float32)

    tiff.SetGeoTransform(frcov_gdal.GetGeoTransform())
    tiff.SetProjection(frcov_gdal.GetProjection())
    tiff.SetMetadataItem("DESCRIPTION","FRACTIONAL COVER")

    # Write bands to file
    for i,band_name in enumerate(band_names):
        in_band = frcov_gdal.GetRasterBand(i+1).ReadAsArray()
        out_band = tiff.GetRasterBand(i+1)
        out_band.WriteArray(in_band)
        out_band.SetDescription(band_name)
        out_band.SetNoDataValue(no_data)
        out_band.SetMetadataItem("UNITS",units[i])
        out_band.SetMetadataItem("DESCRIPTION",descriptions[i])
    del tiff, driver

    os.system(f"gdaladdo -minsize 900 {temp_file}")
    os.system(f"gdal_translate {temp_file} {out_file} -co COMPRESS=LZW -co TILED=YES -co COPY_SRC_OVERVIEWS=YES")

    out_runconfig =  f"output/{frcov_basename}.runconfig.json"
    print(f"Copying runconfig to {out_runconfig}")
    shutil.copyfile(in_file,
                    out_runconfig)


if __name__ == "__main__":
    main()
