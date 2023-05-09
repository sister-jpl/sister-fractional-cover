#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 27 10:43:59 2023

@author: achlus
"""

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
import hytools as ht
import pandas as pd
import numpy as np

try:
    from osgeo import gdal
except:
    import gdal
from PIL import Image


def get_frcov_basename(corfl_basename, crid):
    # Replace product type
    tmp_basename = corfl_basename.replace("L2A_CORFL", "L2B_FRCOV")
    # Split, remove old CRID, and add new one
    tokens = tmp_basename.split("_")[:-1] + [str(crid)]
    return "_".join(tokens)

def generate_metadata(run_config, frcov_met_json_path):
    # Create .met.json file from runconfig for fractional cover
    metadata = run_config["metadata"]
    metadata["product"] = "FRCOV"
    metadata["processing_level"] = "L2B"
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
        if "reflectance_dataset" in file:
            corfl_basename = os.path.basename(file["reflectance_dataset"])
    frcov_basename = get_frcov_basename(corfl_basename, run_config["inputs"]["config"]["crid"])

    corfl_img_path = f"work/{corfl_basename}"
    corfl_hdr_path = f"work/{corfl_basename}.hdr"
    frcov_img_path = f"work/{frcov_basename}"

    # Copy the input files into the work directory (don't use .bin)
    shutil.copyfile(f"input/{corfl_basename}/{corfl_basename}.bin",corfl_img_path)
    shutil.copyfile(f"input/{corfl_basename}/{corfl_basename}.hdr",corfl_hdr_path)

    #Load reflectance im
    rfl = ht.HyTools()
    rfl.read_file(corfl_img_path)

    line_data = (rfl.get_band(0) == rfl.no_data).sum(axis=1)
    start_line = 1+np.argwhere(line_data != rfl.columns)[0][0]
    end_line = rfl.lines - np.argwhere(np.flip(line_data) != rfl.columns)[0][0] -1

    endmember_lib_path = f"{sister_frcov_dir}/data/veg_soil_water_snow_endmembers.csv"
    endmembers = pd.read_csv(endmember_lib_path)

    no_snow = endmembers[endmembers['class'] != 'snow']
    no_snow.to_csv('endmembers_no_snow.csv',index = False)

    no_water = endmembers[endmembers['class'] != 'water']
    no_water.to_csv('endmembers_no_water.csv',index = False)

    # Build command and run unmix.jl
    unmix_exe = f"{specun_dir}/unmix.jl"
    log_path = f"output/{frcov_basename}.log"

    for endmember_file in ['no_snow','no_water']:
        # Add required args
        cmd = ["julia"]

        # Add parallelization if n_cores > 1
        if run_config["inputs"]["config"]["n_cores"] > 1:
            cmd += ["-p", str(run_config["inputs"]["config"]["n_cores"])]

        cmd += [
            unmix_exe,
            corfl_img_path,
            f'endmembers_{endmember_file}.csv',
            "class",
            f"{frcov_img_path}_{endmember_file}",
            "--mode=sma",
            f"--log_file={log_path}",
            "--num_endmembers=3",
            f"--refl_nodata={rfl.no_data}",
            f"--start_line={start_line}",
            f"--end_line={end_line}",
            "--normalization=brightness",
            f"--refl_scale={run_config['inputs']['config']['refl_scale']}"]

        print("Running unmix.jl command: " + " ".join(cmd))
        subprocess.run(" ".join(cmd), shell=True)


    no_snow_frcov_file = f'{frcov_img_path}_no_snow_fractional_cover'
    no_snow_gdal = gdal.Open(no_snow_frcov_file)
    no_snow_frcov  = no_snow_gdal.ReadAsArray()
    no_snow_frcov[no_snow_frcov==rfl.no_data] = np.nan

    no_water_frcov_file = f'{frcov_img_path}_no_water_fractional_cover'
    no_water_gdal = gdal.Open(no_water_frcov_file)
    no_water_frcov  = no_water_gdal.ReadAsArray()
    no_water_frcov[no_water_frcov==rfl.no_data] = np.nan

    water = (rfl.ndi(550,850) > 0) & (rfl.get_wave(900) < .15)

    filter_frcov = np.zeros((rfl.lines,rfl.columns,4))
    filter_frcov[water,:3] =no_snow_frcov[:3,water].T
    filter_frcov[~water,0] =no_water_frcov[0,~water].T
    filter_frcov[~water,1] =no_water_frcov[1,~water].T
    filter_frcov[~water,3] =no_water_frcov[2,~water].T
    filter_frcov[~rfl.mask['no_data']] = -9999

    filter_frcov_file = f'{frcov_img_path}_fractional_cover'
    band_names = ['soil','vegetation','water','snow_ice']
    header = rfl.get_header()
    header['bands'] = 4
    header['band names'] = band_names
    header['wavelength'] = []
    header['fwhm'] = []
    header['bbl'] = []

    writer = ht.io.envi.WriteENVI(filter_frcov_file, header)
    for band_num in range(4):
        writer.write_band(filter_frcov[:,:,band_num], band_num)

    cover_counts = {}

    # Load bands and calculate cover percentile counts
    for band_num in range(0,4):
        frac_dict = {}
        if ("DESIS" in corfl_img_path) and (band_num ==3):
            band_arr = np.zeros((rfl.lines,rfl.columns))
        else:
            band_arr = filter_frcov[:,:,band_num]
        for percent in np.linspace(0,1,11):
            frac_dict[round(percent,1)] = float((band_arr >= percent).sum())
        cover_counts[band_names[band_num]] = frac_dict

    run_config["metadata"]["cover_percentile_counts"] = cover_counts

    # Generate metadata in .met.json file
    frcov_met_json_path = f"output/{frcov_basename}.met.json"
    print(f"Generating metadata from runconfig to {frcov_met_json_path}")
    generate_metadata(run_config, frcov_met_json_path)

    # Generate quicklook
    frcov_ql_path = f"output/{frcov_basename}.png"
    print(f"Generating quicklook to {frcov_ql_path}")

    rgb=  np.array(filter_frcov[:,:,:3])
    rgb[rgb == rfl.no_data] = np.nan
    rgb = (rgb*255).astype(np.uint8)
    im = Image.fromarray(rgb)
    im.save(frcov_ql_path)

    #Convert to COG
    temp_file =  'work/temp_frcover.tif'
    out_file =  f"output/{frcov_basename}.tif"

    print(f"Creating COG {out_file}")

    units = ['PERCENT',
            'PERCENT',
            'PERCENT',
            'PERCENT',]

    descriptions=  ['SOIL PERCENT COVER',
                      'VEGETATION PERCENT COVER',
                      'WATER PERCENTCOVER',
                      'SNOW/ICE PERCENT COVER']

    # Set the output raster transform and projection properties
    driver = gdal.GetDriverByName("GTIFF")
    tiff = driver.Create(temp_file,
                          no_snow_gdal.RasterXSize,
                          no_snow_gdal.RasterYSize,
                          4,
                          gdal.GDT_Float32)

    tiff.SetGeoTransform(no_snow_gdal.GetGeoTransform())
    tiff.SetProjection(no_snow_gdal.GetProjection())
    tiff.SetMetadataItem("DESCRIPTION","FRACTIONAL COVER")

    # Write bands to file
    for i,band_name in enumerate(band_names):
        out_band = tiff.GetRasterBand(i+1)
        out_band.WriteArray(filter_frcov[:,:,i])
        out_band.SetDescription(band_name)
        out_band.SetNoDataValue(rfl.no_data)
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
