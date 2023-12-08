#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SISTER
Space-based Imaging Spectroscopy and Thermal PathfindER
Author: Winston Olson-Duvall
"""

import datetime as dt
import glob
import json
import os
import subprocess
import shutil
import sys
import hytools as ht
import pandas as pd
import numpy as np
import pystac

try:
    from osgeo import gdal
except:
    import gdal
from PIL import Image

def find_nearest_day(day, subset_days):
    distances = []
    for special_day in subset_days:
        diff = (special_day - day) % 365
        distances.append(min(diff, 365 - diff))
    nearest_day = subset_days[distances.index(min(distances))]
    return np.argwhere(nearest_day==subset_days)[0][0]

def get_frcov_basename(corfl_basename, crid):
    # Replace product type
    tmp_basename = corfl_basename.replace("L2A_CORFL", "L2B_FRCOV")
    # Split, remove old CRID, and add new one
    tokens = tmp_basename.split("_")[:-1] + [str(crid)]
    return "_".join(tokens)

def generate_metadata(run_config, frcov_met_json_path, disclaimer):
    # Create .met.json file from runconfig for fractional cover
    metadata = run_config["metadata"]
    metadata["product"] = "FRCOV"
    metadata["processing_level"] = "L2B"
    metadata["description"] = f"{disclaimer}Fractional cover (soil, vegetation, water, snow)"
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
    experimental = run_config['inputs']['experimental']
    if experimental:
        disclaimer = "(DISCLAIMER: THIS DATA IS EXPERIMENTAL AND NOT INTENDED FOR SCIENTIFIC USE) "
    else:
        disclaimer = ""

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

    corfl_basename = os.path.basename(run_config["inputs"]["reflectance_dataset"])
    frcov_basename = get_frcov_basename(corfl_basename, run_config["inputs"]["crid"])
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

    snow_clim_file = f'{os.path.dirname(sister_frcov_dir)}/LIN10A1_snow_climatology_13day.tif'
    snow_clim = gdal.Open(snow_clim_file)
    snow = snow_clim.GetRasterBand(1)
    snow_days = np.arange(1,365,13)

    ulx,pixel_size,_,uly,_,_ = snow_clim.GetGeoTransform()

    bbox_min_x,bbox_min_y = np.min(run_config['metadata']['bounding_box'],axis=0)
    bbox_max_x,bbox_max_y = np.max(run_config['metadata']['bounding_box'],axis=0)

    x_offset = int((bbox_min_x-ulx)/pixel_size)
    width = int((bbox_max_x-ulx)/pixel_size) -x_offset

    y_offset = int((uly-bbox_max_y)/pixel_size)
    height = int((uly-bbox_min_y)/pixel_size) - y_offset

    subset = snow.ReadAsArray(x_offset, y_offset,
                              width, height)

    datetime = dt.datetime.strptime(run_config['metadata']['start_time'], '%Y-%m-%dT%H:%M:%SZ')
    doy = datetime.timetuple().tm_yday
    bit = find_nearest_day(doy, snow_days)
    snow_mask =(subset >> bit) &1
    snow_present = snow_mask.sum() >0

    if snow_present:
        print('Snow expected')
        endmember_files = ['no_water','no_snow']
    else:
        endmember_files = ['no_snow']

    for endmember_file in endmember_files:
        # Add required args
        cmd = ["julia"]

        # Add parallelization if n_cores > 1
        if run_config["inputs"]["n_cores"] > 1:
            cmd += ["-p", str(run_config["inputs"]["n_cores"])]

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
            f"--refl_scale={run_config['inputs']['refl_scale']}"]

        print("Running unmix.jl command: " + " ".join(cmd))
        subprocess.run(" ".join(cmd), shell=True)


    no_snow_frcov_file = f'{frcov_img_path}_no_snow_fractional_cover'
    no_snow_gdal = gdal.Open(no_snow_frcov_file)
    no_snow_frcov  = no_snow_gdal.ReadAsArray()
    no_snow_frcov[no_snow_frcov==rfl.no_data] = np.nan

    filter_frcov = np.zeros((rfl.lines,rfl.columns,4))

    if snow_present:
        no_water_frcov_file = f'{frcov_img_path}_no_water_fractional_cover'
        no_water_gdal = gdal.Open(no_water_frcov_file)
        no_water_frcov  = no_water_gdal.ReadAsArray()
        no_water_frcov[no_water_frcov==rfl.no_data] = np.nan

        water = (rfl.ndi(550,850) > 0) & (rfl.get_wave(900) < .15)

        filter_frcov[water,:3] =no_snow_frcov[:3,water].T
        filter_frcov[~water,0] =no_water_frcov[0,~water].T
        filter_frcov[~water,1] =no_water_frcov[1,~water].T
        filter_frcov[~water,3] =no_water_frcov[2,~water].T

    else:
      filter_frcov[:,:,:3]  = np.moveaxis(no_snow_frcov[:3],0,-1)

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
    fc_description = f"{disclaimer}FRACTIONAL COVER"
    tiff.SetMetadataItem("DESCRIPTION", fc_description)

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

    # If experimental, prefix filenames with "EXPERIMENTAL-"
    if experimental:
        for file in glob.glob(f"output/SISTER*"):
            shutil.move(file, f"output/EXPERIMENTAL-{os.path.basename(file)}")

    # Update the path variables if now experimental
    frcov_file = glob.glob("output/*%s.tif" % run_config['inputs']['crid'])[0]
    out_runconfig = glob.glob("output/*%s.runconfig.json" % run_config['inputs']['crid'])[0]
    log_path = glob.glob("output/*%s.log" % run_config['inputs']['crid'])[0]
    frcov_basename = os.path.basename(frcov_file)[:-4]

    # Generate STAC
    catalog = pystac.Catalog(id=corfl_basename,
                             description=f'{disclaimer}This catalog contains the output data products of the SISTER '
                                         f'fractional cover PGE, including a fractional cover cloud-optimized GeoTIFF. '
                                         f'Execution artifacts including the runconfig file and execution '
                                         f'log file are included with the fractional cover data.')

    # Add items for data products
    tif_files = glob.glob("output/*SISTER*.tif")
    tif_files.sort()
    for tif_file in tif_files:
        metadata = generate_stac_metadata(frcov_basename, fc_description, run_config["metadata"])
        assets = {
            "cog": f"./{os.path.basename(tif_file)}",
        }
        # If it's the fractional cover product, then add png, runconfig, and log
        if os.path.basename(tif_file) == f"{frcov_basename}.tif":
            png_file = tif_file.replace(".tif", ".png")
            assets["browse"] = f"./{os.path.basename(png_file)}"
            assets["runconfig"] = f"./{os.path.basename(out_runconfig)}"
            if os.path.exists(log_path):
                assets["log"] = f"./{os.path.basename(log_path)}"
        item = create_item(metadata, assets)
        catalog.add_item(item)

    # set catalog hrefs
    catalog.normalize_hrefs(f"./output/{frcov_basename}")

    # save the catalog
    catalog.describe()
    catalog.save(catalog_type=pystac.CatalogType.SELF_CONTAINED)
    print("Catalog HREF: ", catalog.get_self_href())

    # Move the assets from the output directory to the stac item directories
    for item in catalog.get_items():
        for asset in item.assets.values():
            fname = os.path.basename(asset.href)
            shutil.move(f"output/{fname}", f"output/{frcov_basename}/{item.id}/{fname}")


def generate_stac_metadata(basename, description, in_meta):

    out_meta = {}
    out_meta['id'] = basename
    out_meta['start_datetime'] = dt.datetime.strptime(in_meta['start_time'], "%Y-%m-%dT%H:%M:%SZ")
    out_meta['end_datetime'] = dt.datetime.strptime(in_meta['end_time'], "%Y-%m-%dT%H:%M:%SZ")
    # Split corner coordinates string into list
    geometry = in_meta['bounding box']
    # Add first coord to the end of the list to close the polygon
    geometry.append(geometry[0])
    out_meta['geometry'] = geometry
    out_meta['properties'] = {
        'sensor': in_meta['sensor'],
        'description': description,
        'product': basename.split('_')[4],
        'processing_level': basename.split('_')[2]
    }
    return out_meta


def create_item(metadata, assets):
    item = pystac.Item(
        id=metadata['id'],
        datetime=metadata['start_datetime'],
        start_datetime=metadata['start_datetime'],
        end_datetime=metadata['end_datetime'],
        geometry=metadata['geometry'],
        bbox=None,
        properties=metadata['properties']
    )
    # Add assets
    for key, href in assets.items():
        item.add_asset(key=key, asset=pystac.Asset(href=href))
    return item


if __name__ == "__main__":
    main()
