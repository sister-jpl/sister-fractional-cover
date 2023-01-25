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
import sys


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


def generate_quicklook(rfl_img_path, output_path):
    # TODO: Add quicklook code
    pass


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
    print("Making work directory")
    if not os.path.exists("work"):
        subprocess.run("mkdir work", shell=True)

    # Make output dir
    print("Making output directory")
    if not os.path.exists("output"):
        subprocess.run("mkdir output", shell=True)

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
    subprocess.run(f"cp input/{corfl_basename}/{corfl_basename}.bin {corfl_img_path}", shell=True)
    subprocess.run(f"cp input/{corfl_basename}/{corfl_basename}.hdr {corfl_hdr_path}", shell=True)

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
    generate_quicklook(frcov_img_path, frcov_ql_path)

    # TODO: Convert to COG
    # The fractional cover ENVI file is located at f"{frcov_img_path}_fractional_cover"

    # Move/rename files to output folders
    # TODO: Move COG files to output folder if needed
    subprocess.run(f"mv runconfig.json output/{frcov_basename}.runconfig.json", shell=True)


if __name__ == "__main__":
    main()
