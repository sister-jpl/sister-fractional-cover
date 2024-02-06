#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SISTER
Space-based Imaging Spectroscopy and Thermal PathfindER
Author: Adam Chlus, Winston Olson-Duvall
"""

import argparse
import json
import os
import sys


def main():
    """
        This function takes as input the path to an inputs.json file and exports a run config json
        containing the arguments needed to run the SISTER ISOFIT PGE.

    """

    parser = argparse.ArgumentParser(description="Parse inputs to create runconfig.json")
    parser.add_argument("--reflectance_dataset", help="Path to reflectance dataset")
    parser.add_argument("--n_cores", help="Number of cores to use for parallelization", default="1")
    parser.add_argument("--refl_scale", help="Reflectance scale value", default="1.0")
    parser.add_argument("--crid", help="CRID value", default="000")
    parser.add_argument("--experimental", help="If true then designates data as experiemntal", default="True")
    args = parser.parse_args()

    run_config = {
        "inputs": {
            "reflectance_dataset": args.reflectance_dataset,
            "n_cores": int(args.n_cores),
            "refl_scale": args.refl_scale,
            "crid": args.crid,
        }
    }

    run_config["inputs"]["experimental"] = True if args.experimental.lower() == "true" else False

    # Add metadata to runconfig
    corfl_basename = os.path.basename(run_config["inputs"]["reflectance_dataset"])

    stac_json_path = os.path.join(run_config["inputs"]["reflectance_dataset"], f"{corfl_basename}.json")
    with open(stac_json_path, "r") as f:
        stac_item = json.load(f)
    run_config["metadata"] = stac_item["properties"]
    run_config["metadata"]["geometry"] = stac_item["geometry"]

    # Write out runconfig.json
    config_file = "runconfig.json"
    with open(config_file, "w") as outfile:
        json.dump(run_config, outfile, indent=4)


if __name__ == "__main__":
    main()
