# SISTER Fractional Cover

## Description

The sister-fractional-cover repository is a wrapper for the L2B spectral unmixing repository called 
[SpectralUnmixing](https://github.com/emit-sds/SpectralUnmixing), which is a general, fast, and flexible spectral 
unmixing package.  It is oriented towards VSWIR imaging spectroscopy data but applicable for different sensor types. 
Also, it includes options for different treatments of endmember library assemblages, including MESMA and bootstrapping 
(aka monte carlo) strategies.

## Dependencies

This repository is built to run on SISTER (Space-based Imaging Spectroscopy and Thermal pathfindER), a data 
processing back-end that allows for the registration of algorithms as executable containers and execution of those 
containers at scale.  The manifest file that configures this repository for registration and describes all of its 
necessary dependencies is called `algorithm_config.yaml`.  In this file you will find:

* The repository URL and version to register
* The base Docker image which this repository gets installed into
* The build script which is used to install this repository into the base Docker image

Specific dependencies for executing the code in this repository can be found in both the Dockerfile and the build 
script.

In addition to the above dependencies, you will need access to the MAAP API via the maap-py library in order to 
register algorithms and submit jobs.  maap-py can be obtained by running:

    git clone --single-branch --branch system-test-8 https://gitlab.com/geospec/maap-py.git

## PGE Arguments

The sister-fractional-cover PGE takes the following arguments:


| Argument      | Type   | Description                                                                        | Default |
|---------------|--------|------------------------------------------------------------------------------------|---------|
| l2a_rfl       | file   | S3 URL to the corrected reflectance dataset folder                                 | -       |
| n_cores       | config | Number of cores for parallelization                                                | 1       |
| refl_nodata   | config | Nodata value expected in input reflectance data                                    | None    |
| refl_scale    | config | Scale image data (divide it by) this amount                                        | None    |
| normalization | config | Flag to indicate the scaling type. Options = none, brightness, specific wavelength | None    |
| crid          | config | Composite Release ID to tag file names                                             | 000     |
| _force_ingest | config | Flag that allows overwriting existing files                                        | True    |

## Outputs

The L2B fractional cover PGE outputs a GeoTIFF containing layers for vegetation, soil, water, and snow and supporting 
files. The outputs of the PGE use the following naming convention:

    SISTER_INSTRUMENT_LEVEL_PRODUCT_YYYYMMDDTHHMMSS_CRID.EXTENSION

| Product                                                       | Format, Units | Example filename                                         |
|---------------------------------------------------------------|---------------|----------------------------------------------------------|
| Fractional cover raster image (vegetation, soil, water, snow) | GeoTIFF       | SISTER_AVNG_L2B_FRCOV_20220814T183137_000.tif            |
| Fractional cover metadata file                                | JSON          | SISTER_AVNG_L2B_FRCOV_20220814T183137_000.met.json       |
| Fractional cover browse image                                 | PNG           | SISTER_AVNG_L2B_FRCOV_20220814T183137_000.png            |
| PGE log file                                                  | Text          | SISTER_AVNG_L2B_FRCOV_20220814T183137_000.log            |
| PGE run config                                                | JSON          | SISTER_AVNG_L2B_FRCOV_20220814T183137_000.runconfig.json |

## Registering the Repository with SISTER

    from maap.maap import MAAP
    
    maap = MAAP(maap_host="34.216.77.111")
    
    algo_config_path = "sister-fractional-cover/algorithm_config.yaml"
    response = maap.register_algorithm_from_yaml_file(file_path=algo_config_path)
    print(response.text)

## Submitting a Job on SISTER

    from maap.maap import MAAP
    
    maap = MAAP(maap_host="34.216.77.111")
    
    job_response = maap.submitJob(
        algo_id="sister-fractional-cover",
        version="1.0.0",
        l2a_rfl="s3://s3.us-west-2.amazonaws.com:80/sister-ops-workspace/LOM/PRODUCTS/AVNG/L2A_CORFL/2022/08/14/SISTER_AVNG_L2A_CORFL_20220814T183137_000",
        n_cores=32,
        crid="000",
        normalization="brightness",
        publish_to_cmr=False,
        cmr_metadata={},
        queue="sister-job_worker-8gb",
        identifier="SISTER_AVNG_L2B_FRCOV_20220814T183137_000")
    
    print(isofit_job_response.id, isofit_job_response.status)
