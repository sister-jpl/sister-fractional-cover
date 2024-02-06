# SISTER Fractional Cover

## Description

The sister-fractional-cover repository is a wrapper for the L2B spectral unmixing repository called 
[SpectralUnmixing](https://github.com/emit-sds/SpectralUnmixing), which is a general, fast, and flexible spectral 
unmixing package.  It is oriented towards VSWIR imaging spectroscopy data but applicable for different sensor types. 
Also, it includes options for different treatments of endmember library assemblages, including MESMA and bootstrapping 
(aka monte carlo) strategies.

## PGE Arguments

The sister-fractional-cover PGE takes the following arguments:


| Argument            | Description                                                                        | Default |
|---------------------|------------------------------------------------------------------------------------|---------|
| reflectance_dataset | S3 URL to the corrected reflectance dataset folder                                 | -       |
| n_cores             | Number of cores for parallelization                                                | 1       |
| refl_scale          | Scale image data (divide it by) this amount                                        | none    |
| normalization       | Flag to indicate the scaling type. Options = none, brightness, specific wavelength | none    |
| crid                | Composite Release ID to tag file names                                             | 000     |
| experimental        | Designates outputs as "experimental"                                               | 'True'  |

## Outputs

The L2B fractional cover PGE outputs a GeoTIFF containing layers for vegetation, soil, water, and snow and supporting 
files. The outputs of the PGE use the following naming convention:

    (EXPERIMENTAL-)SISTER_INSTRUMENT_LEVEL_PRODUCT_YYYYMMDDTHHMMSS_CRID.EXTENSION

Note that the "EXPERIMENTAL-" prefix is optional and is only added when the "experimental" flag is set to True.

The following data products are produced:

| Product                                                       | Format  | Example filename                                         |
|---------------------------------------------------------------|---------|----------------------------------------------------------|
| Fractional cover raster image (vegetation, soil, water, snow) | GeoTIFF | SISTER_AVCL_L2B_FRCOV_20110513T175417_000.tif            |
| Fractional cover metadata file (STAC formatted)               | JSON    | SISTER_AVCL_L2B_FRCOV_20110513T175417_000.json           |
| Fractional cover browse image                                 | PNG     | SISTER_AVCL_L2B_FRCOV_20110513T175417_000.png            |
| PGE log file                                                  | Text    | SISTER_AVCL_L2B_FRCOV_20110513T175417_000.log            |
| PGE run config                                                | JSON    | SISTER_AVCL_L2B_FRCOV_20110513T175417_000.runconfig.json |

Metadata files are [STAC formatted](https://stacspec.org/en) and compatible with tools in the [STAC ecosystem](https://stacindex.org/ecosystem).

## Executing the Algorithm

This algorithm requires [Anaconda Python](https://www.anaconda.com/download)

To install and run the code, first clone the repository and execute the install script:

    git clone https://github.com/sister-jpl/sister-fractional-cover.git
    cd sister-fractional-cover
    ./install.sh
    cd ..

Then, create a working directory and enter it:

    mkdir WORK_DIR
    cd WORK_DIR

Copy input files to the work directory. For each "dataset" input, create a folder with the dataset name, then download 
the data file(s) and STAC JSON file into the folder.  For example, the reflectance dataset input would look like this:

    WORK_DIR/SISTER_AVCL_L2A_CORFL_20110513T175417_000/SISTER_AVCL_L2A_CORFL_20110513T175417_000.bin
    WORK_DIR/SISTER_AVCL_L2A_CORFL_20110513T175417_000/SISTER_AVCL_L2A_CORFL_20110513T175417_000.hdr
    WORK_DIR/SISTER_AVCL_L2A_CORFL_20110513T175417_000/SISTER_AVCL_L2A_CORFL_20110513T175417_000.json

Finally, run the code 

    ../sister-fractional-cover/run.sh --reflectance_dataset SISTER_AVCL_L2A_CORFL_20110513T175417_000
