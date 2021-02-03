[![Documentation Status](https://readthedocs.org/projects/dbsp-drp/badge/?version=latest)](https://dbsp-drp.readthedocs.io/en/latest/?badge=latest)
![Test](https://github.com/finagle29/DBSP_DRP/workflows/Test/badge.svg)

# DBSP_DRP


## Description
DBSP_DRP is a Data Reduction Pipeline for Palomar's workhorse spectrograph DBSP.
It is built on top of [PypeIt](https://github.com/pypeit/PypeIt).
DBSP_DRP automates the reduction, fluxing, telluric correction, and combining of the red and blue sides of one night's
data.
It adds two GUIs to allow for easier control of your reduction:
- select which data to reduce, and verify the correctness of your FITS headers in an editable table GUI
- manually place traces for a sort of manually "forced" spectroscopy with the `-m` option

## Prerequisites
DBSP_DRP's dependencies are detailed in [environment.yml](environment.yml).
You can install all prerequisites by downloading the environment.yml file, navigating to the directory containing it in your terminal window and running
```shell_session
$ conda env create -f environment.yml
```

The telluric correction code provided by PypeIt relies on a large (5 GB) atmospheric model file,
TellFits_Lick_3100_11100_R10000.fits, which can be downloaded
[here](https://drive.google.com/drive/folders/1x5d2_L8pwLDmvvoFUCa-vIoluv3GpowA)
and must be installed into the ``pypeit/data/telluric/`` directory of your PypeIt installation.

An easier alternative is to use the [download_tellfile](bin/download_tellfile) script to download and install the atmospheric model file for you.

## Installation
You can install using `pip`
```shell_session
$ pip install git+https://github.com/finagle29/DBSP_DRP.git
```

Or you can install from source
```shell_session
$ git clone https://github.com/finagle29/DBSP_DRP.git
$ cd DBSP_DRP
$ pip install -e .
```

## Usage
```shell_session
$ dbsp_reduce -r /path/to/data/DBSP_YYYYMMDD -d /path/to/data/DBSP_YYYYMMDD_redux
    [-a --arm RED or BLUE -j --jobs N -p --parameter-file PARAMS.txt]
    [-i --no-interactive -m --manual-extraction -t --skip-telluric]
```
