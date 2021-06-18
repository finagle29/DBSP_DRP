[![Documentation Status](https://readthedocs.org/projects/dbsp-drp/badge/?version=latest)](https://dbsp-drp.readthedocs.io/en/latest/?badge=latest)
[![Build Status](https://travis-ci.org/finagle29/DBSP_DRP.svg?branch=master)](https://travis-ci.org/finagle29/DBSP_DRP)

# DBSP_DRP


## Description
DBSP_DRP is a Data Reduction Pipeline for Palomar's workhorse spectrograph DBSP.
It is built on top of [PypeIt](https://github.com/pypeit/PypeIt) and adds an interactive header validation GUI as well as automating the reduction and fluxing of one night's data.

## Citation
If you use DBSP_DRP in your research, please also reference [PypeIt](https://github.com/pypeit/PypeIt#citation).

## Prerequisites
DBSP_DRP's dependencies are detailed in [environment.yml](environment.yml).
You can install all prerequisites by downloading the environment.yml file, navigating to the directory containing it in your terminal window and running
```shell_session
$ conda env create -f environment.yml
```

The telluric correction code provided by PypeIt relies on a large (5 GB) atmospheric model file, which can be downloaded [here](https://drive.google.com/drive/folders/1x5d2_L8pwLDmvvoFUCa-vIoluv3GpowA)
and must be installed into the ``pypeit/data/telluric/`` directory of your PypeIt installation.

An easier alternative is to use the [download_tellfile](bin/download_tellfile) script to download and install the atmospheric model file for you.

## Installation
You can install using `pip`
```shell_session
$ pip install https://github.com/finagle29/DBSP_DRP/archive/master.tar.gz
```

Or you can install from source
```shell_session
$ cd /path/to/DBSP_DRP
$ python setup.py install
```

## Usage
```shell_session
$ dbsp_reduce -r /path/to/data/DBSP_YYYYMMDD -o /path/to/data/DBSP_YYYYMMDD_redux
```
