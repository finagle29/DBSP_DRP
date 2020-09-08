[![Documentation Status](https://readthedocs.org/projects/dbsp-drp/badge/?version=latest)](https://dbsp-drp.readthedocs.io/en/latest/?badge=latest)
[![Build Status](https://travis-ci.org/finagle29/DBSP_DRP.svg?branch=master)](https://travis-ci.org/finagle29/DBSP_DRP)

# DBSP_DRP


## Description
DBSP_DRP is a Data Reduction Pipeline for Palomar's workhorse spectrograph DBSP.
It is built on top of [PypeIt](https://github.com/pypeit/PypeIt) and adds an interactive header validation GUI as well as automating the reduction and fluxing of one night's data.

## Prerequisites
DBSP_DRP's dependencies are detailed in [environment.yml](environment.yml).
You can install all prerequisites by downloading the environment.yml file, navigating to the directory containing it in your terminal window and running
```shell_session
$ conda env create -f environment.yml
```

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
