[![Documentation Status](https://readthedocs.org/projects/dbsp-drp/badge/?version=latest)](https://dbsp-drp.readthedocs.io/en/latest/?badge=latest)
[![Test](https://github.com/finagle29/DBSP_DRP/actions/workflows/test.yml/badge.svg)](https://github.com/finagle29/DBSP_DRP/actions/workflows/test.yml)

![PyPI version](https://img.shields.io/pypi/v/DBSP_DRP?label=PyPI&logo=python&logoColor=white)
![conda-forge version](https://img.shields.io/conda/vn/conda-forge/dbsp_drp?logo=conda-forge&logoColor=white)
![pip downloads](https://img.shields.io/pypi/dm/DBSP_DRP)
![conda downloads](https://img.shields.io/conda/dn/conda-forge/DBSP_DRP?label=conda%20downloads)

# DBSP_DRP


## Description
DBSP_DRP is a Data Reduction Pipeline for Palomar's workhorse spectrograph DBSP.
It is built on top of [PypeIt](https://github.com/pypeit/PypeIt).
DBSP_DRP automates the reduction, fluxing, telluric correction, and combining of the red and blue sides of one night's
data.
It adds several GUIs to allow for easier control of your reduction:
- select which data to reduce, and verify the correctness of your FITS headers in an editable table GUI
- manually place traces for a sort of manually "forced" spectroscopy with the `-m` option
- after manually placing traces, manually select sky regions and tweak the FWHM of your manual traces

The latest documentation can be found on [Read the Docs](https://dbsp-drp.readthedocs.io/en/latest/index.html).

## Citation
If you use DBSP_DRP in your research, please cite the following publication,
and optionally this repository itself.

```
@misc{dbsp_drp2021_arxiv,
      title={DBSP_DRP: A Python package for automated spectroscopic data reduction of DBSP data}, 
      author={Milan S. Roberson and Christoffer Fremling and Mansi M. Kasliwal},
      year={2021},
      eprint={2107.12339},
      archivePrefix={arXiv},
      primaryClass={astro-ph.IM}
}
```

```
@misc{dbsp_drp2021_github,
    title     = "DBSP_DRP",
    year      = "2021",
    publisher = "GitHub",
    url       = "https://github.com/finagle29/dbsp_drp"}
  }
```

and please also cite [PypeIt](https://github.com/pypeit/PypeIt#citation).

## Prerequisites
DBSP_DRP's dependencies are detailed in [environment.yml](environment.yml).
You can install all prerequisites for a `pip` or source install by downloading the environment.yml file, navigating to the directory containing it in your terminal window and running
```shell_session
$ conda env create -f environment.yml
```
Installing DBSP_DRP using `conda` does not require this step.

The telluric correction code provided by PypeIt relies on a large (5 GB) atmospheric model file,
TellFits_Lick_3100_11100_R10000.fits, which can be downloaded
[here](https://drive.google.com/drive/folders/1FFRWjUZ58HiDuDD33MYqBzMWDQanBRRy)
and must be installed into the ``pypeit/data/telluric/`` directory of your PypeIt installation.

An easier alternative is to use the [download_tellfile](bin/download_tellfile) script to download and install the atmospheric model file for you.

## Installation
You can install using `conda`
```shell_session
$ conda install -c conda-forge dbsp_drp
```

or `pip`
```shell_session
$ pip install dbsp-drp
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
    [-a {red,blue}] [-i] [-m] [--debug] [-j N] [-p PARAMETER_FILE] [-t] [-c]
    [--splicing-interpolate-gaps]
```
