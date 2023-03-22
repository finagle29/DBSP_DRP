[![Documentation Status](https://readthedocs.org/projects/dbsp-drp/badge/?version=latest)](https://dbsp-drp.readthedocs.io/en/latest/?badge=latest)
[![Test](https://github.com/finagle29/DBSP_DRP/actions/workflows/test.yml/badge.svg)](https://github.com/finagle29/DBSP_DRP/actions/workflows/test.yml)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.03612/status.svg)](https://doi.org/10.21105/joss.03612)

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
If you use DBSP_DRP in your research, please cite the following publications, or use the BibTeX provided below.
[![DOI](https://joss.theoj.org/papers/10.21105/joss.03612/status.svg)](https://doi.org/10.21105/joss.03612)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6241526.svg)](https://doi.org/10.5281/zenodo.6241526)

Additionally, please cite [PypeIt](https://github.com/pypeit/PypeIt#citation), with the BibTeX entries provided below (the Zenodo BibTex is for PypeIt 1.6.0, used in this version of DBSP_DRP).

### DBSP_DRP BibTeX
```
@article{dbsp_drp:joss,
  doi = {10.21105/joss.03612},
  url = {https://doi.org/10.21105/joss.03612},
  year = {2022},
  publisher = {The Open Journal},
  volume = {7},
  number = {70},
  pages = {3612},
  author = {Milan Sharma Mandigo-Stoba and Christoffer Fremling and Mansi M. Kasliwal},
  title = {DBSP_DRP: A Python package for automated spectroscopic data reduction of DBSP data},
  journal = {Journal of Open Source Software}
}
@misc{dbsp_drp:arxiv,
      title={DBSP_DRP: A Python package for automated spectroscopic data reduction of DBSP data}, 
      author={Milan Sharma Mandigo-Stoba and Christoffer Fremling and Mansi M. Kasliwal},
      year={2021},
      eprint={2107.12339},
      archivePrefix={arXiv},
      primaryClass={astro-ph.IM}
}
@software{dbsp_drp:zenodo,
  author       = {Mandigo-Stoba, Milan Sharma and
                  Fremling, Christoffer and
                  Kasliwal, Mansi M.},
  title        = {{DBSP\_DRP: A Python package for automated 
                   spectroscopic data reduction of DBSP data}},
  month        = feb,
  year         = 2022,
  publisher    = {Zenodo},
  version      = {v1.0.0},
  doi          = {10.5281/zenodo.6241526},
  url          = {https://doi.org/10.5281/zenodo.6241526}
}
```

### PypeIt BibTeX
```
@article{pypeit:joss_pub,
    doi = {10.21105/joss.02308},
    url = {https://doi.org/10.21105/joss.02308},
    year = {2020},
    publisher = {The Open Journal},
    volume = {5},
    number = {56},
    pages = {2308},
    author = {J. Xavier Prochaska and Joseph F. Hennawi and Kyle B. Westfall and Ryan J. Cooke and Feige Wang and Tiffany Hsyu and Frederick B. Davies and Emanuele Paolo Farina and Debora Pelliccia},
    title = {PypeIt: The Python Spectroscopic Data Reduction Pipeline},
    journal = {Journal of Open Source Software}
}

@software{pypeit:zenodov_v1_6,
  author       = {J. Xavier Prochaska and
                  Joseph Hennawi and
                  Ryan Cooke and
                  Kyle Westfall and
                  Feige Wang and
                  Debora Pelliccia and
                  EmAstro and
                  Milan Roberson and
                  T. E. Pickering and
                  tiffanyhsyu and
                  badpandabear and
                  Asher Wasserman and
                  Timothy Ellsworth Bowers and
                  Nicolas Tejos and
                  Alexa Villaume and
                  Brad Holden and
                  marijana777 and
                  Sunil Simha and
                  JT Schindler and
                  David Young and
                  Andreas Flörs and
                  Matt Wilde and
                  S.Tang and
                  Erik Tollerud and
                  Jacob Isbell and
                  Kristen Thyng and
                  Dan Foreman-Mackey and
                  David Jones and
                  Edward Betts and
                  Zlatan Vasović},
  title        = {pypeit/PypeIt: Version 1.6.0},
  month        = oct,
  year         = 2021,
  publisher    = {Zenodo},
  version      = {1.6.0},
  doi          = {10.5281/zenodo.5548381},
  url          = {https://doi.org/10.5281/zenodo.5548381}
}
```

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
