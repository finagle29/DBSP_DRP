from setuptools import setup

setup(
    name="DBSP_DRP",
    version="0.0.1dev",
    description="Data Reduction Pipeline for Palomar's Double Spectrograph",
    author="Milan S Roberson",
    author_email="mroberso@caltech.edu",
    license="GPL3",
    packages=['dbsp_drp'],
    scripts=['bin/dbsp_reduce'],
    package_data={'dbsp_drp': ['data/*']},
    install_requires=[
        "numpy",
        "astropy",
        "scipy",
        "matplotlib",
        "numba",
        "pyyaml",
        "configobj",
        "scikit-learn",
        "ipython",
        "ginga",
        "packaging",
        "linetools",
        "extension-helpers",
        "pytest",
        "yattag",
        "pypeit"
    ]
)
