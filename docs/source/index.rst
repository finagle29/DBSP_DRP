.. DBSP_DRP documentation master file, created by
   sphinx-quickstart on Fri Aug 21 13:14:00 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. |astropy| image:: http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat
    :target: http://www.astropy.org/

.. |pypi| image:: https://img.shields.io/pypi/v/DBSP_DRP?label=PyPI&logo=python&logoColor=white
    :target: https://pypi.org/project/DBSP-DRP/

.. |conda-forge| image:: https://img.shields.io/conda/vn/conda-forge/dbsp_drp?logo=conda-forge&logoColor=white
    :target: https://anaconda.org/conda-forge/dbsp_drp

DBSP_DRP
========

|conda-forge| |pypi| |astropy|

**Version**: |version|

DBSP_DRP is a Data Reduction Pipeline for Palomar's workhorse spectrograph DBSP.
It is built on top of `PypeIt <https://github.com/pypeit/PypeIt>`__.
DBSP_DRP automates the reduction, fluxing, telluric correction, and combining
of the red and blue sides of one night's data.
It adds several GUIs to allow for easier control of your reduction:

- select which data to reduce, and verify the correctness of your FITS headers in an editable table GUI
- manually place traces for a sort of manually "forced" spectroscopy with the ``-m`` option
- after manually placing traces, manually select sky regions and tweak the FWHM of your manual traces

DBSP_DRP also provides a quicklook script for making real-time decisions during
an observing run, and can open a GUI displaying a minimally reduced exposure in
under 15 seconds.

Statement of Need
-----------------

Palomar Observatory, located near San Diego, CA, is a multinational observatory
with a broad user base.
Users come from large and small institutions, and their observing experience
ranges from novice to expert.
One responsibility for serving such a diverse user base is to provide software
data reduction pipelines for the most frequently used instruments, such as the
Palomar Double Spectrograph (DBSP).
Although DBSP was commissioned in 1982, it remains the workhorse instrument of
the 200” Hale Telescope.
It is used on 42% of the nights in a year, comprising nearly all of the
valuable “dark” (moonless) time.
In previous years, standard astronomical practice left the data reduction up to
the user.
However, attitudes in instrument building have shifted since DBSP was built.
The pipeline is now considered an indispensable component of the astronomical
instrument.
In fact, the difference between a good pipeline and a great pipeline means the
difference between counting some of the photons vs. counting all of the photons.

Spectroscopy is a severe bottleneck in time-domain astronomy; currently less
than 10% of discoveries are spectroscopically classified.
Without a pipeline, data reduction is a difficult process and the standard
method without a pipeline is to use IRAF, a 35 year old program on which
development and maintenance was discontinued in 2013 and whose use is
discouraged by many in the field.
Needless to say, data reduction sans pipeline is extremely time-consuming.
There is a clear need for a modern and stable automated data reduction pipeline
for DBSP.

During observing runs, one would like to be able to quickly inspect data as it
is taken, in order to ensure that it is of sufficient quality to do the desired
science with.
For objects whose brightness may have changed between a
previous observation and the observing run, the observer may have uncertainties
regarding how long of an exposure is needed to produce quality data.
For very faint objects or objects in crowded fields, the observer may not even
be sure that the telescope is pointed at the right object!
A quicklook functionality, that can do a rudimentary reduction to correct for
instrumental signatures and subtract light from the sky, revealing the spectra
of the objects observed, can answer questions of exposure time and whether the
object observed is the right one.


.. toctree::
   :maxdepth: 1
   :caption: Getting Started:

   community
   installing
   usage
   quicklook
   outputs

----

.. toctree::
   :maxdepth: 1
   :caption: Post-Reduction

   viewing
   trim

----

.. toctree::
   :maxdepth: 1
   :caption: Fine-tuning Your Reduction:

   manual_extraction
   qa
   adjust_splicing

----

.. toctree::
   :maxdepth: 1
   :caption: Installing DBSP_DRP On-site

   docker

----

.. toctree::
   :maxdepth: 2
   :caption: For Developers

   modules

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
