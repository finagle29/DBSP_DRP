---
title: 'DBSP_DRP: A Python package for automated spectroscopic data reduction of DBSP data'
tags:
  - Python
  - astronomy
  - data reduction
  - spectroscopy
authors:
  - name: Milan S. Roberson[^1]
    orcid: 0000-0003-1118-3132
    affiliation: "1,2"
  - name: Christoffer Fremling
    orcid: 0000-0002-4223-103X
    affiliation: 2
  - name: Mansi M. Kasliwal
    orcid: 0000-0002-5619-4938
    affiliation: 2
affiliations:
  - name: Schmidt Academy of Software Engineering, California Institute of Technology
    index: 1
  - name: Division of Physics, Mathematics and Astronomy, California Institute of Technology
    index: 2
date: 26 October 2021
bibliography: paper.bib
nocite: |
  @Lunnan2020
---

[^1]: Present address: Department of Physics and Astronomy, University of California, Los Angeles.

# Summary

<!--
A summary describing the high-level functionality and purpose of the software
for a diverse, non-specialist audience.
-->

In astronomy, the spectrum of light emitted from astrophysical sources is of
great use, allowing astronomers to classify objects and measure their
properties.
To measure the spectrum of a source, astronomers use spectrographs, in which
dispersive elements spatially separate the incoming light by wavelength, and
detectors, most commonly CCDs, image this dispersed light.
But to do science with the spectrum, the 2D image in pixel coordinates taken by
the CCD must be converted into a 1D spectrum of flux vs. wavelength.
This process of converting 2D CCD images into 1D spectra is called data
reduction.

To increase the signal-to-noise ratio, astronomers can take multiple exposures
of the same object and coadd their 1D spectra to reveal faint absorption lines
or increase the precision with which an important emission line can be measured.
Many spectrographs have multiple paths that light can go through, and multiple
detectors, each measuring a particular part of the spectrum, to increase the
wavelength range that can be captured in a single exposure, or to allow the
high resolution observation of distinct wavelength ranges.
If two detectors cover an overlapping region, caused by partial reflectance of
a dichroic (wavelength-dependent beam splitter), then the spectra from each
detector need to be spliced together, combining the light collected by each
detector.

DBSP_DRP is a python package that provides fully automated data reduction of
data taken by the Double Spectrograph (DBSP) at the 200-inch Hale Telescope at
Palomar Observatory [@Oke&Gunn1982].
The underlying data reduction functionality to extract 1D spectra, perform flux
calibration and correction for atmospheric absorption, and coadd spectra
together is provided by PypeIt [@Prochaska2020].
The new functionality that DBSP_DRP brings is in orchestrating the complex data
reduction process by making smart decisions so that no user input is required
after verifying the correctness of the metadata in the raw FITS files in a
table-like GUI.
Though the primary function of DBSP_DRP is to automatically reduce an entire
night of data without user input, it has the flexibility for astronomers to
fine-tune the data reduction with GUIs for manually identifying the faintest
objects, as well as exposing the full set of PypeIt parameters to be tweaked
for users with particular science needs.
DBSP_DRP also handles some of the occasional quirks specific to DBSP, such as
swapping FITS header cards, adding (an) extra null byte/s to FITS files making
them not conform to the FITS specification, and not writing the coordinates of
the observation to file.
Additionally, DBSP_DRP contains a quicklook script for making real-time
decisions during an observing run, and it can open a GUI displaying a minimally
reduced exposure in under 15 seconds.
Docker containers are available for ease of deploying DBSP_DRP in its quicklook
configuration (without some large atmospheric model files) or in its full
configuration.

# Statement of Need

<!--
A Statement of Need section that clearly illustrates the research purpose of
the software.
-->

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
However, attitudes in instrument-building have shifted since DBSP was built.
A pipeline is now considered an indispensable component of an astronomical
instrument.
In fact, the difference between a good pipeline and a great pipeline means the
difference between counting some of the photons vs. counting all of the photons.

Spectroscopy is a severe bottleneck in time-domain astronomy; currently less
than 10% of discoveries are spectroscopically classified.
Without a pipeline, data reduction is a difficult process and the standard
method without a pipeline is to use IRAF [@IRAF86;@IRAF93], a 35-year-old
program on which development and maintenance was discontinued in 2013 and whose
use is discouraged by many in the field (e.g. @Ogaz2018).
Needless to say, data reduction sans existing pipeline is extremely
time-consuming.
There is a clear need for a modern and stable automated data reduction pipeline
for DBSP.

During observing runs, one would like to be able to quickly inspect data as it
is taken, in order to ensure that it is of sufficient quality to do the desired
science with.
For objects whose brightness may have changed between a previous observation
and the current observing run, or for nights with highly variable cloud cover,
the observer may be unsure how long of an exposure is needed to produce quality
data.
For very faint objects, objects in crowded fields, or objects with uncertain
positions (e.g. due to high or uncertain motion across the sky), the observer
may not even be sure that the telescope is pointed at the right object!
A quicklook functionality, which can do a rudimentary reduction to correct for
instrumental signatures and subtract light from the sky, revealing the spectra
of the objects observed, can answer questions of exposure time and whether the
object observed is the right one.

DBSP_DRP is currently being used by the ZTF Bright Transient Survey
[@Fremling2020;@Perley2020], the ZTF Census of the Local Universe [@De2020],
and a program investigating ZTF Superluminous Supernovae
(Lunnan et al., 2020; Chen et al., in preparation).
@Ravi2021arXiv is the first (known) publication that used DBSP_DRP for data
reduction.
The development of DBSP_DRP also lays the groundwork towards a fully automated
pipeline for the Next Generation Palomar Spectrograph that is planned to be
deployed on the Palomar 200-inch Hale Telescope in 2023.

# Acknowledgements

MSR acknowledges funding from the Schmidt Academy of Software Engineering,
which is supported by the generosity of Eric and Wendy Schmidt by
recommendation of the Schmidt Futures program.

We thank the following members of the time domain astronomy group at Caltech
for beta-testing and providing valuable feedback during the development of this
pipeline: Andy Tzanidakis, Lin Yan, Aishwarya Dahiwale, Yuhan Yao, Yashvi
Sharma, and Igor Andreoni.

MSR is extremely grateful to the welcoming, friendly, and helpful team of
developers on the PypeIt team, without whom this package would not exist.

This research made use of Astropy,[^2] a community-developed core Python
package for Astronomy [@astropy:2013;@astropy:2018].

[^2]: http://www.astropy.org

# References
