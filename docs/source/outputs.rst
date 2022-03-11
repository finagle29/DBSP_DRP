**********************
Data Reduction Outputs
**********************

Assuming you ran DBSP_DRP with ``dbsp_reduce -r $RAW_PATH -d $OUTPUT_PATH``,
then in the ``$OUTPUT_PATH/Science`` directory you will find:

PypeIt 2D Spectrum files ``spec2d_redNNNN-target_DBSPr_obstimestamp.fits``
described in more detail `here <https://pypeit.readthedocs.io/en/latest/out_spec2D.html>`__.
Briefly, these files hold flat-fielded 2D spectral images, as well as sky, noise, object,
wavelength, and bad pixel mask images.
These files can be visually inspected using the command ``pypeit_show_2dspec SPEC2D_FILE``.

PypeIt 1D Spectrum files ``spec1d_redNNNN-target_DBSPr_obstimstamp.fits``
described in more detail `here <https://pypeit.readthedocs.io/en/latest/out_spec1D.html>`__.
Briefly, these files hold 1D spectra for each object that was traced and extracted from the
raw frame. Each object's spectrum is stored in a separate extension of the FITS file, and the
extension names are of the form ``SPATNNNN-SLITMMMM-DET01`` where the SPAT number describes the
spatial pixel coordinate of the object and the SLIT number is only useful for spectrographs
with multiple slits (and the DET number is only useful for spectrographs with multiple detectors
for the same arm).
These files can be visually inspected using the command ``pypeit_show_1dspec --exten N SPEC1D_FILE``
to view extension N of the file.

PypeIt 1D Coadd files ``redNNNN-redNNNN_target_SPATNNNN.fits`` described in more detail
`here <https://pypeit.readthedocs.io/en/latest/coadd1d.html#current-coadd1d-data-model>`__.
These files exist to separate out multiple objects on the same frame into their own file, and to
coadd consecutive exposures of the same frame.
These files can be visually inspected using the command ``lt_xspec COADD_FILE``.
The filename contains the first and last raw data file and the medial spatial pixel coordinate of the object.
This was done to make the coadd filename of constant length, instead of scaling with the number of coadded frames,
since all operating systems/file systems have maximum allowed filename lengths.

Telluric-corrected 1D Coadd files have ``_tellcorr`` appended to the base filename of the coadd file.

In the directory ``$OUTPUT_PATH/spliced`` are the spliced final data products.
The final data product of DBSP_DRP is a FITS file named ``target_char.fits`` with structure described below,
where ``char`` is a one-letter designation of which object along the slit it is.

.. table:: table
    :widths: 16 7 20 20

    ================ ======== ======================================================= =================================================
    Extension Number Name     Header                                                  Data
    ================ ======== ======================================================= =================================================
    0                PRIMARY  Version info about DBSP_DRP, PypeIt, NumPy and Astropy  None
    1                RED0031  Header from raw red side FITS file                      Fluxed (not telluric-corrected) red side spectrum
    2                RED0032  Header from raw red side FITS file                      Fluxed (not telluric-corrected) red side spectrum
    3                BLUE0032 Header from raw blue side FITS file                     Fluxed blue side spectrum
    4                BLUE0032 Header from raw blue side FITS file                     Fluxed blue side spectrum
    5                RED      Header from raw red side FITS file                      Fluxed (and telluric-corrected) red side spectrum
    6                BLUE     Header from raw blue side FITS file                     Fluxed blue side spectrum
    7                SPLICED  Empty                                                   Final spliced spectrum
    ================ ======== ======================================================= =================================================

The header on the 0th extension also contains cards named ``B_COADD`` and ``R_COADD`` which contain the filename
of the blue and red coadd files, respectively, that were spliced together. This is useful for determining which
traces a particular final output file corresponds to.
The 0th header also contains an ``INTERP_GAPS`` cards, noting whether or not detector gaps were interpolated over
during splicing.
If the splicing has been manually adjusted (see :doc:`adjust_splicing` for more details) then a ``RED_MULT`` card
will also be present, recording the factor multiplied into the red coadd spectrum before splicing with the blue coadd.

If *n* red side files were coadded and *m* blue side files were coadded, then extensions 1 through *n* would contain the
red side raw headers and fluxed spectrum from each individual file, extensions 1+*n* through 1+*n*+*m* would contain the
blue side raw headers and fluxed spectra from each individual file, and the last 3 extensions are the red coadd, blue
coadd, and final spliced spectrum.

If an object was not observed in the blue (red) then there will be no raw blue (red) frames, and the BLUE (RED) extension
would still exist, but contain no data.

Each of the data tables contain columns for `wave`, `flux`, and `sigma`, with wavelength in Angstroms
and both flux and sigma in units of 10\ :sup:`-17`\ erg/s/cm\ :sup:`2`\ /Ang.
