**********************
Data Reduction Outputs
**********************

Assuming you ran DBSP_DRP with ``dbsp_reduce -r $RAW_PATH -d $OUTPUT_PATH``,
then in the ``$OUTPUT_PATH/Science`` directory you will find

PypeIt 2D Spectrum files ``spec2d_redNNNN-target_DBSPr_obstimestamp.fits``
described in more detail `here <https://pypeit.readthedocs.io/en/latest/out_spec2D.html>`_.
Briefly, these files hold flat-fielded 2d spectral images, as well as sky, noise, object,
wavelength, and bad pixel mask images.
These files can be visually inspected using the command ``pypeit_show_2dspec SPEC2D_FILE``.

PypeIt 1D Spectrum files ``spec1d_redNNNN-target_DBSPr_obstimstamp.fits``
described in more detail `here <https://pypeit.readthedocs.io/en/latest/out_spec1D.html>`_.
Briefly, these files hold 1d spectra for each object that was traced and extracted from the
raw frame. Each object's spectrum is stored in a separate extension of the FITS file, and the
extension names are of the form ``SPATNNNN-SLITMMMM-DET01`` where the SPAT number describes the
spatial pixel coordinate of the object and the SLIT number is only useful for spectrographs
with multiple slits (and the DET number is only useful for spectrographs with multiple detectors
for the same arm).
These files can be visually inspected using the command ``pypeit_show_1dspec --exten N SPEC1D_FILE``
to view extension N of the file.

PypeIt 1D Coadd files ``redNNNN-target_SPATNNNN-SLITMMMM-DET01.fits`` described in more detail
`here <https://pypeit.readthedocs.io/en/latest/coadd1d.html#current-coadd1d-data-model>`_.
These files only exist to separate out multiple objects on the same frame into their own file.
These files can be visually inspected using the command ``lt_xspec COADD_FILE``.

Telluric-corrected 1D Coadd files have ``_tellcorr`` appended to the base filename of the coadd file.

The final data product of DBSP_DRP is a FITS file named ``target.fits`` with structure described below.

.. table:: table
    :widths: 16 7 20 20

    ================ ======== ======================================================= =================================================
    Extension Number Name     Header                                                  Data
    ================ ======== ======================================================= =================================================
    0                PRIMARY  Version info about DBSP_DRP, PypeIt, numpy and astropy  None
    1                RED      Header from raw red side fits file                      Fluxed (and telluric-corrected) red side spectrum
    2                BLUE     Header from raw blue side fits file                     Fluxed blue side spectrum
    3                SPLICED  Empty                                                   Final spliced spectrum
    ================ ======== ======================================================= =================================================

Each of the data tables contain columns for `wave`, `flux`, and `sigma`, with wavelength in Angstroms
and both flux and sigma in units of 10\ :sup:`-17`\ erg/s/cm\ :sup:`2`\ /Ang.
