***************************************
DBSP_DRP Quality Assurance (QA) Outputs
***************************************

The QA folder has two main QA pages: ``MF_A.html`` and ``Extraction.html``.

``MF_A.html`` is automatically generated by PypeIt and includes a number of
QA plots regarding the calibration data, which currently only encompass the
wavelength calibration. You can read in more detail about PypeIt's QA plots
`here <https://pypeit.readthedocs.io/en/latest/qa.html>`__.

In ``Extraction.html``, for each target, there is a page that shows various
steps of the reduction for the red and blue sides. From left to right, the
images are flat-fielded frame, sky model, sky-subtracted frame, sky-subtracted
residuals, and sky- and object-subtracted residuals. Above each set of images,
the raw filename, time of observation and the airmass is displayed. Slit edges
are marked in red and green, object traces are marked in orange, and extraction
FWHMs are highlighted in orange. If the extraction FWHM is very different from
what you would expect, it is probably a good idea to use manual tracing on that
target.
