**************
DBSP Quicklook
**************

During an observing run, in order to make time-sensitive decisions about what to observe,
a quicklook script is provided.

Usage
*****

.. code-block :: console

    $ dbsp_ql --help
    usage: dbsp_ql [-h] [--no-show] fname

    Quicklook for P200 DBSP

    positional arguments:
      fname       file to take a quick look at, or else red/blue
                  to just perform rough calibrations

    optional arguments:
      -h, --help  show this help message and exit
      --no-show   Set this flag to suppress opening of plots

First, navigate to a directory you want the output data in.

.. code-block :: console

    $ cd /path/to/workdir

Then, with at least one arc frame and at least one flat frame in ``/path/to/calibs``
run

.. code-block :: console

    $ dbsp_ql /path/to/calibs/red
    $ dbsp_ql /path/to/calibs/blue

to perform quick calibrations for the red and blue sides, respectively.

Then once you have a science frame (either in the same directory, or elsewhere)
in ``/path/to/science/data``, run

.. code-block :: console

    $ dbsp_ql /path/to/science/data/red0030.fits
    $ dbsp_ql /path/to/science/data/blue0030.fits

This step should be very quick (about 15 seconds for the red side, 8 seconds for the blue side),
and each command will pop up a ginga window for you to inspect the sky-subtracted frame.
Also, a GUI will pop up to display the fluxed spectrum for each object identified in the frame.

The quicklook script produces PypeIt 2D Spectrum and 1D Spectrum files, described in :doc:`outputs`.
