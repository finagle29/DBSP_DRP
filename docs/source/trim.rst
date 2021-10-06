**************************
Trim Spectra By Wavelength
**************************

If you want to trim the wavelength range of your spliced spectra, e.g. to ignore data below the
atmospheric cutoff, the command-line tool ``dbsp_trim`` is there to help you!

.. code-block :: console

    $ dbsp_trim --help
    usage: dbsp_trim [-h] [--low LOW] [--high HIGH] fname [fname ...]

    Trim spliced spectrum to wavelength range.

    positional arguments:
    fname        Final data output from DBSP_DRP to trim. Can be multiple files.

    optional arguments:
      -h, --help   show this help message and exit
      --low LOW    Min wavelength after trimming. Default is 3300 Å.
      --high HIGH  Max wavelength after trimming. Default is 10500 Å.

To use ``dbsp_trim`` to trim all of your spliced spectra, located in the folder
``spliced``, you can run the following command:

.. code-block :: console

    $ dbsp_trim --low 3200 --high 11000 spliced/*.fits
