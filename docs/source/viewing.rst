***********************
Viewing Reduced Spectra
***********************

``DBSP_DRP`` comes bundled with ``dbsp_show``, an all-purpose interactive spectrum plotter for
``DBSP_DRP`` and intermediate ``PypeIt`` files, built using Matplotlib.

.. code-block :: console

    $ dbsp_show --help
    usage: dbsp_show [-h] [--extension EXTENSION] [--BOX] [--COUNTS] fname

    All-purpose interactive spectrum plotter for DBSP_DRP and intermediate PypeIt files.

    positional arguments:
      fname                 path to FITS file

    optional arguments:
      -h, --help            show this help message and exit
      --extension EXTENSION
                            Extension name or number
      --BOX                 Use boxcar extraction when plotting PypeIt spec1d files. By default the optimal extraction is plotted.
      --COUNTS              Plot counts when plotting PypeIt spec1d files. By default flux is plotted.
