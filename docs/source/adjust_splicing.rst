*******************************
Adjusting splicing between arms
*******************************

Typically the splicing between arms is quite good, and you only need to adjust
it when the optimal extraction FWHM vary greatly between the red and blue sides.
This can happen when an extended object is near another object on the slit.
You can check the extraction FWHM visually by looking at the Extraction QA page.


.. code-block :: console

    $ dbsp_adjust_splicing --help
    usage: dbsp_adjust_splicing [-h] fname [fname ...]

    Red/Blue splicing adjustment for P200 DBSP

    positional arguments:
      fname       Final data output from DBSP_DRP to adjust the splicing of. Can be one or multiple files.

    optional arguments:
      -h, --help  show this help message and exit


In this GUI, the blue and red spectra are plotted as well as the spliced spectrum.
After using the usual Matplotlib tools to zoom and pan to the region of interest,
you can adjust the red multiplier to multiply the red flux by a constant, and
change if flux will be interpolated across detector gaps, as well as save your
changes to the same file you opened.
After you close the GUI, you will be presented with the next spectrum to adjust
the splicing of.

The default splicing that DBSP_DRP does uses a red multiplier of 1, since the
flux calibration provided by PypeIt is excellent.

.. image:: figures/splicing_adjust.png
