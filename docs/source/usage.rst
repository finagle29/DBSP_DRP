**************
Using DBSP_DRP
**************

After :doc:`installing`, you are ready to reduce some DBSP data!

The basic usage of DBSP_DRP is as follows:

.. code-block :: console

    $ dbsp_reduce -r /path/to/data/DBSP_YYYYMMDD/ -d /output/path/DBSP_YYYYMMDD_redux

When reducing a night of data on a machine with multiple CPU cores, it is highly
recommended to add the ``-j N`` flag to use N jobs for the telluric correction.

*Note:* the default setting is to open a GUI that allows you to verify that the
header data is correct for all of the files. Add the option ``-i`` or
``--no-interactive`` to turn this behavior off.

If you want to only reduce the red arm or the blue arm, add the flag ``--arm red``
or ``-a blue`` for short.

If you want lots of intermediate debugging plots displayed to the screen, add the
flag ``--debug``.

If you want to check that your target traces have been correctly identified, and
manually select them if they were missed, add the flag ``--manual-extraction`` or
``-m`` for short. After the first round of reducing the red and/or blue sides,
``dbsp_reduce`` will open a matplotlib window to display the sky-subtracted spectra,
along with any object traces that were automatically detected. Using the left and
right arrow keys, you can cycle through the spectra. If your science target was not
automatically detected, you can zoom in on the trace using the standard matplotlib
zoom tool and then with your mouse over the trace, press the m key to mark that
trace. If you make a mistake, you can press d with your mouse over a
previously-marked trace to delete the trace. Once you are satisfied with the
identification of traces on all of the spectra, close the matplotlib window to repeat
this process for the other arm. After closing the second window, the manually marked
spectra will be re-reduced.

If you want to fine-tune the reduction parameters, create a parameter file like so:

.. code-block :: cfg

    # params.cfg
    [blue]
    # PypeIt reduction parameters to control the blue side reduction
    [reduce]
        [[skysub]]
            no_local_sky = True
        [[extraction]]
            use_2dmodel_mask = False
    [red]
    # PypeIt parameters to control the red side reduction
    [reduce]
        [[skysub]]
            no_local_sky = True
        [[extraction]]
            use_2dmodel_mask = False

and use the option ``-p params.cfg`` or its longer form ``--parameter-file params.cfg``.

See `PypeIt Parameters <https://pypeit.readthedocs.io/en/stable/pypeit_par.html>`_ for more
details on the full set of available parameters.

Manual Aperture Selection and Sky Selection
*******************************************
Using the ``-m`` flag, after marking traces manually, for each frame you marked a trace on,
a second GUI will pop up showing the collapsed flux, along with any automatically identified
traces (in orange) and your manually placed traces (in blue) along with the FWHM of each shaded
in a lighter color. In this GUI, you can left click and drag your manual traces (in blue) to
adjust their position. Additionally, you can right click and drag the shaded FWHM regions to
adjust their extent. To mark background regions, press b and then left click and drag to mark
background regions by shading them in gray. You can delete a background region by holding your
mouse over the shaded background regions and press d to delete. Once you are finished adjusting
manual traces/FWHMs and marking background regions, close the GUI to be shown the GUI for the
next object you marked a manual trace on.

Make sure to select background regions on either side of your target for the best sky subtraction.

If you are dealing with faint sources, it is a good idea to re-mark in blue any orange
(automatically identified) traces in case parameter changes lose these objects.
