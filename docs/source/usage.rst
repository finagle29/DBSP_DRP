**************
Using DBSP_DRP
**************

After :doc:`installing`, you are ready to reduce some DBSP data!

The basic usage of DBSP_DRP is as follows:

.. code-block:: console

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
