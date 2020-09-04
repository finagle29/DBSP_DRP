**************
Using DBSP_DRP
**************

After :doc:`installing`, you are ready to reduce some DBSP data!

The basic usage of DBSP_DRP is as follows:

.. code-block:: console

    $ dbsp_reduce -r /path/to/data/DBSP_YYYYMMDD/ -d /output/path/DBSP_YYYYMMDD_redux

When reducing a night of data on a machine with multiple CPU cores, it is highly recommended to add the ``-j N`` flag to use N jobs for the telluric correction.

*Note:* the default setting is to open a GUI that allows you to verify that the header data is correct for all of the files.
Add the option ``-i`` or ``--no-interactive`` to turn this behavior off.

If you want to only reduce the red arm or the blue arm, add the flag ``--arm red`` or ``-a blue`` for short.

If you want lots of intermediate debugging plots displayed to the screen, add the flag ``--debug``.
