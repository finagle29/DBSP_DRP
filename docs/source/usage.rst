**************
Using DBSP_DRP
**************

After :doc:`installing`, you are ready to reduce some DBSP data!

.. code-block :: console

    $ dbsp_reduce --help
    usage: dbsp_reduce [-h] [-i] -r ROOT -d OUTPUT_PATH [-a {red,blue}] [-m]
                   [--debug] [-j N] [-p PARAMETER_FILE] [-t] [-c]
                   [--splicing-interpolate-gaps]

    Automatic Data Reduction Pipeline for P200 DBSP

    optional arguments:
        -h, --help            show this help message and exit
        -i, --no-interactive  Interactive file-checking?
        -r ROOT, --root ROOT  File path+root, e.g. /data/DBSP_20200127
        -d OUTPUT_PATH, --output_path OUTPUT_PATH
                                Path to top-level output directory.  Default is the current working directory.
        -a {red,blue}, --arm {red,blue}
                                [red, blue] to only reduce one arm (null splicing)
        -m, --manual-extraction
                                manual extraction
        --debug               debug
        -j N, --jobs N        Number of processes to use
        -p PARAMETER_FILE, --parameter-file PARAMETER_FILE
                                Path to parameter file. The parameter file should be formatted as follows:

                                [blue]
                                ** PypeIt parameters for the blue side goes here **
                                [red]
                                ** PypeIt parameters for the red side goes here **
                                EOF

                                The [red/blue] parameter blocks are optional, and their order does not matter.
        -t, --skip-telluric   Skip telluric correction
        -c, --null-coadd      Don't coadd consecutive exposures of the same target.
                                By default consective exposures will be coadded.
        --splicing-interpolate-gaps
                                Use this option to linearly interpolate across large gaps
                                in the spectrum during splicing. The default behavior is to
                                only use data from one detector in these gaps, which results
                                in a slightly noisier spliced spectrum.


The basic usage of DBSP_DRP is as follows:

.. code-block :: console

    $ dbsp_reduce -r /path/to/data/DBSP_YYYYMMDD/ -d /output/path/DBSP_YYYYMMDD_redux

When reducing a night of data on a machine with multiple CPU cores, it is highly
recommended to add the ``-j N`` flag to use N jobs for the telluric correction.

*Note:* the default setting is to open a GUI that allows you to verify that the
header data is correct for all of the files. Add the option ``-i`` or
``--no-interactive`` to turn this behavior off.

If you only want to reduce a subset of your data, you can select unwanted files in the
header validation table GUI and right click to delete them from the current reduction run.
Be sure to always keep the standard stars you need for fluxing in the data reduction.

If you want to only reduce the red arm or the blue arm, add the flag ``--arm red``
or ``-a blue`` for short.

If you want lots of intermediate debugging plots displayed to the screen, add the
flag ``--debug``.

If you want to check that your target traces have been correctly identified, and
manually select them if they were missed, add the flag ``--manual-extraction`` or
``-m`` for short. Check out :doc:`manual_extraction` for more details.

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
