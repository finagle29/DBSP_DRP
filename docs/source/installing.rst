*******************
Installing DBSP_DRP
*******************

Conda is the recommended pacakage manager used to install ``DBSP_DRP``.

From Source
###########

.. code-block :: console

    $ git clone https://github.com/finagle29/DBSP_DRP.git
    $ cd DBSP_DRP
    $ conda env create -f environment.yml
    $ conda activate dbsp_drp
    $ pip install -e .

This performs an editable install, which allows you to make modifications to the code and immediately see their effects.
Importantly, this can be used in combination with ``git`` branches to test features in development.

Using pip
#########

First download the provided `environment.yml file <https://raw.githubusercontent.com/finagle29/DBSP_DRP/master/environment.yml>`__

Now use the environment.yml file to create a conda environment with the required dependencies.

.. code-block :: console

    $ cd /path/to/Downloads
    $ conda env create -f environment.yml
    $ conda activate dbsp_drp

Now use ``pip`` to install DBSP_DRP

.. code-block :: console

    $ pip install git+https://github.com/finagle29/DBSP_DRP.git

************
Post-Install
************

The telluric correction code provided by PypeIt relies on a large (5 GB) atmospheric model file
(TellFit_Lick_3100_11100_R10000.fits), which can be downloaded
`here <https://drive.google.com/drive/folders/1FFRWjUZ58HiDuDD33MYqBzMWDQanBRRy>`__
and must be installed into the ``pypeit/data/telluric/atm_grids`` directory of your PypeIt installation.

To determine the location of your PypeIt installation, open the Python interpreter and run

.. code-block :: Python

    >>> import pypeit
    >>> import os
    >>> print(os.path.dirname(pypeit.__file__))
    /Users/me/anaconda3/envs/dbsp_drp/lib/python3.7/site-packages/pypeit

An easier alternative is to download and run `this script <https://raw.githubusercontent.com/finagle29/DBSP_DRP/master/bin/download_tellfile>`__,
which will perform the download and install it into the current PypeIt installation.

.. code-block :: console

    $ cd /path/to/Downloads
    $ chmod +x download_tellfile
    $ ./download_tellfile

If you have multiple PypeIt installations on the same machine, you can create a hard link from the one PypeIt installation to the others so you can reuse the atmospheric model file.

.. code-block :: console

    $ ln /path/to/stable/pypeit/data/telluric/atm_grids/TellFit_Lick_3100_11100_R10000.fits /path/to/other/pypeit/data/telluric/atm_grids/TellFit_Lick_3100_11100_R10000.fits

Make sure to use the same filename in both PypeIt installations.
If you're not sure where your PypeIt installations are, run the previous Python snippet in each ``conda`` or ``venv`` environment you want to use ``DBSP_DRP`` in.

Testing your installation
#########################

Make sure your PypeIt installation was successful

.. code-block:: console

    $ run_pypeit -h

Run some built-in tests for DBSP_DRP, including verification that the quicklook script works

.. code-block:: console

    $ cd /path/to/DBSP_DRP
    $ pytest .
