*******************
Installing DBSP_DRP
*******************

Conda is the recommended pacakage manager used to install DBSP_DRP.

From Source
###########

First download the `DBSP_DRP source <https://github.com/finagle29/DBSP_DRP/archive/master.zip>`__ and unzip it.

Now use the environment.yml file to create a conda environment with the required dependencies.

.. code-block:: console

    $ cd /path/to/DBSP_DRP
    $ conda env create -f environment.yml
    $ conda activate dbsp_drp

And use the standard ``setup.py`` invocation to install DBSP_DRP.

.. code-block:: console

    $ python setup.py install


Using pip
#########

First download the provided `environment.yml file <https://raw.githubusercontent.com/finagle29/DBSP_DRP/master/environment.yml>`__

Now use the environment.yml file to create a conda environment with the required dependencies.

.. code-block:: console

    $ cd /path/to/Downloads
    $ conda env create -f environment.yml
    $ conda activate dbsp_drp

Now use ``pip`` to install DBSP_DRP

.. code-block:: console

    $ pip install https://github.com/finagle29/DBSP_DRP/archive/master.tar.gz

*************************
Testing your installation
*************************

Make sure your PypeIt installation was successful

.. code-block:: console

    $ run_pypeit -h
