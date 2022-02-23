*********************
DBSP_DRP Docker Image
*********************

Docker images are provided to allow for portable ``DBSP_DRP`` installations, or for
easily maintaining an on-site ``DBSP_DRP`` installation. The Docker image ``finagle29/dbsp_drp``
contains a full installation of ``DBSP_DRP`` (i.e. including the telluric grid), and the Docker
image ``finagle29/dbsp_ql`` contains a minimal installation of ``DBSP_DRP`` suitable for running
the quicklook script. By appending a tag to the Docker image name, you can control which version
of ``DBSP_DRP`` the downloaded Docker image contains:

.. table::
    :widths: 10 30 30

    ========== ========================================== ============================
    Tag        Version                                    Example Name
    ========== ========================================== ============================
    (blank)    same as ``latest`` tag                     ``finagle29/dbsp_ql``
    ``latest`` latest released version (currently 1.0.0)  ``finagle29/dbsp_ql:latest``
    ``main``   latest commit on ``main`` branch           ``finagle29/dbsp_ql:main``
    ``edge``   latest commit on ``develop`` branch        ``finagle29/dbsp_ql:edge``
    ``1.0.0``  version 1.0.0                              ``finagle29/dbsp_ql:1.0.0``
    ========== ========================================== ============================


The shell script ``run_docker.sh`` is provided for ease of running ``DBSP_DRP`` in a Docker
container. After installing Docker on your computer, with the ``run_docker.sh`` script in
your working directory, you can start a Docker container with ``DBSP_DRP`` pre-installed
using the following command-line invocation:

.. code-block :: console

    $ ./run_docker.sh -v /data/DBSP_20210919:/workdir/data finagle29/dbsp_drp:0.9.0

This will start the Docker container in the terminal window, and allow you to access the
contents of the directory ``/data/DBSP_20210919`` on your computer in the directory
``/workdir/data`` in the Docker container.
