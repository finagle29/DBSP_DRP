FROM condaforge/miniforge3 as dbsp_ql
LABEL Author, Milan Roberson

# create dbsp group
RUN groupadd --gid 10001 dbsp && \
# needed for Qt/PySide2
    apt-get update && apt-get install -y \
    libgl1-mesa-dev \
    libxrender1 \
    xauth
# Uncomment and add to above RUN command, and also uncomment parts of entrypoint.sh
# to give user sudo powers. This should not be necessary, as normal things like
# vim are available via conda.
#    sudo && \
#    echo '%sudo ALL=(ALL) NOPASSWD:ALL' >> \
#        /etc/sudoers

ENV WORKDIR /workdir
WORKDIR $WORKDIR

# set up conda environment
COPY --chown=root:dbsp environment.yml $WORKDIR
RUN conda update --name base conda && \
    conda env create --file environment.yml

# copy repo over and install it
COPY --chown=root:dbsp . $WORKDIR/DBSP_DRP

RUN --mount=source=.git,target=$WORKDIR/DBSP_DRP/.git,type=bind \
    /bin/bash -c ". activate dbsp_drp && \
    pip install --use-feature=in-tree-build DBSP_DRP/" && \
# give dbsp group rwx access to conda installation
    chgrp -R dbsp /opt/conda && \
    chmod 770 -R /opt/conda

CMD [ "/bin/bash" ]

ENTRYPOINT [ "DBSP_DRP/bin/entrypoint.sh" ]

FROM dbsp_ql as dbsp_drp

ARG TELLFILE=""

ENV FETCH_TELLFILE=${TELLFILE:+"1"}
ENV TELLFILE=${TELLFILE:+"foobar"}

RUN if [ "${FETCH_TELLFILE}" = "1" ]; \
    then apt-get update && \
    apt-get install -y curl && \
     /bin/bash -c ". activate dbsp_drp && \
    DBSP_DRP/bin/download_tellfile"; \
    fi;

COPY Dockerfile ${TELLFILE}* /opt/conda/envs/dbsp/lib/python*/site-packages/pypeit/data/telluric/atm_grids/

CMD [ "/bin/bash" ]

ENTRYPOINT [ "DBSP_DRP/bin/entrypoint.sh" ]
