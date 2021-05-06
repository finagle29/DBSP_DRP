FROM condaforge/miniforge3 as dbsp_ql
LABEL Author, Milan Roberson

ARG USER_ID
ARG GROUP_ID

RUN groupadd --non-unique --gid $GROUP_ID user
RUN adduser --disabled-password --gecos '' --uid $USER_ID --gid $GROUP_ID user

# needed for Qt/PySide2
RUN apt-get update && apt-get install -y \
    libgl1-mesa-dev \
    libxrender1 \
    xauth


ENV WORKDIR /workdir
WORKDIR $WORKDIR

# set up conda environment
COPY environment.yml $WORKDIR
RUN conda update --name base conda
USER user
RUN conda info
RUN conda env create --file environment.yml

# copy repo over and install it
COPY . $WORKDIR/DBSP_DRP

RUN /bin/bash -c ". activate dbsp_drp && \
    pip install DBSP_DRP/"

# make dbsp_drp the default conda environment
RUN echo 'conda activate dbsp_drp' >> /home/user/.bashrc

CMD [ "/bin/bash" ]

FROM dbsp_ql as dbsp_drp
USER root
RUN apt-get update && \
    apt-get install -y curl
USER user
RUN /bin/bash -c ". activate dbsp_drp && \
    DBSP_DRP/bin/download_tellfile"
CMD [ "/bin/bash" ]
