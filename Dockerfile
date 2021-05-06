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

# make pip install run in dbsp_drp
SHELL [ "conda", "run", "--name", "dbsp_drp", "/bin/bash", "-c" ]

# copy repo over and install it
COPY . $WORKDIR
RUN pip install -e .

# make dbsp_drp the default conda environment
RUN echo 'conda activate dbsp_drp' >> /root/.bashrc

CMD [ "/bin/bash" ]

FROM dbsp_ql as dbsp_drp
USER root
RUN apt-get update && \
    apt-get install -y curl
USER user
RUN bin/download_tellfile
CMD [ "/bin/bash" ]
