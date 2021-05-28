#! /bin/bash

# exit immediately on non-zero exit code
#set -e

# add group and user
groupadd --non-unique --force --gid ${GROUP_ID} ${GROUP_NAME}
adduser --disabled-password --gecos '' --uid ${USER_ID} --gid ${GROUP_ID} ${USER_NAME}
usermod -a -G dbsp ${USER_NAME}
# See Dockerfile for more details about how to give user sudo powers.
# usermod -a -G sudo ${USER_NAME}

# make dbsp_drp default conda environment
echo 'conda activate dbsp_drp' >> /home/${USER_NAME}/.bashrc

ln -s /.Xauthority /home/${USER_NAME}/.Xauthority
chown ${USER_NAME}:dbsp /.Xauthority

su ${USER_NAME}
