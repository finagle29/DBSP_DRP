#! /bin/bash

last_arg_msg="The last argument to run_docker.sh must be either dbsp_ql or dbsp_drp"

if [[ $# == 0 ]]; then
    echo "Usage: /path/to/run_docker.sh [docker arguments] [dbsp_ql OR dbsp_drp]"
    echo "usually you want [docker arguments] to be something like -v /path/to/your/raw/data:/workdir/data"
    echo "to make /path/to/your/raw/data available in the Docker container in /workdir/data"
    echo $last_arg_msg
    exit 1
fi

for last in $@; do :; done

if [[ $last != "dbsp_drp" && $last != "dbsp_ql" ]]; then
    echo $last_arg_msg
    exit 1
fi

user_name=$(id -un)
group_name=$(id -gn)

if [[ $OSTYPE == "darwin"* ]]; then
    IP=`ifconfig en0 | grep 'inet ' | cut -d " " -f2`
    xhost +$IP
    docker run --rm -t -i --net=host -e DISPLAY=$IP:0 -v /tmp/.X11-unix:/tmp/.X11-unix --volume="$HOME/.Xauthority:/.Xauthority:rw" -e "USER_ID=$(id -u)" -e "GROUP_ID=$(id -g)" -e "USER_NAME=$user_name" -e "GROUP_NAME=$group_name" "$@"
elif [[ $OSTYPE == "linux"* ]]; then
# elif [[ `cat /etc/issue` == "Ubuntu*"]]; then
# not sure if this /needs/ to be ubuntu
    docker run --rm -t -i --net=host -e DISPLAY --volume="$HOME/.Xauthority:/.Xauthority:rw" -e "USER_ID=$(id -u)" -e "GROUP_ID=$(id -g)" -e "USER_NAME=$user_name" -e "GROUP_NAME=$group_name" "$@"
elif [[ $OSTYPE == "win"* ]]; then
    echo "You'll have to try to launch Docker on your own, sorry."
fi
