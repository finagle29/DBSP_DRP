#! /bin/bash
# TODO: test argument passing to the docker run command

last_arg_msg="The last argument to run_docker.sh must be either dbsp_ql or dbsp_drp"

if [[ $# == 0 ]]; then
    echo "Usage: /path/to/run_docker.sh [docker arguments] [dbsp_ql OR dbsp_drp]"
    echo $last_arg_msg
    exit 1
fi

for last in $@; do :; done

if [[ $last != "dbsp_drp" && $last != "dbsp_ql" ]]; then
    echo $last_arg_msg
    exit 1
fi

if [[ $OSTYPE == "darwin"* ]]; then
    IP=`ifconfig en0 | grep 'inet ' | cut -d " " -f2`
    xhost +$IP
    docker run -t -i --net=host -e DISPLAY=$IP:0 -v /tmp/.X11-unix:/tmp/.X11-unix "$@"
elif [[ $OSTYPE == "linux"* ]]; then
# elif [[ `cat /etc/issue` == "Ubuntu*"]]; then
# not sure if this /needs/ to be ubuntu
    docker run -t -i --net=host -e DISPLAY --volume="$HOME/.Xauthority:/root/.Xauthority:rw" "$@"
elif [[ $OSTYPE == "win"* ]]; then
    echo "You'll have to try to launch Docker on your own, sorry."
fi
