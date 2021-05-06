#! /bin/bash

git stash -u

docker build -t dbsp_ql --target dbsp_ql --build-arg USER_ID=$(id -u) --build-arg GROUP_ID=$(id -g) .
docker build -t dbsp_drp --target dbsp_drp --build-arg USER_ID=$(id -u) --build-arg GROUP_ID=$(id -g) .

git stash pop
