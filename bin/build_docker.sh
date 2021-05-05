#! /bin/bash

git stash -u

docker build -t dbsp_ql --target dbsp_ql .
docker build -t dbsp_drp --target dbsp_drp .

git stash pop
