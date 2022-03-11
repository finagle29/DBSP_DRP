#! /bin/bash

git stash -u

chmod -R g+w .

docker build -t dbsp_ql --target dbsp_ql .
docker build -t dbsp_drp --target dbsp_drp .

chmod -R g-w .

git stash pop
