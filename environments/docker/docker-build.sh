#!/bin/bash
docker stop $(docker ps -aq)
docker rm $(docker ps -aq)
docker rmi $(docker images -q) --force
docker build -t gboeing/osmnx .
docker run --rm -v "$PWD":/home/jovyan/work gboeing/osmnx:latest /bin/bash -c "conda env export -n base > /home/jovyan/work/environment.yml"
