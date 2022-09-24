#!/bin/bash
DOCKERUSER=gboeing
PACKAGE=osmnx

# login and remove any existing containers or images
docker login
docker stop $(docker ps -aq)
docker rm $(docker ps -aq)
docker rmi $(docker images -q) --force

# build the image and export the conda env to yml
set -e
docker build -t $DOCKERUSER/$PACKAGE .
docker run --rm -v "$PWD":/home/jovyan/work $DOCKERUSER/$PACKAGE:latest /bin/bash -c "conda env export -n base > /home/jovyan/work/environment.yml"

# get the package version, tag the image with it, then push to hub
VERSION=$(docker run --rm $DOCKERUSER/$PACKAGE:latest /bin/bash -c "python -c \"import $PACKAGE; print($PACKAGE.__version__)\"")
echo "$PACKAGE version $VERSION"
docker tag $DOCKERUSER/$PACKAGE $DOCKERUSER/$PACKAGE:v$VERSION
docker push -a $DOCKERUSER/$PACKAGE
