#!/bin/bash
DOCKERUSER=gboeing
PACKAGE=osmnx

# get the version number for the tag
read -r -p "What is the $PACKAGE version number? " VERSION
TAG=v$VERSION
echo "Image will be tagged $DOCKERUSER/$PACKAGE:$TAG"
read -r -n 1 -p "Ready to proceed? (y/n) " INPUT
echo ""
if [[ $INPUT != "Y" && $INPUT != "y" ]]; then
    exit 1
fi

# remove any existing containers or images
docker stop $(docker ps -aq)
docker rm $(docker ps -aq)
docker rmi $(docker images -q) --force

# login then multi-platform buildx the image and push it to the hub
set -e
docker login
docker buildx build --platform=linux/amd64,linux/arm64 --push -t $DOCKERUSER/$PACKAGE:$TAG .

# import package and print its version as a test, then export conda env to yml
IMPORTED_VERSION=$(docker run --rm $DOCKERUSER/$PACKAGE:$TAG /bin/bash -c "python -c \"import $PACKAGE; print($PACKAGE.__version__)\"")
echo "$Imported PACKAGE version $IMPORTED_VERSION, expected $VERSION"
docker run --rm -v "$PWD":/home/jovyan/work $DOCKERUSER/$PACKAGE:$TAG /bin/bash -c "conda env export -n base > /home/jovyan/work/environment.yml"
