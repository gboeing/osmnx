# OSMnx Docker Image

## Usage

The OSMnx Docker image and container usage instructions are available on [Docker Hub](https://hub.docker.com/r/gboeing/osmnx).

## Update and distribute image

Run `docker-build.sh` in this directory. This script stops/deletes all containers/images on the system, builds the OSMnx image, exports its conda environment to a yml file in the mounted working directory, then tags and pushes the image to Docker Hub for public use.
