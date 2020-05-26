# OSMnx Docker Image

## Usage

The OSMnx docker container image and usage instructions are available on [docker hub](https://hub.docker.com/r/gboeing/osmnx).

## Development

### Build an image from the dockerfile

```
docker build -t gboeing/osmnx .
```

### Run bash in this container and export final conda environment to a yml file

```
docker run --rm -it -v "$PWD":/home/jovyan/work gboeing/osmnx /bin/bash
conda env export -n base > /home/jovyan/work/environment.yml
```

### Push the built image to hub so others can pull/run it

```
docker login
docker tag gboeing/osmnx gboeing/osmnx:v0.0.0
docker push gboeing/osmnx
```

### Stop/delete all local docker containers/images

```
docker stop $(docker ps -aq)
docker rm $(docker ps -aq)
docker rmi $(docker images -q) --force
```
