#!/bin/bash
set -e
docker login
docker buildx build --no-cache --pull --load --platform=linux/arm64 -f environments/docker/Dockerfile -t gboeing/osmnx:test .
IMPORTED_VERSION=$(docker run --rm gboeing/osmnx:test /bin/bash -c "ipython -c \"import osmnx; print(osmnx.__version__)\"")
echo "Imported $IMPORTED_VERSION"
