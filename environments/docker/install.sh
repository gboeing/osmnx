#!/bin/bash
set -euo pipefail

# rasterio doesn't provide linux/arm64 wheels, so if the target platform is
# linux/arm64, don't install this optional dependency (attempting to build
# it rather than install the wheel will also fail).
# see https://github.com/rasterio/rasterio-wheels/issues/69
if [[ "$TARGETPLATFORM" == "linux/arm64" ]]
then
  NOEXTRA="--no-extra=raster --no-extra=all"
else
  NOEXTRA=""
fi

# install all requirements into the existing system environment
uv export --no-cache --no-build --all-extras $NOEXTRA --group examples > requirements-temp.txt
uv pip install --no-cache --no-build --system --compile-bytecode --strict -r requirements-temp.txt
rm -f requirements-temp.txt
uv cache clean
python --version
ipython -c "import osmnx; print('OSMnx version', osmnx.__version__)"
