#!/bin/bash
set -euo pipefail

# rasterio doesn't have linux/arm64 wheels, so if this platform target is
# linux/arm64 then don't install this optional dependency (attempting to
# build it rather than install the wheel also fails).
if [[ "$TARGETPLATFORM" == "linux/arm64" ]]
then
  NOEXTRA="--no-extra=raster"
else
  NOEXTRA=""
fi

wget -qO- https://astral.sh/uv/0.7.9/install.sh | sh
source $HOME/.local/bin/env
uv export --no-cache --no-build --no-dev --all-extras $NOEXTRA > reqs.txt
uv pip install --no-cache --no-build --system --strict -r reqs.txt -r requirements-extras.txt
uv cache clean
uv pip check
uv pip list
uv pip show osmnx
python --version
python -m ipykernel install --sys-prefix --name ox --display-name "Python (ox)"
rm -f -r -v /opt/conda/share/jupyter/kernels/python3
jupyter kernelspec list
ipython -c "import osmnx; print('OSMnx version', osmnx.__version__)"
