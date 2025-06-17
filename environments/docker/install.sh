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

wget -qO- https://astral.sh/uv/0.7.12/install.sh | sh
source $HOME/.local/bin/env
uv sync --no-cache --no-build --system --no-dev --all-extras
uv cache clean
uv pip check
uv pip list
uv pip show osmnx
python --version
python -m ipykernel install --sys-prefix --name ox --display-name "Python (ox)"
rm -f -r -v /opt/conda/share/jupyter/kernels/python3
jupyter kernelspec list
ipython -c "import osmnx; print('OSMnx version', osmnx.__version__)"
