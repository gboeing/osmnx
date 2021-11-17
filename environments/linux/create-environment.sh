#!/bin/bash
set -e
eval "$(conda shell.bash hook)"
conda deactivate
mamba env remove -n ox --yes
mamba clean --all --yes --quiet
mamba create -c conda-forge --strict-channel-priority -n ox --file "../docker/requirements.txt" --yes
eval "$(conda shell.bash hook)"
conda activate ox
pip uninstall osmnx --yes
pip install -e ../../.
python -m ipykernel install --sys-prefix --name ox --display-name "Python (ox)"
mamba clean --all --yes --quiet
mamba env export -n ox > environment.yml
mamba list
jupyter kernelspec list
ipython -c "import osmnx; print('OSMnx version', osmnx.__version__)"
