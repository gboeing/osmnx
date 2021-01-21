#!/bin/bash
set -e
eval "$(conda shell.bash hook)"
conda deactivate
conda env remove -n ox --yes
conda clean --all --yes --quiet
conda config --set show_channel_urls true
conda config --set channel_priority strict
conda config --prepend channels conda-forge
mamba create -n ox --file "../docker/requirements.txt" --yes
eval "$(conda shell.bash hook)"
conda activate ox
pip uninstall osmnx --yes
pip install -e ../../.
python -m ipykernel install --sys-prefix --name ox --display-name "Python (ox)"
conda clean --all --yes --quiet
conda env export -n ox > environment.yml
conda list
jupyter kernelspec list
ipython -c "import osmnx; print('OSMnx version', osmnx.__version__)"
