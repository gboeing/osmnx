#!/bin/bash
eval "$(conda shell.bash hook)"
conda deactivate
conda env remove -n ox --yes
conda clean --all --yes
conda config --set show_channel_urls true
conda config --set channel_priority strict
conda config --prepend channels conda-forge
conda update conda -n base --yes
conda create -n ox --file "../requirements.txt" python-igraph --yes
eval "$(conda shell.bash hook)"
conda activate ox
pip uninstall osmnx --yes
pip install -e ../../.
python -m ipykernel install --user --name ox --display-name "Python (ox)"
conda clean --all --yes
conda env export -n ox > environment-linux.yml
conda list
jupyter kernelspec list
ipython -c "import osmnx; print(osmnx.__version__)"
