#!/bin/bash
set -e
ENV=ox
PACKAGE=osmnx
eval "$(conda shell.bash hook)"
conda deactivate
mamba env remove -n $ENV --yes
mamba clean --all --yes --quiet
mamba create -c conda-forge --strict-channel-priority -n $ENV --file "../docker/requirements.txt" --yes
eval "$(conda shell.bash hook)"
conda activate $ENV
pip uninstall $PACKAGE --yes
pip install -e ../../.
python -m ipykernel install --sys-prefix --name $ENV --display-name "Python ($ENV)"
mamba clean --all --yes --quiet
mamba env export -n $ENV > environment.yml
mamba list
jupyter kernelspec list
ipython -c "import $PACKAGE; print('$PACKAGE version', $PACKAGE.__version__)"
