#!/bin/bash
set -e
ENV=ox
ENV_PATH=$(conda info --base)/envs/$ENV
PACKAGE=osmnx
eval "$(conda shell.bash hook)"
conda activate base
conda env remove -n $ENV --yes
mamba create -n $ENV --yes -c conda-forge --strict-channel-priority --file "../docker/requirements.txt"
eval "$(conda shell.bash hook)"
conda activate $ENV
python -m pip --python $ENV_PATH uninstall $PACKAGE --yes
python -m pip --python $ENV_PATH install -e ../../.
python -m ipykernel install --prefix $ENV_PATH --name $ENV --display-name "Python ($ENV)"
conda env export -n $ENV > environment.yml
conda list -n $ENV
jupyter kernelspec list
ipython -c "import $PACKAGE; print('$PACKAGE version', $PACKAGE.__version__)"
