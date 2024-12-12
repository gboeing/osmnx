#!/bin/bash
set -euo pipefail
echo "Run conda deactivate before running this script."
ENV=ox
ENV_PATH=$(conda info --base)/envs/$ENV
PACKAGE=osmnx
eval "$(conda shell.bash hook)"
conda activate base
conda env remove --yes -n $ENV || true
conda update --yes -c conda-forge --strict-channel-priority -n base conda mamba
mamba create --yes -c conda-forge --strict-channel-priority -n $ENV --file requirements.txt
conda activate $ENV
python -m pip --python $ENV_PATH uninstall $PACKAGE --yes
python -m pip --python $ENV_PATH install -e ../.
python -m ipykernel install --prefix $ENV_PATH --name $ENV --display-name "Python ($ENV)"
conda list -n $ENV
python -m pip --python $ENV_PATH check
jupyter kernelspec list
ipython -c "import $PACKAGE; print('$PACKAGE version', $PACKAGE.__version__)"
