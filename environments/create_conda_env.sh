#!/bin/bash
set -euo pipefail
echo "Run conda deactivate before running this script."
ENV=ox
ENV_PATH=$(conda info --base)/envs/$ENV
PACKAGE=osmnx
uv --version
eval "$(conda shell.bash hook)"
conda deactivate
conda env remove --yes -n $ENV || true
conda create --yes -c conda-forge --strict-channel-priority -n $ENV python
eval "$(conda shell.bash hook)"
conda activate $ENV
uv export --no-build --all-extras --all-groups > ./environments/requirements-temp.txt
uv pip install --no-build --strict -r ./environments/requirements-temp.txt
rm -f ./environments/requirements-temp.txt
python -m pip --python "$ENV_PATH" uninstall $PACKAGE --yes
python -m pip --python "$ENV_PATH" install -e .
python -m ipykernel install --prefix "$ENV_PATH" --name $ENV --display-name "Python ($ENV)"
conda list -n $ENV
python -m pip --python "$ENV_PATH" check
jupyter kernelspec list
ipython -c "import $PACKAGE; print('$PACKAGE version', $PACKAGE.__version__)"
