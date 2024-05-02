SET ENV="ox"
SET PACKAGE="osmnx"
CALL %USERPROFILE%\miniforge3\Scripts\activate.bat && ^
conda deactivate && ^
mamba env remove -n %ENV% --yes && ^
mamba clean --all --yes --quiet --no-banner && ^
mamba create -c conda-forge --strict-channel-priority -n %ENV% --file "../docker/requirements.txt" --yes --no-banner && ^
conda activate %ENV% && ^
uv pip uninstall %PACKAGE% && ^
uv pip install -e ../../. && ^
python -m ipykernel install --sys-prefix --name %ENV% --display-name "Python (%ENV%)" && ^
mamba clean --all --yes --quiet --no-banner && ^
mamba env export -n %ENV% > environment.yml && ^
mamba list --no-banner && ^
jupyter kernelspec list && ^
ipython -c "import %PACKAGE%; print('%PACKAGE% version', %PACKAGE%.__version__)"
