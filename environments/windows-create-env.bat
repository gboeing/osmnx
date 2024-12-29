SET CONDA_ROOT=%USERPROFILE%\miniforge3
SET ENV=ox
SET ENV_PATH=%CONDA_ROOT%\envs\%ENV%
SET PACKAGE=osmnx
CALL %CONDA_ROOT%\Scripts\activate.bat && ^
conda activate base && ^
conda env remove --yes -n %ENV% && ^
mamba create --yes -c conda-forge --strict-channel-priority -n %ENV% --file requirements\requirements-all.txt && ^
conda activate %ENV% && ^
python -m pip --python %ENV_PATH%\python.exe uninstall %PACKAGE% --yes && ^
python -m pip --python %ENV_PATH%\python.exe install -e ../. && ^
python -m ipykernel install --prefix %ENV_PATH% --name %ENV% --display-name "Python (%ENV%)" && ^
conda list && ^
python -m pip --python %ENV_PATH%\python.exe check && ^
jupyter kernelspec list && ^
ipython -c "import %PACKAGE%; print('%PACKAGE% version', %PACKAGE%.__version__)"
