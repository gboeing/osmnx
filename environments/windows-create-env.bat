SET CONDA_ROOT=%USERPROFILE%\miniforge3
SET ENV=ox
SET ENV_PATH=%CONDA_ROOT%\envs\%ENV%
SET PACKAGE=osmnx
CALL %CONDA_ROOT%\Scripts\activate.bat && ^
conda activate base && ^
conda env remove -n %ENV% --yes && ^
mamba create -n %ENV% --yes -c conda-forge --strict-channel-priority --file requirements.txt && ^
conda activate %ENV% && ^
python -m pip --python %ENV_PATH%\python.exe uninstall %PACKAGE% --yes && ^
python -m pip --python %ENV_PATH%\python.exe install -e ../. && ^
python -m pip --python %ENV_PATH%\python.exe check && ^
python -m ipykernel install --prefix %ENV_PATH% --name %ENV% --display-name "Python (%ENV%)" && ^
conda list && ^
jupyter kernelspec list && ^
ipython -c "import %PACKAGE%; print('%PACKAGE% version', %PACKAGE%.__version__)"
