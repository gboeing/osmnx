CALL %USERPROFILE%\mambaforge\Scripts\activate.bat
CALL mamba deactivate
CALL mamba env remove -n ox --yes
CALL mamba clean --all --yes --quiet
CALL mamba create -c conda-forge --strict-channel-priority -n ox --file "../docker/requirements.txt" --yes
CALL conda activate ox
CALL pip uninstall osmnx --yes
CALL pip install -e ../../.
CALL python -m ipykernel install --sys-prefix --name ox --display-name "Python (ox)"
CALL mamba clean --all --yes --quiet
CALL mamba env export > environment.yml
CALL mamba list
CALL jupyter kernelspec list
CALL ipython -c "import osmnx; print('OSMnx version', osmnx.__version__)"
