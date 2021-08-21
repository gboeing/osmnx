CALL %USERPROFILE%\miniconda3\Scripts\activate.bat
CALL conda deactivate
CALL conda env remove -n ox --yes
CALL conda clean --all --yes --quiet
CALL conda config --set show_channel_urls true
CALL conda config --set channel_priority strict
CALL conda config --prepend channels conda-forge
CALL mamba create -n ox --file "../docker/requirements.txt" --yes
CALL conda activate ox
CALL pip uninstall osmnx --yes
CALL pip install -e ../../.
CALL python -m ipykernel install --sys-prefix --name ox --display-name "Python (ox)"
CALL conda clean --all --yes --quiet
CALL conda env export > environment.yml
CALL conda list
CALL jupyter kernelspec list
CALL ipython -c "import osmnx; print('OSMnx version', osmnx.__version__)"
