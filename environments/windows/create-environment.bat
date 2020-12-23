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
CALL pip install https://ocf.berkeley.edu/~gboeing/share/python_igraph-0.8.3-cp39-cp39-win_amd64.whl
CALL python -m ipykernel install --sys-prefix --name ox --display-name "Python (ox)"
CALL conda clean --all --yes --quiet
CALL conda env export > environment.yml
CALL conda list
CALL jupyter kernelspec list
CALL ipython -c "import osmnx; print('OSMnx version', osmnx.__version__)"
