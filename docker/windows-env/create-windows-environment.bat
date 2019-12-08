CALL conda deactivate
CALL conda clean --all --yes
CALL conda config --set show_channel_urls true
CALL conda config --set channel_priority strict
CALL conda config --prepend channels conda-forge
CALL conda update conda -n base --yes
CALL conda env remove -n ox --yes
CALL conda create -n ox --file "../requirements.txt" --yes
CALL conda activate ox
CALL pip uninstall osmnx --yes
CALL pip install -e ../../.
CALL pip install https://ocf.berkeley.edu/~gboeing/share/python_igraph-0.7.1.post6-cp37-cp37m-win_amd64.whl
CALL python -m ipykernel install --user --name ox --display-name "Python (ox)"
CALL conda clean --all --yes
CALL conda env export > environment-windows.yml
CALL conda list
CALL jupyter kernelspec list
CALL python -c "import osmnx; print(osmnx.__version__)"
