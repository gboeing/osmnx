CALL conda update --strict-channel-priority conda -n base --yes
CALL conda deactivate
CALL conda env remove -n ox --yes
CALL conda create -n ox -c conda-forge --file "../../requirements.txt" --file "../../requirements-dev.txt" --file "requirements-more.txt" --yes
CALL conda activate ox
CALL python -m ipykernel install --user --name ox --display-name "Python (ox)"
CALL jupyter labextension install @jupyter-widgets/jupyterlab-manager
CALL jupyter labextension install jupyter-leaflet
CALL pip install https://ocf.berkeley.edu/~gboeing/share/python_igraph-0.7.1.post6-cp37-cp37m-win_amd64.whl
CALL conda clean --all --yes
CALL conda env export > environment-windows.yml
CALL conda list
CALL jupyter kernelspec list
CALL python -c "import osmnx; print(osmnx.__version__)"
