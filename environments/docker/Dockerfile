########################################################################
# OSMnx Dockerfile
# License: MIT, see full license in LICENSE.txt
# Web: https://github.com/gboeing/osmnx
########################################################################

FROM jupyter/base-notebook
LABEL maintainer="Geoff Boeing <boeing@usc.edu>"
LABEL url="https://github.com/gboeing/osmnx"
LABEL description="OSMnx is a Python package to retrieve, model, analyze, and visualize OpenStreetMap networks and other spatial data."

COPY requirements.txt /tmp/

# install packages in one RUN to keep image tidy
RUN mamba update --yes -c conda-forge --strict-channel-priority -n base mamba && \
    mamba install --update-all --force-reinstall --yes -c conda-forge --strict-channel-priority --file /tmp/requirements.txt && \
    rm -f -r -v /opt/conda/share/jupyter/kernels/python3 && \
    python -m ipykernel install --sys-prefix --name ox --display-name "Python (ox)" && \
    mamba clean --all --yes && \
    mamba info --all && \
    mamba list && \
    jupyter kernelspec list && \
    ipython -c "import osmnx; print('OSMnx version', osmnx.__version__)"

# copy default jupyterlab settings, then set jupyter working directory to map to mounted volume
COPY overrides.json /opt/conda/share/jupyter/lab/settings/
WORKDIR /home/jovyan/work

# set default command to launch when container is run
CMD ["jupyter", "lab", "--ip='0.0.0.0'", "--port=8888", "--no-browser", "--NotebookApp.token=''", "--NotebookApp.password=''"]
