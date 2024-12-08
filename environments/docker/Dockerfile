FROM jupyter/base-notebook
LABEL maintainer="Geoff Boeing <boeing@usc.edu>"
LABEL url="https://osmnx.readthedocs.io"
LABEL description="OSMnx is a Python package to easily download, model, analyze, and visualize street networks and other geospatial features from OpenStreetMap."

# copy the package files needed for installation
COPY --chmod=0755 ./environments/docker/requirements.txt ./osmnx/
COPY --chmod=0755 ./osmnx/ ./osmnx/osmnx/
COPY --chmod=0755 ./pyproject.toml ./osmnx/
COPY --chmod=0755 ./README.md ./osmnx/

# install and configuration everything in one RUN to keep image tidy
RUN conda update --yes -c conda-forge --strict-channel-priority -n base mamba && \
    mamba install --update-all --force-reinstall --yes -c conda-forge --strict-channel-priority --file ./osmnx/requirements.txt && \
    python -m pip install --no-cache-dir -e ./osmnx/ && \
    python -m ipykernel install --sys-prefix --name ox --display-name "Python (ox)" && \
    rm -f -r -v /opt/conda/share/jupyter/kernels/python3 && \
    python -m pip cache purge && \
    conda clean --all --yes --quiet && \
    conda info --all && \
    conda list && \
    jupyter kernelspec list && \
    ipython -c "import osmnx; print('OSMnx version', osmnx.__version__)"

# copy default jupyterlab settings, then set jupyter working directory to map to mounted volume
COPY --chmod=0755 ./environments/docker/overrides.json /opt/conda/share/jupyter/lab/settings/
WORKDIR /home/jovyan/work

# set default command to launch when container is run
CMD ["jupyter", "lab", "--ip='0.0.0.0'", "--port=8888", "--no-browser", "--NotebookApp.token=''", "--NotebookApp.password=''"]
