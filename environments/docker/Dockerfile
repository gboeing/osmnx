########################################################################
# OSMnx Dockerfile
# License: MIT, see full license in LICENSE.txt
# Web: https://github.com/gboeing/osmnx
########################################################################

FROM continuumio/miniconda3
LABEL maintainer="Geoff Boeing <boeing@usc.edu>"
LABEL url="https://github.com/gboeing/osmnx"
LABEL description="OSMnx is a Python package to retrieve, model, analyze, and visualize OpenStreetMap networks and other spatial data."

COPY requirements.txt /tmp/

# configure conda and install packages in one RUN to keep image tidy
RUN conda config --set show_channel_urls true && \
    conda config --set channel_priority strict && \
    conda config --prepend channels conda-forge && \
    conda update --yes -n base conda && \
    conda install --update-all --force-reinstall --yes --file /tmp/requirements.txt python-igraph && \
    python -m ipykernel install --name ox --display-name "Python (ox)" && \
    conda clean --all --yes && \
    conda info --all && \
    conda list

# launch jupyter in the local working directory that we mount
WORKDIR /home/jovyan/work

# set default command to launch when container is run
CMD ["jupyter", "lab", "--ip='0.0.0.0'", "--port=8888", "--no-browser", "--allow-root", "--NotebookApp.token=''", "--NotebookApp.password=''"]

# to test, import OSMnx and print its version
RUN ipython -c "import osmnx; print(osmnx.__version__)"
