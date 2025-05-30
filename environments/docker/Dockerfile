FROM jupyter/base-notebook
LABEL maintainer="Geoff Boeing <boeing@usc.edu>"
LABEL url="https://osmnx.readthedocs.io"
LABEL description="OSMnx is a Python package to easily download, model, analyze, and visualize street networks and other geospatial features from OpenStreetMap."

# expose $TARGETPLATFORM to the install.sh script
ARG TARGETPLATFORM

# copy the package files needed for installation
COPY --chown=jovyan --chmod=0755 ./environments/docker/install.sh ./osmnx/
COPY --chown=jovyan --chmod=0755 ./environments/requirements-extras.txt ./osmnx/
COPY --chown=jovyan --chmod=0755 ./osmnx/ ./osmnx/osmnx/
COPY --chown=jovyan --chmod=0755 ./pyproject.toml ./osmnx/
COPY --chown=jovyan --chmod=0755 ./README.md ./osmnx/

# install and configure everything in one RUN to keep image tidy
RUN cd ./osmnx && bash install.sh

# set jupyter working directory to map to mounted volume
WORKDIR /home/jovyan/work

# set default command to launch when container is run
CMD ["jupyter", "lab", "--ip='0.0.0.0'", "--port=8888", "--no-browser", "--NotebookApp.token=''", "--NotebookApp.password=''"]
