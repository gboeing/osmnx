########################################################################
# OSMnx Dockerfile
# License: MIT, see full license in LICENSE.txt
# Web: https://github.com/gboeing/osmnx
#
# Build an image from the dockerfile:
# >>> docker build -t gboeing/osmnx .
#
# Push the built image to hub so others can pull/run it:
# >>> docker tag gboeing/osmnx gboeing/osmnx:v0.0.0
# >>> docker login
# >>> docker push gboeing/osmnx
#
# Run bash in this container and export final conda environment to a yml file:
# >>> docker run --rm -it -v %cd%:/home/jovyan/work gboeing/osmnx /bin/bash
# >>> conda env export -n base > /home/jovyan/work/environment.yml
#
# Run jupyter lab in this container:
# >>> docker run --rm -it -p 8888:8888 -v %cd%:/home/jovyan/work gboeing/osmnx
#
# Stop/delete all local docker containers/images:
# >>> docker stop $(docker ps -aq)
# >>> docker rm $(docker ps -aq)
# >>> docker rmi $(docker images -q)
########################################################################

FROM continuumio/miniconda3
LABEL maintainer="Geoff Boeing <g.boeing@northeastern.edu>"
LABEL url="https://github.com/gboeing/osmnx"
LABEL description="OSMnx is a Python package to retrieve, model, analyze, and visualize OpenStreetMap networks and other spatial data."

# configure conda and install packages in one RUN to keep image tidy
RUN conda config --set show_channel_urls true && \
    conda config --prepend channels conda-forge && \
    conda update --strict-channel-priority --yes -n base conda && \
    conda install --strict-channel-priority --update-all --force-reinstall --yes jupyterlab osmnx python-igraph && \
    conda clean --yes --all && \
    conda info --all && \
    conda list

# launch jupyter in the local working directory that we mount
WORKDIR /home/jovyan/work

# set default command to launch when container is run
CMD ["jupyter", "lab", "--ip='0.0.0.0'", "--port=8888", "--no-browser", "--allow-root", "--NotebookApp.token=''", "--NotebookApp.password=''"]

# to test, import OSMnx and print its version
RUN ipython -c "import osmnx; print(osmnx.__version__)"
