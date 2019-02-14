########################################################################
# OSMnx Dockerfile
# License: MIT, see full license in LICENSE.txt
# Web: https://github.com/gboeing/osmnx
#
# Build an image from the dockerfile:
# >>> docker build -t gboeing/osmnx .
#
# Push the built image to hub so others can pull/run it:
# >>> docker tag gboeing/osmnx gboeing/osmnx
# >>> docker login
# >>> docker push gboeing/osmnx
#
# Run bash in this container as root (via -u 0):
# >>> docker run --rm -it -u 0 --name osmnx -v %cd%:/home/jovyan/work gboeing/osmnx /bin/bash
#
# Run jupyter lab in this container:
# >>> docker run --rm -it --name osmnx -p 8888:8888 -v %cd%:/home/jovyan/work gboeing/osmnx
#
# Stop/delete all local docker containers/images:
# >>> docker stop $(docker ps -aq)
# >>> docker rm $(docker ps -aq)
# >>> docker rmi $(docker images -q)
########################################################################

FROM jupyter/base-notebook
LABEL maintainer="Geoff Boeing <g.boeing@northeastern.edu>"

# symlink and permissions
USER root
RUN ln -s /opt/conda/bin/jupyter /usr/local/bin
USER $NB_UID

# configure conda and install packages in one RUN to keep image tidy
RUN conda config --set show_channel_urls true && \
    conda config --prepend channels conda-forge && \
    conda update --override-channels -c conda-forge --yes -n base conda && \
    conda update --override-channels -c conda-forge --yes --all && \
    conda install --override-channels -c conda-forge --yes osmnx python-igraph && \
    conda clean --yes --all && \
    conda info --all && \
    conda list

# launch notebook in the local working directory that we mount
WORKDIR /home/jovyan/work

# set default command to launch when container is run
CMD ["jupyter", "lab", "--no-browser", "--NotebookApp.token=''", "--NotebookApp.password=''"]
