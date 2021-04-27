FROM ubuntu

# Needed for mpich so doesnt prompt questions
ENV DEBIAN_FRONTEND noninteractive
RUN apt update && apt install -y fcm
RUN apt update && apt install -y vim
RUN apt update && apt install -y libnetcdf-dev libnetcdff-dev
RUN apt update && apt install -y gcc gfortran make ksh


WORKDIR /usr/src
# COPY MOZAIC.tar.gz .
# RUN tar -xvzf MOZAIC.tar.gz MOZAIC && rm MOZAIC.tar.gz

RUN svn co http://forge.ipsl.jussieu.fr/igcmg/svn/TOOLS/MOZAIC MOZAIC && \
    bash -l -c 'echo export NETCDF_LIB_DIR=\"$(nc-config --libs)\" >> /usr/src/setup_env.sh' && \
    bash -l -c 'echo export NETCDF_INC_DIR=\"$(nc-config --cflags)\" >> /usr/src/setup_env.sh' && \
    bash -l -c 'echo export NETCDFC_LDFLAGS=\"$(nc-config --libs)\" >> /usr/src/setup_env.sh' && \
    bash -l -c 'echo export NETCDFFORTRAN_LDFLAGS=\"$(nc-config --flibs)\" >> /usr/src/setup_env.sh' && \
    bash -l -c 'echo export NETCDFC_CFLAGS=\"$(nc-config --cflags)\" >> /usr/src/setup_env.sh' && \
    bash -l -c 'echo export NETCDFFORTRAN_FFLAGS=\"$(nc-config --fflags)\" >> /usr/src/setup_env.sh' && \
    bash -l -c 'echo source /usr/src/setup_env.sh >> /etc/bash.bashrc'

SHELL ["/bin/bash", "-c"]

# Copy the build config that should work with gfortran
COPY bld.cfg /usr/src/MOZAIC
WORKDIR /usr/src/MOZAIC
RUN source /usr/src/setup_env.sh && fcm build -f

WORKDIR /usr/src
# Modified version of DOMSK to compile with gfortran
# SRC files modified for variable character lenghts
# SRC files modified for implicit logical * Real calculations
COPY DOMSK.tar.gz .
RUN tar -xvzf DOMSK.tar.gz DOMSK && rm DOMSK.tar.gz
WORKDIR /usr/src/DOMSK
RUN source /usr/src/setup_env.sh && fcm build -f

