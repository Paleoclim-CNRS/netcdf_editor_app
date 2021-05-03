FROM ubuntu

# Needed for mpich so doesnt prompt questions
ENV DEBIAN_FRONTEND noninteractive
RUN apt update && apt install -y fcm libnetcdf-dev libnetcdff-dev gcc gfortran make ksh nco netcdf-bin rsync


WORKDIR /usr/src
# COPY MOZAIC.tar.gz .
# RUN tar -xvzf MOZAIC.tar.gz MOZAIC && rm MOZAIC.tar.gz

# RUN svn co http://forge.ipsl.jussieu.fr/igcmg/svn/TOOLS/MOZAIC MOZAIC && \
RUN bash -l -c 'echo export NETCDF_LIB_DIR=\"$(nc-config --libs)\" >> /usr/src/setup_env.sh' && \
    bash -l -c 'echo export NETCDF_INC_DIR=\"$(nc-config --cflags)\" >> /usr/src/setup_env.sh' && \
    bash -l -c 'echo export NETCDFC_LDFLAGS=\"$(nc-config --libs)\" >> /usr/src/setup_env.sh' && \
    bash -l -c 'echo export NETCDFFORTRAN_LDFLAGS=\"$(nc-config --flibs)\" >> /usr/src/setup_env.sh' && \
    bash -l -c 'echo export NETCDFC_CFLAGS=\"$(nc-config --cflags)\" >> /usr/src/setup_env.sh' && \
    bash -l -c 'echo export NETCDFFORTRAN_FFLAGS=\"$(nc-config --fflags)\" >> /usr/src/setup_env.sh' && \
    bash -l -c 'echo source /usr/src/setup_env.sh >> /etc/bash.bashrc'

SHELL ["/bin/bash", "-c"]

# Copy the build config that should work with gfortran
# COPY bld.cfg /usr/src/MOZAIC
# COPY MOZAIC_FROM_IRENE /usr/src/MOZAIC
COPY MOZAIC_DOCKER.tar.gz /usr/src
RUN tar -xvzf MOZAIC_DOCKER.tar.gz MOZAIC && rm MOZAIC_DOCKER.tar.gz
WORKDIR /usr/src/MOZAIC
RUN source /usr/src/setup_env.sh && fcm build -f

WORKDIR /usr/src
# Modified version of DOMSK to compile with gfortran
# SRC files modified for variable character lenghts
# SRC files modified for implicit logical * Real calculations
# COPY DOMSK.tar.gz .
COPY DOMSK_DOCKER.tar.gz /usr/src
# COPY DOMSK ./DOMSK
RUN tar -xvzf DOMSK_DOCKER.tar.gz DOMSK && rm DOMSK_DOCKER.tar.gz
# The Logname is used in different places and written as an attribute to file 
# See calls to flioputa with LOGNAME
# Dont know if this is important
# ENV LOGNAME DOCKER
WORKDIR /usr/src/DOMSK
RUN source /usr/src/setup_env.sh && fcm build -f

WORKDIR /usr/src
COPY util_docker.tar.gz /usr/src
RUN tar -xvzf util_docker.tar.gz util && rm util_docker.tar.gz

RUN apt install -y 

WORKDIR /tmp
# COPY orca2.3.nc /tmp/data/orca2.3.nc
COPY PALEORCA2.30MaTotV1/*.nc  /tmp/data/
COPY mosaic.py /tmp