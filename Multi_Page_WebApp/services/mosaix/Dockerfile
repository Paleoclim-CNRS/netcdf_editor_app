FROM ubuntu

# Needed for mpich so doesnt prompt questions
ENV DEBIAN_FRONTEND noninteractive \
    H5DIR /usr/local/hdf5 \
    NCDIR /usr/local/netcdf \
    NFDIR /usr/local/netcdff

RUN apt update && \
    apt install -y perl make m4 mpich subversion wget curl zlib1g-dev libcurl4-gnutls-dev \
    nco python3 python3-pip uuid && python3 -m pip install netcdf4 scipy\
    && ln -s /usr/bin/make /usr/bin/gmake

RUN apt update && apt install -y git ksh \
    && cpan install URI

# Download source files
RUN mkdir /src && cd /src && \
    wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.12/hdf5-1.12.0/src/hdf5-1.12.0.tar.gz \
    && cd /src \
    && tar -xvzf hdf5-1.12.0.tar.gz \
    && rm hdf5-1.12.0.tar.gz \
    && cd /src/hdf5-1.12.0 \
    && CC=/usr/bin/mpicc  ./configure --prefix=${H5DIR} --enable-parallel --with-zlib \
    && make \
    && make install

WORKDIR /src
RUN wget https://github.com/Unidata/netcdf-c/archive/v4.7.4.tar.gz \
    && tar -xvzf  v4.7.4.tar.gz \
    && rm -rf  v4.7.4.tar.gz \
    && cd netcdf-c-4.7.4 \
    && CC=/usr/bin/mpicc CPPFLAGS=-I${H5DIR}/include LDFLAGS=-L${H5DIR}/lib   ./configure --enable-parallel-tests --prefix=${NCDIR} --enable-remote-fortran-bootstrap\
    && make install \
    && make build-netcdf-fortran \
    && cd /src/netcdf-c-4.7.4/netcdf-fortran \
    && CC=/usr/bin/mpicc FC=/usr/bin/mpif90 CPPFLAGS=-I${H5DIR}/include\ -I${NCDIR}/include LDFLAGS=-L${H5DIR}/lib   ./configure --prefix=${NFDIR}\
    && make install 

RUN svn co http://forge.ipsl.jussieu.fr/ioserver/svn/XIOS/trunk XIOS
# COPY XIOS_arch_file.env XIOS/arch/arch-GCC_LINUX.env
# COPY arch-linux_ubuntu.fcm XIOS/arch/arch-GCC_LINUX.fcm
# COPY arch-linux_ubuntu.path XIOS/arch/arch-GCC_LINUX.path
COPY services/mosaix/ARCH/arch-linux_ubuntu* XIOS/arch/

RUN cd XIOS && \
     ./make_xios --arch linux_ubuntu --jobs 8

RUN svn co http://forge.ipsl.jussieu.fr/igcmg/svn/TOOLS/MOSAIX MOSAIX
RUN sed -i '/-L%{XIOS_DIR}\/lib/c\%LD_FLAGS         -L%{XIOS_DIR}/lib -lxios %BASE_LD %ARCH_LD' MOSAIX/bld.cfg
RUN apt update && apt install -y fcm
RUN cd MOSAIX && ./make_mosaix
# Needed by arch-X64_IRENE.fcm -> arch-linux_ubuntu.fcm
ENV NETCDFC_INCDIR=${NCDIR}/include \
    NETCDFC_LIBDIR=${NCDIR}/lib \
    NETCDFFORTRAN_INCDIR=${NFDIR}/include \
    NETCDFFORTRAN_LIBDIR=${NFDIR}/lib \
    LD_LIBRARY_PATH=$NETCDFFORTRAN_LIBDIR:$NCDIR/lib:$H5DIR/lib:/src/XIOS/lib\
    PATH=$NCDIR/bin:$PATH \
    HOME=/root \
    R_IN=/root/Scratch/IGCM \
    TMPDIR=/root/Scratch/TMP \
    SUBMIT_DIR=$(pwd) \
    MpiRun="/usr/bin/mpirun -n 1" \
    PyRun=""
COPY services/mosaix/PATCHS/create_weights.patch services/mosaix/PATCHS/Runoff.patch services/mosaix/PATCHS/Runoff_periodicity.patch services/mosaix/PATCHS/Mosaix.patch /src/MOSAIX/
RUN cd MOSAIX && patch < Runoff.patch \
    && patch < create_weights.patch \
    && patch < Runoff_periodicity.patch \
    && patch < Mosaix.patch

COPY services/mosaix/DATA/ORCA2.3_coordinates_mask.nc ${HOME}/Scratch/IGCM/OCE/NEMO/ORCA2.3/
COPY services/mosaix/DATA/LMD9695_grid.nc ${HOME}/Scratch/IGCM/ATM/GRID/

RUN mkdir -p /src/NEMO
WORKDIR /src/NEMO

# Download Modipsl

COPY services/mosaix/MODIFIED_MY_SRC/ /tmp/SOURCES

# RUN svn co http://forge.ipsl.jussieu.fr/igcmg/svn/modipsl/trunk modipsl \
#     && cd modipsl/util \
#     # Don't Download XIOS
#     && sed -i "/-C- NEMO_v6.1  XIOS\/branchs\/xios-2.5/d" mod.def \
#     # Download Nemo
#     && ./model NEMO_v6.1 \
#     && ./ins_make -t gfortran \
#     # Use already built XIOS
#     && ln -s /src/XIOS /src/NEMO/modipsl/modeles/XIOS \
#     # change arch from gfortran to linux_ubuntu
#     && sed -i 's/FCM_ARCH\ =\ gfortran\ /FCM_ARCH\ =\ linux_ubuntu/' /src/NEMO/modipsl/config/NEMO_v6/Makefile \
#     # Use the X64_IRENE arch-file but rename to linux_ubuntu
#     && cp /src/NEMO/modipsl/config/NEMO_v6/ARCH/arch-X64_IRENE.fcm /src/NEMO/modipsl/config/NEMO_v6/ARCH/arch-linux_ubuntu.fcm \
#     # Replace intel compile options with gfortran equivalents
#     && sed -i "s/%FCFLAGS             -i4 -r8 -O3 -fp-model precise/%FCFLAGS             -O3 -fdefault-real-8 -ffree-line-length-none/" /src/NEMO/modipsl/config/NEMO_v6/ARCH/arch-linux_ubuntu.fcm \
#     # Invert the order of the flags so that Nemo links correctly
#     && sed -i 's/%LDFLAGS %USER_LIB/%USER_LIB %LDFLAGS/' /src/NEMO/modipsl/modeles/NEMOGCM/TOOLS/COMPILE/bldxag.cfg \
#     # DO NOT (re)build XIOS
#     && sed -i 's/xios orca2lim3pisces/orca2lim3pisces/' /src/NEMO/modipsl/config/NEMO_v6/Makefile \
#     # Rename MY_SRC files to ensure copy of sources
#     && cd /src/NEMO/modipsl/modeles/NEMOGCM/CONFIG/ORCA1_LIM3_PISCES/MY_SRC \
#     && mv iom.F90 iom.F90.old \
#     && mv p4zsms.F90 p4zsms.F90.old \
#     && mv traadv.F90 traadv.F90.old \
#     # Tell makenemo to use the modified MY_SRC files
#     && sed -i "s/ORCA2_LIM3_PISCES -r ORCA1_LIM3_PISCES -j16/ORCA2_LIM3_PISCES -r ORCA1_LIM3_PISCES -e \"\$\(MYSRC_DIR\)\/SOURCES\" -j16/" /src/NEMO/modipsl/config/NEMO_v6/Makefile \
#     && mv /tmp/SOURCES/* /src/NEMO/modipsl/config/NEMO_v6/SOURCES \
#     # Build Nemo
#     && cd /src/NEMO/modipsl/config/NEMO_v6 && gmake ORCA2LIM3PISCES

# COPY DATA/OR2L3P-CLIM-CORE.85237.tar.gz /tmp

# Python messaging and running
COPY services/mosaix/receive.py services/mosaix/steps.py services/mosaix/mosaix.py /usr/src/app/
COPY flask_app /usr/src/app/flask_app
RUN pip3 install /usr/src/app/flask_app && pip3 install pika

ENV FLASK_APP=/usr/src/app/netcdf_editor_app
ENTRYPOINT [ "python3", "/usr/src/app/receive.py"]