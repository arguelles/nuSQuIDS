# build with: `docker build -t nusquids .`

# Use the official Ubuntu base image
FROM ubuntu:20.04

# Set non-interactive frontend for automatic installation
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    build-essential \
    g++ \
    gcc \
    wget \
    python3 \
    python3-pip \
    python3-dev \
    git \
    openssh-client \
    pkg-config \
    libgsl-dev \
    libboost-all-dev \
    libboost-python-dev \
    libboost-python1.71.0 \
    libhdf5-dev \
    hdf5-tools \
    cmake \
    python3-sphinx \
    && rm -rf /var/lib/apt/lists/* \
    && ln -s /usr/bin/python3 /usr/bin/python \
#    && ln -s /usr/lib/x86_64-linux-gnu/libpython3.8.so /usr/lib/x86_64-linux-gnu/libpython3.so \
    && pip install numpy h5py

ENV CXX=g++
ENV CC=gcc
ENV SPACE=/nusquids-space
ENV BUILDPATH=$SPACE/local
ENV SOURCEPATH=$SPACE/sources
ENV PREFIX=$BUILDPATH
ENV PATH=$BUILDPATH/bin:$PATH
ENV LD_LIBRARY_PATH=$BUILDPATH/lib/:$LD_LIBRARY_PATH
ENV C_INCLUDE_PATH=$BUILDPATH/include/:$C_INCLUDE_PATH
ENV CPLUS_INCLUDE_PATH=/usr/include/hdf5/serial/:$BUILDPATH/include/:$CPLUS_INCLUDE_PATH
ENV CXX_INCLUDE_PATH=$BUILDPATH/include/:$CXX_INCLUDE_PATH
ENV PKG_CONFIG_PATH=$BUILDPATH/lib/pkgconfig:$PKG_CONFIG_PATH
ENV PYTHONPATH=$BUILDPATH/lib/python3.8/site-packages:$PYTHONPATH
ENV HDF5_DISABLE_VERSION_CHECK=1

RUN mkdir -p $SPACE $BUILDPATH $SOURCEPATH && \
mkdir -p $PREFIX/bin $PREFIX/include $PREFIX/lib $PREFIX/lib64

WORKDIR $SOURCEPATH

# The following can be uncommented to compile boost from source. That did not help so far with our problem. The env variable CPLUS_INCLUDE_PATH has to be set after the ./bootstrap.sh command, otherwise it fails for some reason

# # Download and extract Boost 1.71.0
# RUN wget https://boostorg.jfrog.io/artifactory/main/release/1.71.0/source/boost_1_71_0.tar.bz2 \
# && tar --bzip2 -xf boost_1_71_0.tar.bz2 \
# && rm boost_1_71_0.tar.bz2

# # Set working directory
# WORKDIR $SOURCEPATH/boost_1_71_0
# RUN ./bootstrap.sh --prefix=$PREFIX --with-python=/usr/bin/python3 && ./b2 install

# ENV CPLUS_INCLUDE_PATH=$GOLEMBUILDPATH/include/:$CPLUS_INCLUDE_PATH

# WORKDIR $SOURCEPATH

RUN git clone https://github.com/jsalvado/SQuIDS.git SQuIDS && \
    git clone https://github.com/arguelles/nuSQuIDS.git nuSQuIDS

WORKDIR $SOURCEPATH/SQuIDS
RUN ./configure --prefix=/usr && make && make install

WORKDIR $SOURCEPATH/nuSQuIDS
RUN ./configure --prefix=/usr --with-python-bindings \
    && make \
    && make python \
    && make install

WORKDIR $SOURCEPATH
