#!/bin/bash

for var in "$@"
do
    GSLINC+=`expr "$var" : '\(-GSL_INC=.*\)'`
    GSLLIB+=`expr "$var" : '\(-GSL_LIB=.*\)'`
    SQPATH+=`expr "$var" : '\(-SQUIDS_PATH=.*\)'`
    HDF5LIB+=`expr "$var" : '\(-HDF5_LIB=.*\)'`
    HDF5INC+=`expr "$var" : '\(-HDF5_INC=.*\)'`
done

GSLINC_PATH=${GSLINC##-GSL_INC=}
GSLLIB_PATH=${GSLLIB##-GSL_LIB=}
HDF5LIB_PATH=${HDF5INC##-HDF5_LIB=}
HDF5INC_PATH=${HDF5LIB##-HDF5_INC=}
SQUIDS_PATH=${SQPATH##-SQUIDS_PATH=}

echo '
# Compiler
CC=gcc
CXX=g++
FC=gfortran
AR=ar
LD=ld

PATH_nuSQUIDS= $(shell pwd)
PATH_SQUIDS=$(PATH_nuSQUIDS)/../SQuIDS

SOURCES = $(wildcard src/*.cpp)
OBJECTS = $(SOURCES:.cpp=.o)

EXAMPLES_SRC=$(wildcard examples/*.cpp)
EXAMPLES=$(EXAMPLES_SRC:.cpp=.exe)

'  > ./Makefile

if [ -n "$GSLLIB" ]; then
echo "LIBDIR=-L${GSLLIB_PATH}" >> ./Makefile
echo "INCDIR=-I${GSLINC_PATH}" >> ./Makefile
fi

if [ -n "$HDF5LIB" ]; then
echo "LIBDIR+=-L${HDF5LIB_PATH}" >> ./Makefile
echo "INCDIR+=-I${HDF5INC_PATH}" >> ./Makefile
fi

# Directories
if [ -n "$SQUIDS_PATH" ]; then
echo "PATH_SQUIDS=$(SQUIDS_PATH)" >> ./Makefile
fi

echo '
LIBSQUIDS=$(PATH_SQUIDS)/lib
INCSQUIDS=$(PATH_SQUIDS)/inc

INCnuSQUIDS=$(PATH_nuSQUIDS)/inc
LIBnuSQUIDS=$(PATH_nuSQUIDS)/lib

# FLAGS
LDFLAGS+= $(LIBDIR) -L$(LIBSQUIDS)
LDFLAGS+= -lgsl -lgslcblas
LDFLAGS+= -lSQUIDS
LDFLAGS+= -lhdf5 -lhdf5_hl -lhdf5_hl_cpp
INCCFLAGS = $(INCDIR) -I$(INCSQUIDS) -I$(INCnuSQUIDS)
CXXFLAGS= -O3 -fPIC -std=c++11 $(INCCFLAGS)

# Project files
NAME=nuSQUIDS
STAT_PRODUCT=lib$(NAME).a
DYN_PRODUCT=lib$(NAME)$(DYN_SUFFIX)

OS_NAME=$(shell uname -s)
ifeq ($(OS_NAME),Linux)
	DYN_SUFFIX=.so
	DYN_OPT=-shared -Wl,-soname,$(DYN_PRODUCT)
endif
ifeq ($(OS_NAME),Darwin)
  CC=clang
	CXX=clang++
	LD=clang++
	DYN_SUFFIX=.dylib
	DYN_OPT=-dynamiclib -install_name $(LIBnuSQUIDS)/$(DYN_PRODUCT)
endif


# Compilation rules
all: $(STAT_PRODUCT) $(DYN_PRODUCT) $(EXAMPLES)

%.exe : %.cpp
	$(CXX) $(CXXFLAGS) $(LIBDIR) -L$(LIBSQUIDS) -L$(LIBnuSQUIDS) $< -o $@ -lSQUIDS -lnuSQUIDS
	mv $@ bin/

$(DYN_PRODUCT) : $(OBJECTS)
	@echo Linking $(DYN_PRODUCT)
	@$(CXX) $(DYN_OPT)  $(LDFLAGS) -o $(DYN_PRODUCT) $(OBJECTS)
	mv $(DYN_PRODUCT) $(PATH_nuSQUIDS)/lib/$(DYN_PRODUCT)

$(STAT_PRODUCT) : $(OBJECTS)
	@echo Linking $(STAT_PRODUCT)
	@$(AR) -rcs $(STAT_PRODUCT) $(OBJECTS)
	mv $(STAT_PRODUCT) $(PATH_nuSQUIDS)/lib/$(STAT_PRODUCT)

%.o : %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

.PHONY: clean
clean:
	rm -f src/*.o examples/*.exe lib/* bin/*
' >> ./Makefile

nusqpath=`pwd`

echo '
#ifndef __GLOBAL_H
#define __GLOBAL_H
' > inc/global.h
echo "
#define XSECTION_LOCATION '$nusqpath/data/xsections/'
#define SUN_MODEL_LOCATION  '$nusqpath/data/astro/bs05_agsop.dat'
#define SUN_MODEL_NELECTRON_LOCATION '$nusqpath/data/astro/nele_bs05op.dat'
" >> inc/global.h
echo '
#endif
' >> inc/global.h
