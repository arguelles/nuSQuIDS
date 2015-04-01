#!/bin/bash

function check_pkgconfig(){
	if [ "$CHECKED_PKGCONFIG" ]; then return; fi
	echo "Looking for pkg-config..."
	which pkg-config 2>&1 > /dev/null
	if [ "$?" -ne 0 ]; then
		echo "Error: pkg-config not found; you will need to specify library locations manually" 1>&2
		exit 1
	fi
	CHECKED_PKGCONFIG=1
}

function find_package(){
	PKG=$1
	VAR_PREFIX=`echo $PKG | tr [:lower:] [:upper:]`
	TMP_FOUND=`eval echo "$"${VAR_PREFIX}_FOUND`
	if [ "$TMP_FOUND" ]; then return; fi
	check_pkgconfig
	echo "Looking for $PKG..."

	pkg-config --exists $PKG
	if [ "$?" -ne 0 ]; then
		echo "Error: $PKG not installed or not registered with pkg-config" 1>&2
    lowername=`echo $PKG | tr [A-Z] [a-z]`
    echo "Please specify location using the --with-"$lowername"-incdir and --with-"$lowername"-libdir flags" 1>&2
		exit 1
	fi
	if [ $# -ge 2 ]; then
		MIN_VERSION=$2
		pkg-config --atleast-version $MIN_VERSION $PKG
		if [ "$?" -ne 0 ]; then
			echo "Error: installed $PKG verson ("`pkg-config --modversion $PKG`") is too old; version >=$MIN_VERSION is required" 1>&2
			exit 1
		fi
	fi
	eval ${VAR_PREFIX}_FOUND=1
	eval ${VAR_PREFIX}_CFLAGS=\"`pkg-config --cflags $PKG`\"
	eval ${VAR_PREFIX}_LDFLAGS=\"`pkg-config --libs $PKG`\"
}

PREFIX=/usr/local

VERSION_NUM=100000
VERSION=`echo $VERSION_NUM | awk '{
	major = int($1/100000);
	minor = ($1/100)%1000;
	patch = $1%100;
	print major"."minor"."patch;
}'`

OS_NAME=`uname -s`

GUESS_CC=gcc
GUESS_CXX=g++
GUESS_AR=ar
GUESS_LD=ld
if [ "$OS_NAME" = Linux ]; then
	DYN_SUFFIX=.so
	DYN_OPT='-shared -Wl,-soname,$(shell basename $(DYN_PRODUCT))'
fi
if [ "$OS_NAME" = Darwin ]; then
	GUESS_CC=clang
	GUESS_CXX=clang++
	GUESS_LD=clang++
	DYN_SUFFIX=.dylib
  DYN_OPT='-dynamiclib -install_name $(LIBnuSQUIDS)/$(DYN_PRODUCT) -compatibility_version $(VERSION) -current_version $(VERSION)'
fi

CC=${CC-$GUESS_CC}
CXX=${CXX-$GUESS_CXX}
AR=${AR-$GUESS_AR}
LD=${LD-$GUESS_LD}

HELP="Usage: ./config.sh [OPTION]... 

Installation directories:
  --prefix=PREFIX         install files in PREFIX
                          [$PREFIX]

By default, \`make install' will install all the files in
\`$PREFIX/bin', \`$PREFIX/lib' etc.  You can specify
an installation prefix other than \`$PREFIX' using \`--prefix',
for instance \`--prefix=\$HOME'.

The following options can be used to maunally specify the 
locations of dependencies:
  --with-gsl-incdir=DIR    use the copy of gsl in DIR
  --with-gsl-libdir=DIR    use the copy of gsl in DIR
  --with-hdf5-incdir=DIR   use the copy of HDF5 in DIR
  --with-hdf5-libdir=DIR   use the copy of HDF5 in DIR
  --with-squids-incdir=DIR        use the copy of SQuIDS in DIR
  --with-squids-libdir=DIR        use the copy of SQuIDS in DIR

Some influential environment variables:
CC          C compiler command
CXX         C++ compiler command
AR          Static linker command
LD          Dynamic linker command
" #`

for var in "$@"
do
	if [ "$var" = "--help" -o "$var" = "-h" ]; then
		echo "$HELP"
		exit 0
	fi

	TMP=`echo "$var" | sed -n 's/^--prefix=\(.*\)$/\1/p'`
	if [ "$TMP" ]; then PREFIX="$TMP"; continue; fi

	TMP=`echo "$var" | sed -n 's/^--with-gsl-incdir=\(.*\)$/\1/p'`
	if [ "$TMP" ]; then GSL_INCDIR="$TMP"; continue; fi

	TMP=`echo "$var" | sed -n 's/^--with-gsl-libdir=\(.*\)$/\1/p'`
	if [ "$TMP" ]; then GSL_LIBDIR="$TMP"; continue; fi

	TMP=`echo "$var" | sed -n 's/^--with-hdf5-incdir=\(.*\)$/\1/p'`
	if [ "$TMP" ]; then HDF5_INCDIR="$TMP"; continue; fi

	TMP=`echo "$var" | sed -n 's/^--with-hdf5-libdir=\(.*\)$/\1/p'`
	if [ "$TMP" ]; then HDF5_LIBDIR="$TMP"; continue; fi

	TMP=`echo "$var" | sed -n 's/^--with-squids-incdir=\(.*\)$/\1/p'`
	if [ "$TMP" ]; then SQUIDS_INCDIR="$TMP"; continue; fi

	TMP=`echo "$var" | sed -n 's/^--with-squids-libdir=\(.*\)$/\1/p'`
	if [ "$TMP" ]; then SQUIDS_LIBDIR="$TMP"; continue; fi

	TMP=`echo "$var" | sed -n 's/^--with-python/true/p'`
	if [ "$TMP" ]; then PYTHON_BINDINGS=true; continue; fi
done

if [ "$GSL_INCDIR" -a "$GSL_LIBDIR" ]; then
	echo "Checking manually specified GSL..."
	if [ -d "$GSL_INCDIR/gsl" \
         -a -f "$GSL_INCDIR/gsl/gsl_version.h" \
         -a -d "$GSL_LIBDIR" \
         -a -f "$GSL_LIBDIR/libgsl.a" ]; then
		GSL_FOUND=1
		GSL_CFLAGS="-I$GSL_INCDIR"
		GSL_LDFLAGS="-L$GSL_LIBDIR -lgsl -lgslcblas -lm"
	else
		echo "Warning: manually specifed GSL not found; will attempt auto detection"
	fi
fi

find_package gsl 1.15

if [ "$HDF5_INCDIR" -a "$HDF5_LIBDIR" ]; then
	echo "Checking manually specified HDF5..."
	if [ -d "$HDF5_INCDIR" \
         -a -f "$HDF5_INCDIR/H5version.h" \
         -a -d "$HDF5_LIBDIR" \
         -a -f "$HDF5_LIBDIR/libhdf5.a" \
         -a -f "$HDF5_LIBDIR/libhdf5_hl.a" ]; then
		HDF5_FOUND=1
		HDF5_CFLAGS="-I$HDF5_INCDIR"
		HDF5_LDFLAGS="-L$HDF5_LIBDIR -lhdf5 -lhdf5_hl"
	else
		echo "Warning: manually specifed HDF5 not found; will attempt auto detection"
	fi
fi

find_package HDF5 1.8

if [ "$SQUIDS_INCDIR" -a "$SQUIDS_LIBDIR" ]; then
	echo "Checking manually specified SQUIDS..."
	if [ -d "$SQUIDS_INCDIR" \
         -a -d "$SQUIDS_LIBDIR" \
         -a -f "$SQUIDS_LIBDIR/libSQuIDS.a" ]; then
		SQUIDS_FOUND=1
		SQUIDS_CFLAGS="-I$SQUIDS_INCDIR"
		SQUIDS_LDFLAGS="-L$SQUIDS_LIBDIR -lSQuIDS"
	else
		echo "Warning: manually specifed SQUIDS not found; will attempt auto detection"
	fi
fi

find_package squids 1.2

if [ ! -d ./lib/ ]; then
    mkdir lib;
fi


echo "Generating config file..."
echo "prefix=$PREFIX" > nusquids.pc
echo '
libdir=${prefix}/lib
includedir=${prefix}/inc

Name: nuSQuIDS
Description: Toolbox for neutrino oscillation experiments
URL: https://github.com/arguelles/nuSQuIDS' >> nusquids.pc
echo "Version: $VERSION" >> nusquids.pc
echo 'Requires: gsl >= 1.15
Requires: hdf5 >= 1.8
Requires: squids >= 1.2.0
Libs: -L${libdir} -lnuSQuIDS
Cflags: -I${includedir}
' >> nusquids.pc

echo "Generating version header..."
sed -e "s|__NUSQUIDS_VERSION__|$VERSION_NUM|g" \
    -e "s|__NUSQUIDS_VERSION_STR__|$VERSION|g" \
    < resources/version.h.in > inc/version.h

echo "Generating makefile..."
echo "# Compiler
CC=$CC
CXX=$CXX
AR=$AR
LD=$LD

DYN_SUFFIX=$DYN_SUFFIX
DYN_OPT=$DYN_OPT

VERSION=$VERSION
PREFIX=$PREFIX
" > ./Makefile

echo '
PATH_nuSQUIDS=$(shell pwd)
PATH_SQUIDS=$(SQUIDS_DIR)

SOURCES = $(wildcard src/*.cpp)
OBJECTS = $(SOURCES:.cpp=.o)

EXAMPLES_SRC=$(wildcard examples/*.cpp)
EXAMPLES=$(EXAMPLES_SRC:.cpp=.exe)

CXXFLAGS= -std=c++11

# Directories
'  >> ./Makefile
echo "GSL_CFLAGS=$GSL_CFLAGS" >> ./Makefile
echo "GSL_LDFLAGS=$GSL_LDFLAGS" >> ./Makefile

echo "HDF5_CFLAGS=$HDF5_CFLAGS" >> ./Makefile
echo "HDF5_LDFLAGS=$HDF5_LDFLAGS" >> ./Makefile

echo "SQUIDS_CFLAGS=$SQUIDS_CFLAGS" >> ./Makefile
echo "SQUIDS_LDFLAGS=$SQUIDS_LDFLAGS" >> ./Makefile
echo '

INCnuSQUIDS=$(PATH_nuSQUIDS)/inc
LIBnuSQUIDS=$(PATH_nuSQUIDS)/lib

# FLAGS
CFLAGS= -O3 -fPIC -I$(INCnuSQUIDS) $(SQUIDS_CFLAGS) $(GSL_CFLAGS) $(HDF5_CFLAGS)
LDFLAGS= -Wl,-rpath -Wl,$(LIBnuSQUIDS) -L$(LIBnuSQUIDS) -lnuSQuIDS
LDFLAGS+= $(SQUIDS_LDFLAGS) $(GSL_LDFLAGS) $(HDF5_LDFLAGS)

# Project files
NAME=nuSQuIDS
STAT_PRODUCT=lib$(NAME).a
DYN_PRODUCT=lib$(NAME)$(DYN_SUFFIX)

# Compilation rules
all: $(STAT_PRODUCT) $(DYN_PRODUCT)

examples : $(EXAMPLES)

%.exe : %.cpp
	$(CXX) $(CXXFLAGS) $(CFLAGS) $< $(LDFLAGS) -o $@
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
	$(CXX) $(CXXFLAGS) -c $(CFLAGS) $< -o $@

.PHONY: clean
clean:
	rm -f src/*.o examples/*.exe lib/* bin/*

doxygen:
	@doxygen src/doxyfile
docs:
	@doxygen src/doxyfile

test: $(DYN_PRODUCT) $(STAT_PRODUCT)
	@cd ./test ; ./run_tests

install: $(DYN_PRODUCT) $(STAT_PRODUCT)
	@echo Installing headers in $(PREFIX)/include/nuSQuIDS
	@mkdir -p $(PREFIX)/include/nuSQuIDS
	@cp $(INCnuSQUIDS)/*.h $(PREFIX)/include/nuSQuIDS
	@echo Installing libraries in $(PREFIX)/lib
	@mkdir -p $(PREFIX)/lib
	@cp $(LIBnuSQUIDS)/$(DYN_PRODUCT) $(LIBnuSQUIDS)/$(STAT_PRODUCT) $(PREFIX)/lib
	@echo Installing config information in $(PREFIX)/lib/pkgconfig
	@mkdir -p $(PREFIX)/lib/pkgconfig
	@cp nusquids.pc $(PREFIX)/lib/pkgconfig
' >> ./Makefile

if [ $PYTHON_BINDINGS ]; then
  echo "Generating Python bindings makefile..."
  PYTHONVERSION=`python -c 'import sys; print str(sys.version_info.major)+"."+str(sys.version_info.minor)'`
  PYTHONLIBPATH=`python -c 'import sys; import re; print [ y for y in sys.path if re.search("\/lib\/python'$PYTHONVERSION'$",y)!=None ][0];'`
  PYTHONINCPATH=$PYTHONLIBPATH/../../include/python$PYTHONVERSION
  echo "
# Compiler
CC=$CC
CXX=$CXX
AR=$AR
LD=$LD

DYN_SUFFIX=.so
DYN_OPT=$DYN_OPT

PYTHON_VERSION = ${PYTHONVERSION}
LDFLAGS+= -L${PYTHONLIBPATH}
LDFLAGS+= -lpython${PYTHON_VERSION} -lboost_python
LDFLAGS+= ${SQUIDS_LDFLAGS} ${GSL_LDFLAGS} ${HDF5_LDFLAGS}
INCCFLAGS+= -I${PYTHONINCPATH}  -I../inc/
" > resources/python/src/Makefile

echo '
PATH_nuSQUIDS= $(shell pwd)/../../../
PATH_nuSQUIDSpy= $(shell pwd)/../

SOURCES = $(wildcard *.cpp)
OBJECTS = $(SOURCES:.cpp=.o)

# Directories
LIBDIR=$(PATH_nuSQUIDS)/lib
INCDIR=$(PATH_nuSQUIDS)/inc
LIBnuSQUIDS=$(PATH_nuSQUIDS)/lib
INCnuSQUIDS=$(PATH_nuSQUIDS)/inc

LDFLAGS+= -L$(LIBnuSQUIDS) -lnuSQuIDS
INCCFLAGS+= -I$(INCnuSQUIDS)
CXXFLAGS= -O3 -fPIC -std=c++11 $(INCCFLAGS)
' >> resources/python/src/Makefile

echo '

# Project files
NAME=nuSQUIDSpy
STAT_PRODUCT=$(NAME).a
DYN_PRODUCT=$(NAME)$(DYN_SUFFIX)

OS_NAME=$(shell uname -s)
ifeq ($(OS_NAME),Linux)
	DYN_SUFFIX=.so
	DYN_OPT=-shared -Wl,-soname,$(DYN_PRODUCT)
endif
ifeq ($(OS_NAME),Darwin)
	DYN_SUFFIX=.so
	DYN_OPT=-dynamiclib -install_name $(PATH_nuSQUIDSpy)/bindings/$(DYN_PRODUCT)
endif

# Compilation rules
all: $(STAT_PRODUCT) $(DYN_PRODUCT)

$(DYN_PRODUCT) : $(OBJECTS)
	@echo Linking $(DYN_PRODUCT)
	@$(CXX) $(DYN_OPT)  $(LDFLAGS) -o $(DYN_PRODUCT) $(OBJECTS)
	mv $(DYN_PRODUCT) $(PATH_nuSQUIDSpy)/bindings/$(DYN_PRODUCT)

$(STAT_PRODUCT) : $(OBJECTS)
	@echo Linking $(STAT_PRODUCT)
	@$(AR) -rcs $(STAT_PRODUCT) $(OBJECTS)
	mv $(STAT_PRODUCT) $(PATH_nuSQUIDSpy)/bindings/$(STAT_PRODUCT)

%.o : %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

.PHONY: clean
clean:
	rm -f *.o ../lib/*.so ../lib/*.a
' >> resources/python/src/Makefile
fi

nusqpath=`pwd`

echo '
#ifndef __GLOBAL_H
#define __GLOBAL_H
' > inc/global.h
echo "
#define XSECTION_LOCATION \"$nusqpath/data/xsections/\"
#define SUN_MODEL_LOCATION  \"$nusqpath/data/astro/bs05_agsop.dat\"
#define SUN_MODEL_NELECTRON_LOCATION \"$nusqpath/data/astro/nele_bs05op.dat\"
#define EARTH_MODEL_LOCATION \"$nusqpath/data/astro/EARTH_MODEL_PREM.dat\"
" >> inc/global.h
echo '
#endif
' >> inc/global.h

echo "Done."
echo "To build library, run the following: make
After, to build examples: make examples"
