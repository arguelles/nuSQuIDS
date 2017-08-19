#!/bin/bash

check_pkgconfig(){
	if [ "$CHECKED_PKGCONFIG" ]; then return; fi
	echo "Looking for pkg-config..."
	which pkg-config 2>&1 > /dev/null
	if [ "$?" -ne 0 ]; then
		echo "Error: pkg-config not found; you will need to specify library locations manually" 1>&2
		exit 1
	fi
	CHECKED_PKGCONFIG=1
}

find_package(){
	PKG=$1
	VAR_PREFIX=`echo $PKG | tr [:lower:] [:upper:]`
	TMP_FOUND=`eval echo "$"${VAR_PREFIX}_FOUND`
	if [ "$TMP_FOUND" ]; then return; fi
	check_pkgconfig
	echo "Looking for $PKG..."

	pkg-config --exists $PKG
	if [ "$?" -ne 0 ]; then
		echo " $PKG not found with pkg-config"
		return
	fi
	if [ $# -ge 2 ]; then
		MIN_VERSION=$2
		pkg-config --atleast-version $MIN_VERSION $PKG
		if [ "$?" -ne 0 ]; then
			echo "Error: installed $PKG version ("`pkg-config --modversion $PKG`") is too old; version >=$MIN_VERSION is required" 1>&2
			exit 1
		fi
	fi
	echo " Found $PKG version `pkg-config --modversion $PKG`"
	eval ${VAR_PREFIX}_FOUND=1
	eval ${VAR_PREFIX}_VERSION=\"`pkg-config --modversion $PKG`\"
	eval ${VAR_PREFIX}_CFLAGS=\"`pkg-config --cflags $PKG`\"
	eval ${VAR_PREFIX}_LDFLAGS=\"`pkg-config --libs $PKG`\"
	eval ${VAR_PREFIX}_INCDIR=\"`pkg-config --variable=includedir $PKG`\"
	eval ${VAR_PREFIX}_LIBDIR=\"`pkg-config --variable=libdir $PKG`\"
}

find_hdf5(){
	PKG=hdf5
	echo "Looking for $PKG..."
	VAR_PREFIX=`echo $PKG | tr [:lower:] [:upper:]`
	TMP_FOUND=`eval echo "$"${VAR_PREFIX}_FOUND`
	if [ "$TMP_FOUND" ]; then return; fi

	which h5cc 2>&1 > /dev/null
	if [ "$?" -ne 0 ]; then return; fi

	which h5ls 2>&1 > /dev/null
	if [ "$?" -eq 0 ]; then
		HDF5_VERSION=`h5ls --version | sed 's/.* \([0-9.]*\)/\1/'`
		echo " Found $PKG version $HDF5_VERSION via executables in \$PATH"
		if [ $# -ge 1 ]; then
			MIN_VERSION=$1
			#TODO: actually check version
		fi
	else
		echo " h5ls not found; cannot check $PKG version"
		echo " Proceeding with unknown version and hoping for the best"
	fi
	HDF5_COMPILE_COMMAND=`h5cc -show`
	for item in $HDF5_COMPILE_COMMAND; do
		item=`echo "$item" | sed 's| |\n|g' | sed -n 's/.*-L\([^ ]*\).*/\1/p'`
		if [ -n "$item" ]; then
			POSSIBLE_HDF5_LIBDIRS="$POSSIBLE_HDF5_LIBDIRS
				$item"
		fi
	done
	for HDF5_LIBDIR in $POSSIBLE_HDF5_LIBDIRS; do
		if [ -d $HDF5_LIBDIR -a \( -e $HDF5_LIBDIR/libhdf5.a -o -e $HDF5_LIBDIR/libhdf5.so \) ]; then
			break
		fi
	done
	if [ ! -d $HDF5_LIBDIR -o ! \( -e $HDF5_LIBDIR/libhdf5.a -o -e $HDF5_LIBDIR/libhdf5.so \) ]; then
		echo " Unable to guess $PKG library directory"
		return
	fi
	POSSIBLE_HDF5_INCDIRS=`echo "$HDF5_COMPILE_COMMAND" | sed 's| |\n|g' | sed -n 's/.*-I\([^ ]*\).*/\1/p'`
	POSSIBLE_HDF5_INCDIRS="$POSSIBLE_HDF5_INCDIRS ${HDF5_LIBDIR}/../include"
	for HDF5_INCDIR in $POSSIBLE_HDF5_INCDIRS; do
		if [ -d $HDF5_INCDIR -a -e $HDF5_INCDIR/H5version.h ]; then
			break
		fi
	done
	if [ ! -d $HDF5_INCDIR -o ! $HDF5_INCDIR/H5version.h ]; then
		echo " Unable to guess $PKG include directory"
		return
	fi

	HDF5_CFLAGS="-I${HDF5_INCDIR}"
	HDF5_LDFLAGS=`echo "$HDF5_COMPILE_COMMAND" | \
	sed 's/ /\\
	/g' | \
	sed -n -E \
	-e '/^[[:space:]]*-l/p' \
	-e '/^[[:space:]]*-L/p' \
	-e '/^[[:space:]]*-Wl,/p' \
	-e 's/^[[:space:]]*.*lib([^.]*)\.a/-l\1/p' \
	-e 's/^[[:space:]]*.*lib([^.]*)\.so/-l\1/p' \
	-e 's/^[[:space:]]*.*lib([^.]*)\.dylib/-l\1/p' `
	HDF5_LDFLAGS=`echo $HDF5_LDFLAGS` # collapse to single line

	HDF5_FOUND=1
}

try_find_boost(){
	GUESS_DIR=$1
	PKG=boost
	VAR_PREFIX=`echo $PKG | tr [:lower:] [:upper:]`
	TMP_FOUND=`eval echo "$"${VAR_PREFIX}_FOUND`
	if [ "$TMP_FOUND" ]; then return; fi
	echo "Looking for $PKG in $GUESS_DIR..."
	POSSIBLE_BOOST_LIBDIRS="${GUESS_DIR}/lib ${GUESS_DIR}/lib64 ${GUESS_DIR}/lib/x86_64-linux-gnu"
	POSSIBLE_BOOST_INCDIRS="${GUESS_DIR}/include"
	for BOOST_LIBDIR in $POSSIBLE_BOOST_LIBDIRS; do
		if [ -d $BOOST_LIBDIR -a \( -e $BOOST_LIBDIR/libboost_python.a -o -e $BOOST_LIBDIR/libboost_python.so \) ]; then
			break
		fi
	done
	if [ ! -d $BOOST_LIBDIR -o ! \( -e $BOOST_LIBDIR/libboost_python.a -o -e $BOOST_LIBDIR/libboost_python.so \) ]; then
		echo " Unable to locate the boost_python libray in $GUESS_DIR"
		return
	fi
	for BOOST_INCDIR in $POSSIBLE_BOOST_INCDIRS; do
		if [ -d $BOOST_INCDIR -a -e $BOOST_INCDIR/boost/python.hpp ]; then
			break
		fi
	done
	if [ ! -d $BOOST_INCDIR -o ! $BOOST_INCDIR/boost/python.hpp ]; then
		echo " Unable to locate boost/python.hpp in $GUESS_DIR"
		return
	fi
	BOOST_CFLAGS="-I${BOOST_INCDIR}"
	BOOST_LDFLAGS="-Wl,-rpath -Wl,${BOOST_LIBDIR} -L${BOOST_LIBDIR} -lboost_python"
	BOOST_FOUND=1
	echo " Found boost in $GUESS_DIR"
}

ensure_found(){
	PKG=$1
	VAR_PREFIX=`echo $PKG | tr [:lower:] [:upper:]`
	TMP_FOUND=`eval echo "$"${VAR_PREFIX}_FOUND`
	if [ "$TMP_FOUND" ]; then return; fi
	#not found
	echo "Error: $PKG not installed or not registered with pkg-config" 1>&2
	lowername=`echo $PKG | tr [A-Z] [a-z]`
	echo "Please specify location using the --with-"$lowername" flag" 1>&2
	exit 1
}

PREFIX=/usr/local

VERSION_NUM=101001
VERSION=`echo $VERSION_NUM | awk '{
	major = int($1/100000);
	minor = int($1/100)%1000;
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
  DYN_OPT='-dynamiclib -compatibility_version $(VERSION) -current_version $(VERSION)'
fi

CC=${CC-$GUESS_CC}
CXX=${CXX-$GUESS_CXX}
AR=${AR-$GUESS_AR}
LD=${LD-$GUESS_LD}

PYTHON_EXE="python"

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
  --with-gsl=DIR           use the copy of GSL in DIR
                           assuming headers are in DIR/include
                           and libraries in DIR/lib
  --with-gsl-incdir=DIR    use the copy of GSL in DIR
  --with-gsl-libdir=DIR    use the copy of GSL in DIR
  --with-hdf5=DIR          use the copy of HDF5 in DIR
                           assuming headers are in DIR/include
                           and libraries in DIR/lib
  --with-hdf5-incdir=DIR   use the copy of HDF5 in DIR
  --with-hdf5-libdir=DIR   use the copy of HDF5 in DIR
  --with-squids=DIR        use the copy of SQuIDS in DIR
                           assuming headers are in DIR/include
                           and libraries in DIR/lib
  --with-squids-incdir=DIR        use the copy of SQuIDS in DIR
  --with-squids-libdir=DIR        use the copy of SQuIDS in DIR
For the python bindings the following flags are used:
  --with-python-bindings         enable python binding compilation
  --with-boost-incdir=DIR        use the copy of Boost in DIR
  --with-boost-libdir=DIR        use the copy of Boost in DIR
  --with-boost=DIR               use the copy of Boost in DIR
                                 assuming headers are in DIR/include
                                 and libraries in DIR/lib
  --python-bin=PYTHON_EXECUTABLE use this python executable
                                 (default is 'python')

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

	TMP=`echo "$var" | sed -n 's/^--with-gsl=\(.*\)$/\1/p'`
	if [ "$TMP" ]; then
		GSL_INCDIR="${TMP}/include";
		GSL_LIBDIR="${TMP}/lib";
	continue; fi

	TMP=`echo "$var" | sed -n 's/^--with-gsl-incdir=\(.*\)$/\1/p'`
	if [ "$TMP" ]; then GSL_INCDIR="$TMP"; continue; fi

	TMP=`echo "$var" | sed -n 's/^--with-gsl-libdir=\(.*\)$/\1/p'`
	if [ "$TMP" ]; then GSL_LIBDIR="$TMP"; continue; fi

	TMP=`echo "$var" | sed -n 's/^--with-hdf5=\(.*\)$/\1/p'`
	if [ "$TMP" ]; then
		HDF5_INCDIR="${TMP}/include";
		HDF5_LIBDIR="${TMP}/lib";
	continue; fi

	TMP=`echo "$var" | sed -n 's/^--with-hdf5-incdir=\(.*\)$/\1/p'`
	if [ "$TMP" ]; then HDF5_INCDIR="$TMP"; continue; fi

	TMP=`echo "$var" | sed -n 's/^--with-hdf5-libdir=\(.*\)$/\1/p'`
	if [ "$TMP" ]; then HDF5_LIBDIR="$TMP"; continue; fi

	TMP=`echo "$var" | sed -n 's/^--with-squids=\(.*\)$/\1/p'`
	if [ "$TMP" ]; then
		SQUIDS_INCDIR="${TMP}/include";
		SQUIDS_LIBDIR="${TMP}/lib";
	continue; fi

	TMP=`echo "$var" | sed -n 's/^--with-squids-incdir=\(.*\)$/\1/p'`
	if [ "$TMP" ]; then SQUIDS_INCDIR="$TMP"; continue; fi

	TMP=`echo "$var" | sed -n 's/^--with-squids-libdir=\(.*\)$/\1/p'`
	if [ "$TMP" ]; then SQUIDS_LIBDIR="$TMP"; continue; fi

	TMP=`echo "$var" | sed -n 's/^--with-python-bindings/true/p'`
	if [ "$TMP" ]; then PYTHON_BINDINGS=true; continue; fi

	TMP=`echo "$var" | sed -n 's/^--with-boost=\(.*\)$/\1/p'`
	if [ "$TMP" ]; then
		BOOST_INCDIR="${TMP}/include";
		BOOST_LIBDIR="${TMP}/lib";
	continue; fi

	TMP=`echo "$var" | sed -n 's/^--with-boost-libdir=\(.*\)$/\1/p'`
	if [ "$TMP" ]; then BOOST_LIBDIR="$TMP"; continue; fi
	TMP=`echo "$var" | sed -n 's/^--with-boost-incdir=\(.*\)$/\1/p'`
	if [ "$TMP" ]; then BOOST_INCDIR="$TMP"; continue; fi

	TMP=`echo "$var" | sed -n 's/^--python-bin=\(.*\)$/\1/p'`
	if [ "$TMP" ]; then PYTHON_EXE=$TMP; continue; fi

	echo "config.sh: Unknown or malformed option '$var'" 1>&2
	exit 1
done

if [ "$GSL_INCDIR" -a "$GSL_LIBDIR" ]; then
	echo "Checking manually specified GSL..."
	if [ -d "$GSL_INCDIR/gsl" \
         -a -e "$GSL_INCDIR/gsl/gsl_version.h" \
         -a -d "$GSL_LIBDIR" \
         -a -e "$GSL_LIBDIR/libgsl.a" ]; then
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
         -a -e "$HDF5_INCDIR/H5version.h" \
         -a -d "$HDF5_LIBDIR" \
         -a -e "$HDF5_LIBDIR/libhdf5.a" \
         -a -e "$HDF5_LIBDIR/libhdf5_hl.a" ]; then
		HDF5_FOUND=1
		HDF5_CFLAGS="-I$HDF5_INCDIR"
		HDF5_LDFLAGS="-L$HDF5_LIBDIR -lhdf5 -lhdf5_hl"
	else
		echo "Warning: manually specifed HDF5 not found; will attempt auto detection"
	fi
fi

#Do not use this due to broken Ubuntu package
#find_package hdf5 1.8
find_hdf5

if [ "$SQUIDS_INCDIR" -a "$SQUIDS_LIBDIR" ]; then
	echo "Checking manually specified SQUIDS..."
	if [ -d "$SQUIDS_INCDIR" \
         -a -d "$SQUIDS_LIBDIR" \
         -a -e "$SQUIDS_LIBDIR/libSQuIDS.a" ]; then
		SQUIDS_FOUND=1
		SQUIDS_CFLAGS="-I$SQUIDS_INCDIR"
		SQUIDS_LDFLAGS="-L$SQUIDS_LIBDIR -lSQuIDS"
		if $CXX --version | grep -q "Free Software Foundation"; then
			SQUIDS_CFLAGS="$SQUIDS_CFLAGS -Wno-abi"
		fi
	else
		echo "Warning: manually specifed SQUIDS not found; will attempt auto detection"
	fi
fi

find_package squids 1.2

ensure_found gsl
ensure_found hdf5
ensure_found squids

if [ ! -d ./build/ ]; then
    mkdir build;
fi
if [ ! -d ./lib/ ]; then
    mkdir lib;
fi

echo "Generating config file..."

# Somewhat evil: HDF5 does not register with pkg-config, which causes the latter
# to error out because it cannot find nuSQuIDS dependencies.
# Solution: Since we found HDF5 (hopefully correctly), register it ourselves.
echo "# WARNING: This configuration file was heutristically generated by nuSQuIDS
# and may not be complete or correct
libdir=${HDF5_LIBDIR}
includedir=${HDF5_INCDIR}" > lib/hdf5.pc
echo '
Name: HDF5
Description: "A data model, library, and file format for storing and managing data."
URL: https://www.hdfgroup.org/HDF5/' >> lib/hdf5.pc
echo "Version: ${HDF5_VERSION}" >> lib/hdf5.pc
echo "Cflags: ${HDF5_CFLAGS}
Libs: ${HDF5_LDFLAGS}
" >> lib/hdf5.pc

echo "prefix=$PREFIX" > lib/nusquids.pc
echo '
libdir=${prefix}/lib
includedir=${prefix}/inc

Name: nuSQuIDS
Description: Toolbox for neutrino oscillation experiments
URL: https://github.com/arguelles/nuSQuIDS' >> lib/nusquids.pc
echo "Version: $VERSION" >> lib/nusquids.pc
echo 'Requires: gsl >= 1.15 hdf5 >= 1.8 squids >= 1.2.0
Libs: -L${libdir} -lnuSQuIDS
Cflags: -I${includedir}
' >> lib/nusquids.pc

echo "Generating version header..."
sed -e "s|__NUSQUIDS_VERSION__|$VERSION_NUM|g" \
    -e "s|__NUSQUIDS_VERSION_STR__|$VERSION|g" \
    < resources/version.h.in > include/nuSQuIDS/version.h

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
OBJECTS = $(patsubst src/%.cpp,build/%.o,$(SOURCES))

EXAMPLES := examples/Single_energy/single_energy \
            examples/Multiple_energy/multiple_energy \
            examples/Atm_default/atm_default \
            examples/Bodies/bodies \
            examples/Xsections/xsections \
            examples/NSI/nsi \
            examples/Atm_NSI/atm_nsi \
            examples/HDF5_Write_Read/write \
            examples/HDF5_Write_Read/read \
						examples/Constant_density_layers/const_dens_layers

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

INCnuSQUIDS=$(PATH_nuSQUIDS)/include
LIBnuSQUIDS=$(PATH_nuSQUIDS)/lib

# FLAGS
CFLAGS= -O3 -fPIC -I$(INCnuSQUIDS) $(SQUIDS_CFLAGS) $(GSL_CFLAGS) $(HDF5_CFLAGS)
LDFLAGS= -Wl,-rpath -Wl,$(LIBnuSQUIDS) -L$(LIBnuSQUIDS)
LDFLAGS+= $(SQUIDS_LDFLAGS) $(GSL_LDFLAGS) $(HDF5_LDFLAGS)
EXMAPLES_FLAGS=-I$(INCnuSQUIDS) $(CXXFLAGS) $(CFLAGS)

# Project files
NAME=nuSQuIDS
STAT_PRODUCT=lib/lib$(NAME).a
DYN_PRODUCT=lib/lib$(NAME)$(DYN_SUFFIX)

# Compilation rules
all: $(STAT_PRODUCT) $(DYN_PRODUCT)

examples : $(EXAMPLES)

$(DYN_PRODUCT) : $(OBJECTS)
	@echo Linking $(DYN_PRODUCT)
	@$(CXX) $(DYN_OPT)  $(LDFLAGS) -o $(DYN_PRODUCT) $(OBJECTS)

$(STAT_PRODUCT) : $(OBJECTS)
	@echo Linking $(STAT_PRODUCT)
	@$(AR) -rcs $(STAT_PRODUCT) $(OBJECTS)

build/%.o : src/%.cpp
	@echo Compiling $< to $@
	@$(CXX) $(CXXFLAGS) -c $(CFLAGS) $< -o $@

examples/Single_energy/single_energy : $(DYN_PRODUCT) examples/Single_energy/main.cpp
	@echo Compiling single energy example
	@$(CXX) $(EXMAPLES_FLAGS) examples/Single_energy/main.cpp -lnuSQuIDS $(LDFLAGS) -o $@

examples/Multiple_energy/multiple_energy : $(DYN_PRODUCT) examples/Multiple_energy/main.cpp
	@echo Compiling multiple energy example
	@$(CXX) $(EXMAPLES_FLAGS) examples/Multiple_energy/main.cpp -lnuSQuIDS $(LDFLAGS) -o $@

examples/HDF5_Write_Read/write : $(DYN_PRODUCT) examples/HDF5_Write_Read/write.cpp
	@echo Compiling HDF5 write
	@$(CXX) $(EXMAPLES_FLAGS) examples/HDF5_Write_Read/write.cpp -lnuSQuIDS $(LDFLAGS) -o $@

examples/HDF5_Write_Read/read : $(DYN_PRODUCT) examples/HDF5_Write_Read/read.cpp
	@echo Compiling HDF5 read
	@$(CXX) $(EXMAPLES_FLAGS) examples/HDF5_Write_Read/read.cpp -lnuSQuIDS $(LDFLAGS) -o $@

examples/Atm_default/atm_default : $(DYN_PRODUCT) examples/Atm_default/main.cpp
	@echo Compiling atmospheric example
	@$(CXX) $(EXMAPLES_FLAGS) examples/Atm_default/main.cpp -lnuSQuIDS $(LDFLAGS) -o $@

build/exBody.o : examples/Bodies/exBody.h examples/Bodies/exBody.cpp
	@$(CXX) $(EXMAPLES_FLAGS) -c examples/Bodies/exBody.cpp -o $@

build/ex_bodies_main.o : examples/Bodies/exBody.h examples/Bodies/main.cpp
	@$(CXX) $(EXMAPLES_FLAGS) -c examples/Bodies/main.cpp -o $@

examples/Bodies/bodies : $(DYN_PRODUCT) build/ex_bodies_main.o build/exBody.o
	@echo Compiling bodies example
	@$(CXX) $(EXMAPLES_FLAGS) build/ex_bodies_main.o build/exBody.o -lnuSQuIDS $(LDFLAGS) -o $@

build/exCross.o : examples/Xsections/exCross.h examples/Xsections/exCross.cpp
	@$(CXX) $(EXMAPLES_FLAGS) -c examples/Xsections/exCross.cpp -o $@

examples/Xsections/xsections : $(DYN_PRODUCT) examples/Xsections/main.cpp build/exCross.o
	@echo Compiling cross section example
	@$(CXX) $(EXMAPLES_FLAGS) examples/Xsections/main.cpp build/exCross.o -lnuSQuIDS $(LDFLAGS) -o $@

examples/NSI/nsi : $(DYN_PRODUCT) examples/NSI/main.cpp
	@echo Compiling non-standard interaction example
	@$(CXX) $(EXMAPLES_FLAGS) examples/NSI/main.cpp -lnuSQuIDS $(LDFLAGS) -o $@

examples/Atm_NSI/atm_nsi : $(DYN_PRODUCT) examples/Atm_NSI/main.cpp
	@echo Compiling atmospheric non-standard interaction example
	@$(CXX) $(EXMAPLES_FLAGS) examples/Atm_NSI/main.cpp -lnuSQuIDS $(LDFLAGS) -o $@

examples/Constant_density_layers/const_dens_layers : $(DYN_PRODUCT) examples/Constant_density_layers/main.cpp
	@echo Compiling constant density layer example
	@$(CXX) $(EXMAPLES_FLAGS) examples/Constant_density_layers/main.cpp -lnuSQuIDS $(LDFLAGS) -o $@

.PHONY: clean test docs
clean:
	rm -rf $(PATH_nuSQUIDS)/build/*.o $(PATH_nuSQUIDS)/examples/*.exe $(PATH_nuSQUIDS)/$(STAT_PRODUCT) $(PATH_nuSQUIDS)/$(DYN_PRODUCT) $(PATH_nuSQUIDS)/bin/*

doxygen:
	@doxygen src/doxyfile
docs:
	@doxygen src/doxyfile

test: $(DYN_PRODUCT) $(STAT_PRODUCT)
	@cd ./test ; ./run_tests

install: $(DYN_PRODUCT) $(STAT_PRODUCT)
	@echo Installing headers in $(PREFIX)/include/nuSQuIDS
	@mkdir -p $(PREFIX)/include/nuSQuIDS
	@cp $(INCnuSQUIDS)/nuSQuIDS/*.h $(PREFIX)/include/nuSQuIDS
	@echo Installing libraries in $(PREFIX)/lib
	@mkdir -p $(PREFIX)/lib
	@cp $(DYN_PRODUCT) $(STAT_PRODUCT) $(PREFIX)/lib
	@echo Installing config information in $(PREFIX)/lib/pkgconfig
	@mkdir -p $(PREFIX)/lib/pkgconfig
	@cp lib/nusquids.pc $(PREFIX)/lib/pkgconfig' >> ./Makefile
if `pkg-config hdf5`; then
:
else
	echo '	@cp lib/hdf5.pc $(PREFIX)/lib/pkgconfig' >> ./Makefile
fi

echo "
export CXX=\"${CXX}\"
export CFLAGS=\"${CFLAGS} ${SQUIDS_CFLAGS} ${GSL_CFLAGS} ${HDF5_CFLAGS}\"
export CXXFLAGS=\"${CXXFLAGS} -std=c++11\"
export LDFLAGS=\"${LDFLAGS} ${SQUIDS_LDFLAGS} ${GSL_LDFLAGS} ${HDF5_LDFLAGS} -lnuSQuIDS\"
" > test/env_vars.sh
if uname | grep -q 'Darwin' ; then
	printf "export DYLD_LIBRARY_PATH=\"" >> test/env_vars.sh
	if [ "$DYLD_LIBRARY_PATH" ]; then
		printf "${DYLD_LIBRARY_PATH}:" >> test/env_vars.sh
	fi
	printf "/lib:/usr/lib:\${HOME}/lib:/usr/local/lib:" >> test/env_vars.sh
else
	printf "export LD_LIBRARY_PATH=\"" >> test/env_vars.sh
fi
echo "../lib:${SQUIDS_LIBDIR}:${GSL_LIBDIR}:${HDF5_LIBDIR}\"" >> test/env_vars.sh

if [ $PYTHON_BINDINGS ]; then
	echo "Generating Python bindings makefile..."

	if [ "$BOOST_INCDIR" -a "$BOOST_LIBDIR" ]; then
		echo "Checking manually specified boost..."
		if [ -d "$BOOST_INCDIR" \
		  -a -d "$BOOST_LIBDIR" \
		  -a -e "$BOOST_INCDIR/boost/python.hpp" \
		  -a -e "$BOOST_LIBDIR/libboost_python.a" ]; then
			BOOST_FOUND=1
			BOOST_CFLAGS="-I$BOOST_INCDIR"
			BOOST_LDFLAGS="-Wl,-rpath -Wl,${BOOST_LIBDIR} -L$BOOST_LIBDIR -lboost_python"
		else
			echo "Warning: manually specifed boost not found; will attempt auto detection"
		fi
	fi

	try_find_boost /usr
	try_find_boost /usr/local

	if [ -z $BOOST_LIBDIR -a  -z $BOOST_INCDIR ]; then
		echo "Error: Specify BOOST library path using --with-boost-libdir and BOOST include path using --with-boost-incdir."
		exit 1
	fi
	if [ -z $BOOST_LIBDIR ]; then
		echo "Error: Specify BOOST library path using  --with-boost-libdir."
		exit 1
	fi
	if [ -z $BOOST_INCDIR ]; then
		echo "Error: Specify BOOST include path using  --with-boost-incdir."
		exit 1
	fi

	PYTHONVERSION=`${PYTHON_EXE} -c 'import sys; print(str(sys.version_info.major)+"."+str(sys.version_info.minor))'`
	PYTHONLIBPATH=`${PYTHON_EXE} -c 'import sys; import re; print([ y for y in sys.path if re.search("\/lib\/python'$PYTHONVERSION'$",y)!=None ][0]);'`
	PYTHONINCPATH=`find $PYTHONLIBPATH/../../include/ -name "python${PYTHONVERSION}*"`
	PYTHONNUMPYCHECK=`${PYTHON_EXE} ./resources/check_numpy.py`
	if [ -z "$PYTHONNUMPYCHECK" ]; then
		echo "Error: numpy not found. Specify numpy installation location with --with-numpy-incdir."
		exit 1
	fi
	PYTHONNUMPYINC=$PYTHONNUMPYCHECK
	echo "
# Compiler
CC=$CC
CXX=$CXX
AR=$AR
LD=$LD

DYN_SUFFIX=.so
DYN_OPT=$DYN_OPT

PYTHON_VERSION = ${PYTHONVERSION}
LDFLAGS+= -L${PYTHONLIBPATH}  -L${PYTHONLIBPATH}/..
LDFLAGS+= -lpython${PYTHONVERSION}
LDFLAGS+= ${BOOST_LDFLAGS}
LDFLAGS+= ${SQUIDS_LDFLAGS} ${GSL_LDFLAGS} ${HDF5_LDFLAGS}
INCCFLAGS+= -I${PYTHONINCPATH} ${SQUIDS_CFLAGS} ${GSL_CFLAGS} ${HDF5_CFLAGS} -I../include/
INCCFLAGS+= ${BOOST_CFLAGS}
INCCFLAGS+= -I${PYTHONNUMPYINC}
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
STAT_PRODUCT=$(PATH_nuSQUIDSpy)/bindings/$(NAME).a
DYN_PRODUCT=$(PATH_nuSQUIDSpy)/bindings/$(NAME)$(DYN_SUFFIX)

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

$(STAT_PRODUCT) : $(OBJECTS)
	@echo Linking $(STAT_PRODUCT)
	@$(AR) -rcs $(STAT_PRODUCT) $(OBJECTS)

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
' > include/nuSQuIDS/global.h
echo "
#define XSECTION_LOCATION \"$nusqpath/data/xsections/\"
#define SUN_MODEL_LOCATION  \"$nusqpath/data/astro/bs05_agsop.dat\"
#define SUN_MODEL_NELECTRON_LOCATION \"$nusqpath/data/astro/nele_bs05op.dat\"
#define EARTH_MODEL_LOCATION \"$nusqpath/data/astro/EARTH_MODEL_PREM.dat\"
" >> include/nuSQuIDS/global.h
echo '
#endif
' >> include/nuSQuIDS/global.h

echo "Done."
echo "To build library, run the following: make
After, to build examples: make examples"
