# Compiler
CC=clang
CXX=clang++
AR=ar
LD=clang++

DYN_SUFFIX=.dylib
DYN_OPT=-dynamiclib -install_name $(LIBnuSQUIDS)/$(DYN_PRODUCT) -compatibility_version $(VERSION) -current_version $(VERSION)

VERSION=1.0.0
PREFIX=/usr/local


PATH_nuSQUIDS=$(shell pwd)
PATH_SQUIDS=$(SQUIDS_DIR)

SOURCES = $(wildcard src/*.cpp)
OBJECTS = $(SOURCES:.cpp=.o)

EXAMPLES_SRC=$(wildcard examples/*.cpp)
EXAMPLES=$(EXAMPLES_SRC:.cpp=.exe)

CXXFLAGS= -std=c++11

# Directories

GSL_CFLAGS=-I/usr/local/Cellar/gsl/1.15/include 
GSL_LDFLAGS=-L/usr/local/Cellar/gsl/1.15/lib -lgsl -lgslcblas -lm 
HDF5_CFLAGS=-I/usr/local/Cellar/hdf5/1.8.13/include/
HDF5_LDFLAGS=-L/usr/local/Cellar/hdf5/1.8.13/lib -lhdf5 -lhdf5_hl -lhdf5_hl_cpp
SQUIDS_CFLAGS=-I/usr/local/include -I/usr/local/Cellar/gsl/1.15/include 
SQUIDS_LDFLAGS=-L/usr/local/lib -L/usr/local/Cellar/gsl/1.15/lib -lSQUIDS -lgsl -lgslcblas -lm 


INCnuSQUIDS=$(PATH_nuSQUIDS)/inc
LIBnuSQUIDS=$(PATH_nuSQUIDS)/lib

# FLAGS
CFLAGS= -O3 -fPIC -I$(INCnuSQUIDS) $(SQUIDS_CFLAGS) $(GSL_CFLAGS) $(HDF5_CFLAGS)
LDFLAGS= -Wl,-rpath -Wl,$(LIBnuSQUIDS) -L$(LIBnuSQUIDS) -lnuSQUIDS
LDFLAGS+= $(SQUIDS_LDFLAGS) $(GSL_LDFLAGS) $(HDF5_LDFLAGS)

# Project files
NAME=nuSQUIDS
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
	doxygen

test: $(DYN_PRODUCT) $(STAT_PRODUCT)
	@cd ./test ; ./run_tests

install: $(DYN_PRODUCT) $(STAT_PRODUCT)
	@echo Installing headers in $(PREFIX)/include/nuSQuIDS
	@mkdir -p $(PREFIX)/include/nuSQuIDS
	@cp $(INCDIR)/*.h $(PREFIX)/include/nuSQuIDS
	@mkdir -p $(PREFIX)/include/nuSQuIDS/detail
	@cp $(INCDIR)/detail/*.h $(PREFIX)/include/nuSQuIDS/detail
	@echo Installing libraries in $(PREFIX)/lib
	@cp $(DYN_PRODUCT) $(STAT_PRODUCT) $(PREFIX)/lib
	@echo Installing config information in $(PREFIX)/lib/pkgconfig
	@mkdir -p $(PREFIX)/lib/pkgconfig
	@cp nusquids.pc $(PREFIX)/lib/pkgconfig

