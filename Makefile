
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
	$(CXX) $(CXXFLAGS) $(LIBDIR) -L$(LIBSQUIDS) -L$(LIBnuSQUIDS) $< -o $@ -lSQUIDS -lnuSQUIDS -lgsl -lgslcblas

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

