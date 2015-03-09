# Compiler
CXX=g++
AR=ar
LD=ld

PYTHON_VERSION = 2.7

PATH_nuSQUIDS= $(shell pwd)/../../../
PATH_nuSQUIDSpy= $(shell pwd)/../

SOURCES = $(wildcard *.cpp)
OBJECTS = $(SOURCES:.cpp=.o)

# Directories
LIBDIR=$(PATH_nuSQUIDS)/lib
INCDIR=$(PATH_nuSQUIDS)/inc
LIBnuSQUIDS=$(PATH_nuSQUIDS)/lib
INCnuSQUIDS=$(PATH_nuSQUIDS)/inc

# FLAGS
LDFLAGS= -lgsl -lgslcblas
LDFLAGS+= -lSQuIDS
LDFLAGS+= -L$(LIBnuSQUIDS) -lnuSQuIDS
LDFLAGS+= -lpython$(PYTHON_VERSION) -lboost_python
INCCFLAGS+= -I/usr/include/python$(PYTHON_VERSION) -I../inc/
INCCFLAGS+= -I$(INCnuSQUIDS)
#INCCFLAGS+= -I$(INCHDF5)
CXXFLAGS= -g -O3 -fPIC -std=c++11 $(INCCFLAGS)


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
