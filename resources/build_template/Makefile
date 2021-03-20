# Set this to the name you want for your compiled program
PROGRAM:=

# Set this variable to a space-separated list of all C++ 
# implementation files for your program, e.g.:
# IMPLEMENTATION_FILES=main.cpp functions.cpp
IMPLEMENTATION_FILES:=

# Put in this variable the names of any header files you write 
# which are used by your implementation files, e.g.:
# HEADERS:=functions.h
HEADERS:=

# c++11 is the minimum version which works with nuSQuIDS.
# However, if you have a newer compiler and wish to use newer
# language features you can change this to more recent versions, 
# e.g. c++14, c++17, etc.
CXX_VERSION:=c++11

# This is the file extension assumed for implementation files.
# E.g. if you name your files main.C, foo.C, etc., set this to C
CXX_SUFFIX:=cpp

# Add here any additional flags you want to use when compiling,
# e.g. -O3 for maximum optimization, or -g to include debug info
CXXFLAGS+=

# Add here any additional flags you want to use when linking,
# e.g. -lfftw3 to link against the 'FFTW' library
LDFLAGS+=




################################################################
#   DON"T CHANGE THE PARTS BELOW (UNLESS YOU REALLY NEED TO)   #
################################################################
CXX?=c++

OBJECTS=${IMPLEMENTATION_FILES:.${CXX_SUFFIX}=.o}

CXXFLAGS+=-std=${CXX_VERSION} $$(cat .nusquids_cxxflags)
LDFLAGS+=$$(cat .nusquids_ldflags)

all : ${PROGRAM}
	@( [ "${PROGRAM}" ] || echo "Error: You have not set a program name in the Makefile") && [ "${PROGRAM}" ]

${PROGRAM} : ${OBJECTS} Makefile .nusquids_ldflags
	@( [ "${OBJECTS}" ] || echo "Error: You have not defined any implementation files in the Makefile") && [ "${OBJECTS}" ]
	@( [ "${CXX_SUFFIX}" ] || echo "Error: CXX_SUFFIX is not set in the Makefile") && [ "${CXX_SUFFIX}" ]
	@( echo "${OBJECTS}" | (! grep -q '\.[^o]') || echo "Error: At least one file in your IMPLEMENTATION_FILES does not have the suffix specified by CXX_SUFFIX") && (echo "${OBJECTS}" | (! grep -q '\.[^o]'))
	${CXX} ${OBJECTS} ${LDFLAGS} -o ${PROGRAM}

# This pattern rule works for GNU make
%.o : %.${CXX_SUFFIX} Makefile ${HEADERS} .nusquids_cxxflags
	@( [ "${CXX_VERSION}" ] || echo "Error: CXX_VERSION is not set in the Makefile") && [ "${CXX_VERSION}" ]
	${CXX} -c ${CXXFLAGS} $< -o $@

# This old-fashioned suffix rule works for BSD make
.SUFFIXES : .${CXX_SUFFIX} .o
.${CXX_SUFFIX}.o : Makefile ${HEADERS} .nusquids_cxxflags
	@( [ "${CXX_VERSION}" ] || echo "Error: CXX_VERSION is not set in the Makefile") && [ "${CXX_VERSION}" ]
	${CXX} -c ${CXXFLAGS} $< -o $@

.nusquids_cxxflags : Makefile
	pkg-config nusquids squids --cflags > .nusquids_cxxflags

.nusquids_ldflags : Makefile
	pkg-config nusquids squids --libs > .nusquids_ldflags

.PHONY: all clean %.o
clean : 
	@( echo "${OBJECTS}" | (! grep -q '\.[^o]') || echo "Error: At least one file in your IMPLEMENTATION_FILES does not have the suffix specified by CXX_SUFFIX") && (echo "${OBJECTS}" | (! grep -q '\.[^o]'))
	rm -rf ${PROGRAM} ${OBJECTS} .nusquids_cxxflags .nusquids_ldflags