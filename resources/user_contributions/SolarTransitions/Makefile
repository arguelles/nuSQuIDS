# FLAGS
CFLAGS= -O3 -fPIC -std=c++11
CFLAGS+= -I./include `pkg-config --cflags squids nusquids hdf5`
LDFLAGS+= `pkg-config --libs squids nusquids hdf5` -lhdf5_hl -lpthread

all: solar_transitions

build/SolarModel.o : SolarModel.h SolarModel.cpp
	@$(CXX) $(CFLAGS) -c SolarModel.cpp -o $@

build/SolarProbabilities.o : SolarProbabilities.h SolarProbabilities.cpp
	@$(CXX) $(CFLAGS) -c SolarProbabilities.cpp -o $@

solar_transitions: main.cpp build/SolarModel.o build/SolarProbabilities.o
	@echo Compiling solar transitions calculation
	@$(CXX) $(CFLAGS) main.cpp build/SolarModel.o build/SolarProbabilities.o -lnuSQuIDS $(LDFLAGS) -o $@

.PHONY: clean
clean:
	rm -rf ./build/*.o solar_transitions
