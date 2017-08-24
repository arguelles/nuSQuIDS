LIBS = `pkg-config --cflags nusquids` `pkg-config --libs nusquids`

CXX= g++ -std=c++11
all: main.o exe

main.o: main.cpp
	$(CXX) -c $*.cpp
exe: $(OBJFILES)
	$(CXX)  ./main.o $(LIBS) -o atm_benchmark.exe 
clean:
	rm -f *.o *.mod *.exe









