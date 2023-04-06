all: main

main:
	g++ main.cpp -o main -Wno-abi -I/usr/local/include -I/usr/include/hdf5/serial -L/usr/local/lib -L/usr/lib/x86_64-linux-gnu/hdf5/serial -lnuSQuIDS -lhdf5 -lSQuIDS -lgsl -lgslcblas -lm -g -O0 -pthread -lhdf5_serial_hl

clean:
	rm -f main
