# Template for Building Simple C++ Programs

Since nuSQuIDS is a library, to do calculations with it one must generally write a program linked against it. 
One option is to use the python bindings, which allow use of the library either through a python script, or an interactive python interpreter session. 
However, not all features of the library can be efficiently exposed to python, and the python interface has overheads. 
It is often useful to program against the library's native C++ API, but preparing a C++ program for compilation is a task many people find tedious, so this template is provided for greater convenience. 

## Structure of the Template

The template itself is the single `Makefile` in this directory. 
It has two parts:
First, a set of variables with descriptions which the user should set to control what will be built and how. 
Second, a set of variables and rules which form the implementation for compiling the user's program and ensuring that it is suitably linked with nuSQuIDS. 
The template rules assume that the nuSQuIDS library has been compiled and installed so that the included `pkg-config` definition is available. 
Note that, depending on your system and where you chose to install the library, this may require adding the install location to your `PKG_CONFIG_PATH` environment variable. 
The two halves of the makefile are separated by a comment which advises against modifying the second half, as doing so should not be necessary in most circumstances. 
The first half contains comments to guide usage, while the second half is only minimally commented. 

## Example Use

Suppose that we want to write a program which will use nuSQuIDS to calculate the oscillation probability for a unit beam of muon neutrinos of a given energy through the Earth. One might write it as follows:

	#include <iostream>
	#include "nuSQuIDS/nuSQuIDS.h"
	int main(int argc, char* argv[]){
		using namespace nusquids;
		if(argc<2){
			std::cout << "Usage: earth_osc energy" << std::endl;
			return 0;
		}
		nuSQUIDS nus(3,neutrino);
		squids::Const units;
		nus.Set_E(std::stod(argv[1])*units.GeV);
		nus.Set_Body(std::make_shared<EarthAtm>());
		nus.Set_Track(std::make_shared<EarthAtm::Track>(acos(-1.0)));
		nus.Set_initial_state(marray<double,1>({3},{0,1,0}),flavor);
		nus.EvolveState();
		for(int i=0; i<3; i++)
			std::cout << nus.EvalFlavor(i) << ' ';
		std::cout << std::endl;
		return 0;
	}

This program takes an energy, assumed to be in GeV, as its single argument, and propagates a suitable beam of neutrinos through the full diameter of the Earth. 
It could be generalized and made more useful in myriad ways, but we will use this simple version for brevity. 
Let us assume that we save this file as `earth_osc.cpp`, in some new directory where we wish to compile it. 
Next, we would copy the template Makefile to the same directory:

	# In the directory containing earth_osc.cpp:
	cp ${NUSQUIDS_SOURCE_PATH}/resources/build_template/Makefile ./

Assuming that NUSQUIDS_SOURCE_PATH is the path where the nuSQuIDS source code was placed. 
Then, we edit the copied `Makefile` to set the necessary variables. 

First, we set `PROGRAM` to the name of our program, which will be `earth_osc` for this example:

	PROGRAM:=earth_osc

Next, we specify the implementation files to be compiled. We have only one, `earth_osc.cpp`:

	IMPLEMENTATION_FILES:=earth_osc.cpp

This is all that is required. It should now be sufficient to run `make` to compile the program, which we can then run. 
Doing so might look something like this:

	$ make                                                                                                  
	pkg-config nusquids squids --cflags > .nusquids_cxxflags
	/usr/local/bin/clang++ -c  -std=c++11 $(cat .nusquids_cxxflags) earth_osc.cpp -o earth_osc.o
	pkg-config nusquids squids --libs > .nusquids_ldflags
	/usr/local/bin/clang++ earth_osc.o  $(cat .nusquids_ldflags) -o earth_osc
	$ ./earth_osc 10                                                                                        
	0.0835682 0.507911 0.408521 

To clean up automatically generated files, a `clean` target is also provided which will erase them:

	$ make clean                                                                                            
	rm -rf earth_osc earth_osc.o .nusquids_cxxflags .nusquids_ldflags

## More Advanced Programs

Your program need not have only a single implementation file. 
For example, you might implement a new type of `Body`, and for neatness place it in its own header and implementation file:

	/* NewBody.h */
	#ifndef NEWBODY_H
	#define NEWBODY_H
	#include <nuSQuIDS/nuSQuIDS.h>

	class NewBody : public nusquids::Body{
		NewBody(double density);
		//Declarations of other necessary functions and data types . . . 
	};
	#endif //NEWBODY_H

	/* NewBody.cpp */
	#include "NewBody.h"
	
	NewBody::NewBody(double density){ /* implementation. . . */ }

	// Implementations of other functions. . . 

Naturally, to use this new type of `Body` you would `#include "NewBody.h"` in your main implementation file. But how to ensure that it is properly compiled and dependencies are accounted for, etc.?
First, since you have a second implementation file, add it to `IMPLEMENTATION_FILES`:

	IMPLEMENTATION_FILES:=earth_osc.cpp NewBody.cpp

Files in this variable are just separated by spaces. 
Next, add your new header file to the `HEADERS` variable:

	HEADERS:=NewBody.h

Like `IMPLEMENTATION_FILES`, this is just a space separated list, which can contain as many files as you need. 

With these changes, running `make` should compile both implementation files and link both into the program executable. 
If the header file is changed, both implementation files will be recompiled, etc. 

## Other Variables and Settings

The `CXX_VERSION` variable may be set if you want to compile your program with a different language standard. nuSQuIDS requires C++11, but should also work with C++14, C++17, C++20 and future versions. 

The `CXX_SUFFIX` variable may be changed if you prefer to name your C++ implementation files with a different file extension, such as `.cxx`, `.c++`, or `.C`. Note that no dot should be included in the value of the variable, so if you name your files lie `my_file.cxx`, just set it to `cxx`. 

`CXXFLAGS` can be set to include any additional compiler settings you want to use. 
Common additions are the `-g` option for debug information, or `-O2`/`-O3` optimization options. 

`LDFLAGS` can be likewise set to add additional options to the linker; this is most commonly used to request linking against additional libraries besides nuSQuIDS. 

## Portability and Limitations

The template makefile has been written and tested to work for both GNU make and BSD make, despite their differences in dialect. 
It also assumes the availability of some standard unix programs (the test program `[`, `cat`, `echo`, `grep`) and Bourne shell features, as well as the `pkg-config` tool. 
Naturally, a C++ compiler is required, whose name is obtained from the `CXX` environment variable, and is assumed to be `c++` if this is not set. 

The rules of the Makefile are written assuming that the compiler has a command line interface approximately like GNU `g++`. 
As such, they have been tested with `g++` and `clang++`. 
Broadly, they should also work with other Unix compilers, such as the Intel compiler (`icc`) and the Nvidia HPC C++ compiler (`nvc++`), however these have not been tested and may require additional, non-default flag settings. 

Due to a general limitation of `make`, the template makefile will not work correctly if any file name contains whitespace, or the paths to the working directory for building the program, or the install locations of any libraries, etc. contain whitespace. 

This simple build system does not encompass use of python, and as such will not build any bindings to expose custom code for use in python. 
In future, it may be extended, or a separate template may be provided to handle this more complex task. 