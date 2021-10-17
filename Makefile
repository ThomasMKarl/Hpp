CXX=g++
NVCC=nvcc
CXXFLAGS=-Wall -Wextra -O3 -std=c++17 -fPIC -fopenmp

INCLUDE_DIR=include
INCLUDE_CUDA_DIR=/usr/local/cuda/include
INCLUDE_QT_DIR=/usr/include/qt
INCLUDE_TEST_DIR=test/include

LIB_DIR=lib
LIB_CUDA_DIR=/usr/local/cuda/lib64

LDFLAGS=-lhdf5
LDFLAGS_TEST_DIR=-lgtest

all:
	make lib
	make test
	make doc
	make examples
	make qt

lib: simulation.o grid.o ising.o
	$(CXX) $(CXXFLAGS) -I $(INCLUDE_DIR) -I $(INCLUDE_CUDA_DIR) -L $(LIB_CUDA_DIR) $(LDFLAGS) src/model/ising/ising.o src/grid/grid.o src/simulation.o -shared -fPIC -o lib/libheisenberg.so
	$(CXX) $(CXXFLAGS) -I $(INCLUDE_DIR) -I $(INCLUDE_CUDA_DIR) -L $(LIB_CUDA_DIR) $(LDFLAGS) src/model/ising/ising.o src/grid/grid.o src/simulation.o -fPIC -c --static -o lib/libheisenberg.a
	ar crv src/libheisenberg.a
	mv src/libheisenberg.a lib/libheisenberg.a
	ranlib lib/libheisenberg.a

simulation.o:
	$(CXX) $(CXXFLAGS) -I $(INCLUDE_DIR) -I $(INCLUDE_CUDA_DIR) -I $(INCLUDE_QT_DIR) -L $(LIB_CUDA_DIR) $(LDFLAGS) src/simulation.cpp -c -o src/simulation.o
grid.o:
	$(CXX) $(CXXFLAGS) -I $(INCLUDE_DIR) -I $(INCLUDE_CUDA_DIR) -L $(LIB_CUDA_DIR) $(LDFLAGS) src/grid/grid.cpp -c -o src/grid/grid.o
ising.o:
	$(CXX) $(CXXFLAGS) -I $(INCLUDE_DIR) -I $(INCLUDE_CUDA_DIR) -L $(LIB_CUDA_DIR) $(LDFLAGS) src/model/ising/ising.cpp -c -o src/model/ising/ising.o

test: test/main.cpp include/grid/grid_test.h include/model/ising/ising_test.h
	$(CXX) $(CXXFLAGS) -I $(INCLUDE_DIR) -I $(INCLUDE_TEST_DIR) -I $(INCLUDE_CUDA_DIR) -I $(INCLUDE_QT_DIR) -L $(LIB_CUDA_DIR) -L $(LIB_DIR) $(LDFLAGS) main.cpp lib/libheisenberg.a -o bin/test
	export GTEST_OUTPUT="xml:./test/test_log.xml"
	./bin/test

doc: include/simulation.h include/grid/grid.h include/model/model.h include/model/ising/ising.h
	doxygen Doxyfile

examples: lib
	$(CXX) $(CXXFLAGS) -I $(INCLUDE_DIR) -I $(INCLUDE_QT_DIR) -I $(INCLUDE_CUDA_DIR) $(LDFLAGS) src/simulateAll3d.cpp -o tools/examples/simulateAll3d lib/libheisenberg.a

qt:
	qmake -o src/qt_application/Makefile src/qt_application/
	make -C src/qt_application/


clean:
	rm -rf src/*.o
realclean:
	make clean
	rm -f lib/*
