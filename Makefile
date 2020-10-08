CXX=g++-8
NVCC=nvcc
CXXFLAGS=-I include/ -I /usr/local/cuda/include -Wall -Wextra -O3 -std=c++17 -fPIC -fopenmp
LDFLAGS=-lblas -lgsl -lgslcblas -lm

lib: algorithms.o grid.o ising.o #src/model/potts/potts.o src/model/xy/xy.o src/model/hb/hb.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) src/model/ising/ising.o src/grid/grid.o src/algorithm.o -shared -fPIC -o lib/libheisenberg.so
	$(CXX) $(CXXFLAGS) $(LDFLAGS) src/model/ising/ising.o src/grid/grid.o src/algorithm.o -fPIC -c -o src/libheisenberg.o
	ar crv src/libheisenberg.o lib/libheisenberg.a
	ranlib lib/libheisenberg.a
algorithms.o:
	$(CXX) $(CXXFLAGS) src/algorithm.cpp $(LDFLAGS) -c -o src/algorithm.o
grid.o:
	$(CXX) $(CXXFLAGS) src/grid/grid.cpp $(LDFLAGS) -c -o src/grid/grid.o
ising.o:
	$(CXX) $(CXXFLAGS) src/model/ising/ising.cpp $(LDFLAGS) -c -o src/model/ising/ising.o

test: main.cpp lib/libheisenberg.a
	$(CXX) $(CXXFLAGS) main.cpp lib/libheisenberg.a $(LDFLAGS) -o bin/test
doc: include/algorithm.h include/grid/grid.h include/model/model.h include/model/ising/ising.h
	doxygen Doxyfile

clean:
	rm -rf src/*.o
realclean:
	make clean
	rm -f lib/*
